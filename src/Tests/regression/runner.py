"""Regression tier runner of the TNL-SPH unified test suite.

Executes configurations ONE case at a time (always serial, D13) through the
pinned phase pipeline::

    snapshot templates -> init -> build -> restore (finally)
        -> run -> evaluate + comp_time -> [--store-results]

Phase semantics:

* PREFLIGHT (Metis G18/H5) is FIRST in :func:`run` -- and runs before any
  ``--list`` output in the CLI, too:
  :func:`suite.runner_core.check_templates_clean` over the full
  ``examples/**/template/config.h`` set; a non-empty answer aborts with
  ``SystemExit`` naming every dirty tracked template.  This is the SIGKILL
  leftover net: a crashed tier can leave init-rewritten templates behind,
  and no subsequent run may proceed until they are restored.  (The task
  brief showed the call as ``check_templates_clean([])``; an empty list
  short-circuits to "nothing dirty" by the pinned API contract, so the
  runner passes no argument -- the full-repo scan.  QA (d) depends on it.)
* The tracked ``template/config.h`` rewritten by init is snapshotted BEFORE
  init and restored in ``finally`` after the build phase (D8) -- exceptions
  included.  run / evaluate / store edit no tracked files and therefore run
  AFTER the restore.
* init: :func:`suite.runner_core.init_case` with the per-conf
  ``init-script`` support and non-reserved kwarg forwarding.
* build (skipped with ``build=False`` / ``--no-build``):
  :func:`suite.build.build_target` -- leaf-dir rule, ``conf["target"]``
  override wins (Metis G12) -- followed by an existence assert on
  :func:`suite.build.case_binary`; a missing binary after a zero-exit build
  fails the phase with "binary missing after successful build".
* run: :func:`suite.runner_core.run_case`; a non-zero exit code fails the
  phase.
* evaluate: :func:`suite.runner_core.evaluate` (lazy
  ``regression.metrics.REGISTRY`` string resolution) into
  ``entry.tests_passed/tests_total`` -- an absent ``evaluation-function``
  yields ``(0, 0)``; ``entry.comp_time`` comes from
  :func:`suite.runner_core.parse_comp_time`.
* ``--store-results`` renames ``case_dir/results`` to
  ``results_<sanitize_tag(case-tag)>`` (port of
  ``src/tools/multirunner.py:123-128`` with the D6 sanitizer
  ``regression.metrics.references.sanitize_tag``) ONLY when the run phase
  passed AND the evaluation succeeded.  An existing destination is NEVER
  overwritten -- the entry gets a ``FileExistsError`` error instead
  (warn-overwrite protection).  Sanitized tags keep ``/``, so missing
  intermediate parents of the destination are created (the legacy code did
  not and crashed on slash-carrying tags).

Failure policy (port of ``multirunner.run_cases``'s try/except,
``multirunner.py:129-133``): any exception in any phase FAILs that phase,
sets ``entry.error`` to ``"<ClassName>: <message>"`` truncated at 240
characters, and the tier CONTINUES with the next configuration -- it never
aborts mid-tier.  Per-phase outcomes, return codes, wall seconds and error
tails are recorded in :class:`suite.results.PhaseResult`.

CLI::

    PYTHONPATH=src/Tests python3 -m regression.runner \
        --conf <module-path-or-builtin-name> [--group SUBSTR] \
        [--store-results] [--no-store] [--jobs N] [--no-build] \
        [--verbose] [--list]

``--conf`` accepts an explicit file path (absolute or CWD-relative), a path
relative to ``src/Tests/regression/`` (e.g. ``configurations/qa_smoke.py``),
or a built-in group name (default ``suite_all``).  ``--group`` is a case-tag
substring filter.

Exit codes: ``0`` = no FAIL entry; ``1`` = at least one FAIL entry; ``2`` =
harness error (configuration load failure, PREFLIGHT abort, empty ``--group``
match, unexpected exception).  ``--list`` prints the loaded entries (after
the PREFLIGHT) and exits 0 without touching anything.

Session logs: at the END of a completed tier run the CLI writes
``src/Tests/logs/log_<yyyymmdd-HHMMSS>.json`` and ``.md`` via
:func:`suite.results.write_json_log`/:func:`write_markdown_log`, with the
session document from :func:`suite.results.session_log`.  Machine info:
``{"gpu": regression.perf.detect_gpu(), "host": platform.node()}`` -- the
``regression.perf`` import is lazy and its absence is tolerated (empty gpu
string).
"""

from __future__ import annotations

import argparse
import platform
import sys
import time
from pathlib import Path

from suite import build as suite_build
from suite import config_loader, paths, reporting, results, runner_core
from suite.results import PhaseResult, ResultEntry

TIER = "regression"

REGRESSION_DIR = paths.TESTS_DIR / "regression"
LOG_DIR = paths.TESTS_DIR / "logs"

_ERROR_MAX = 240
_LOG_TAIL = 2000
_TIMESTAMP_FORMAT = "%Y%m%d-%H%M%S"


def _error_text(exc: BaseException) -> str:
    """``"<ClassName>: <message>"`` truncated at ``_ERROR_MAX`` characters."""
    text = f"{type(exc).__name__}: {exc}"
    return text if len(text) <= _ERROR_MAX else text[: _ERROR_MAX - 3] + "..."


def _log_tail(exc: BaseException) -> str:
    """Keep the LAST ``_LOG_TAIL`` characters of the exception message."""
    text = str(exc)
    return text if len(text) <= _LOG_TAIL else "...\n" + text[-_LOG_TAIL:]


def _exc_returncode(exc: BaseException) -> int:
    """Return code for a phase FAIL: the causing subprocess rc when known."""
    rc = getattr(exc.__cause__, "returncode", None)
    if isinstance(rc, int) and rc != 0:
        return int(rc)
    return 1


def _preflight() -> None:
    """Metis G18/H5 net: abort naming every dirty tracked template file.

    Runs ``git diff`` over the full ``examples/**/template/config.h`` set via
    :func:`suite.runner_core.check_templates_clean` (no argument = full repo
    scan; an explicit empty list would short-circuit to "nothing dirty" and
    defeat the crash-leftover net).
    """
    dirty = runner_core.check_templates_clean()
    if dirty:
        raise SystemExit(
            "regression.runner: PREFLIGHT abort - dirty tracked template file(s)\n"
            "detected (crash leftover or manual edit; Metis H5). Restore them\n"
            "first, e.g. `git checkout -- <file>`:\n"
            + "\n".join(f"  {path}" for path in dirty)
        )
    print("preflight: tracked examples/**/template/config.h clean", flush=True)


def _filter_entries(configurations: list[dict], group: str | None) -> list[dict]:
    """Substring ``--group`` filter on ``case-tag``; aborts when nothing matches."""
    matched = [
        conf
        for conf in configurations
        if group is None or group in str(conf.get("case-tag", ""))
    ]
    if not matched:
        available = ", ".join(str(conf.get("case-tag", "")) for conf in configurations)
        raise SystemExit(
            f"regression.runner: no entries matched group filter {group!r} "
            f"(available: {available or '<none>'})"
        )
    return matched


def _store_results(entry: ResultEntry, case_dir: Path) -> None:
    """Rename ``results/`` to ``results_<sanitize_tag(case-tag)>`` (D6 sanitizer).

    Port of ``multirunner.py:123-128``.  Never overwrites: an existing
    destination raises ``FileExistsError`` naming it (warn-overwrite
    protection).  Sanitized tags keep ``/`` (e.g. ``WCSPH-DBC/...``), so
    missing intermediate parents of the destination are created first.
    """
    from regression.metrics.references import sanitize_tag  # lazy: heavy import chain

    src = case_dir / "results"
    if not src.is_dir():
        return  # the run produced nothing to store
    dest = case_dir / ("results_" + sanitize_tag(entry.tag))
    if dest.exists():
        raise FileExistsError(
            f"store-results destination already exists (refusing to overwrite): {dest}"
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    src.rename(dest)


def _run_one(
    index: int,
    total: int,
    conf: dict,
    *,
    store_results: bool,
    verbose: bool,
    do_build: bool,
) -> ResultEntry:
    """Run one configuration through the pinned phase pipeline; never raises."""
    tag = str(conf.get("case-tag", ""))
    case_rel = str(conf.get("case", ""))
    case_dir = paths.EXAMPLES_DIR / case_rel
    entry = runner_core.results_entry_for(conf, tier=TIER)

    print(f"\n[{index}/{total}] {tag}", flush=True)

    snapshot = runner_core.snapshot_tracked_templates(case_dir)
    started = time.perf_counter()
    try:
        print(f"Initializing case: {case_rel} in {case_dir}.", flush=True)
        runner_core.init_case(case_dir, conf, verbose)
        entry.phases["init"] = PhaseResult(
            status="PASS", returncode=0, seconds=round(time.perf_counter() - started, 3)
        )
        reporting.print_phase(f"init {tag} ({entry.phases['init'].seconds:.2f}s)")

        if do_build:
            print(f"Compiling case: {case_rel}.", flush=True)
            started = time.perf_counter()
            suite_build.build_target(case_rel, target=conf.get("target"))
            binary = suite_build.case_binary(case_rel, conf.get("target"))
            if not binary.is_file():
                raise RuntimeError(f"binary missing after successful build: {binary}")
            entry.phases["build"] = PhaseResult(
                status="PASS",
                returncode=0,
                seconds=round(time.perf_counter() - started, 3),
            )
            reporting.print_phase(f"build {tag} ({entry.phases['build'].seconds:.2f}s)")
    except Exception as exc:
        name = "init" if "init" not in entry.phases else "build"
        entry.phases[name] = PhaseResult(
            status="FAIL",
            returncode=_exc_returncode(exc),
            seconds=round(time.perf_counter() - started, 3),
            log=_log_tail(exc),
        )
        entry.error = _error_text(exc)
        reporting.print_phase(f"{name} {tag}", ok=False, error_detail=entry.error)
        return entry
    finally:
        runner_core.restore_tracked_templates(case_dir, snapshot)

    print(f"Executing case: {case_rel}.", flush=True)
    started = time.perf_counter()
    rc = runner_core.run_case(case_dir, verbose)
    entry.phases["run"] = PhaseResult(
        status="PASS" if rc == 0 else "FAIL",
        returncode=rc,
        seconds=round(time.perf_counter() - started, 3),
    )
    reporting.print_phase(
        f"run {tag} (exit {rc}, {entry.phases['run'].seconds:.2f}s)", ok=(rc == 0)
    )
    if rc != 0:
        entry.error = f"run.py finished with exit code {rc}"
        return entry

    try:
        entry.tests_passed, entry.tests_total = runner_core.evaluate(conf, case_dir)
        entry.comp_time = runner_core.parse_comp_time(case_dir)
        if store_results:
            _store_results(entry, case_dir)
    except Exception as exc:
        entry.error = _error_text(exc)
        reporting.print_phase(f"evaluate {tag}", ok=False, error_detail=entry.error)
    return entry


def run(
    configurations: list[dict],
    *,
    store_results: bool = False,
    verbose: bool = False,
    jobs: int | None = None,
    build: bool = True,
    group: str | None = None,
) -> tuple[int, list[ResultEntry]]:
    """Run the regression tier over ``configurations``; return (FAIL count, entries).

    PREFLIGHT first (Metis H5), then the ``--group`` substring filter, then
    strictly SERIAL per-case execution through the phase pipeline documented
    in the module docstring.  Per-case failures never abort the tier.
    """
    _preflight()
    matched = _filter_entries(configurations, group)
    if jobs is not None:
        suite_build.set_jobs(jobs)
    if build:
        suite_build.ensure_configured()
    entries = [
        _run_one(
            index,
            len(matched),
            conf,
            store_results=store_results,
            verbose=verbose,
            do_build=build,
        )
        for index, conf in enumerate(matched, 1)
    ]
    failed = sum(1 for entry in entries if results.result_failed(entry))
    return failed, entries


def _machine_info() -> dict:
    """``{"gpu": perf.detect_gpu(), "host": platform.node()}`` -- perf lazy/tolerated."""
    try:
        from regression import perf  # noqa: PLC0415 - lazy import (todo 11)
    except Exception:
        gpu = ""
    else:
        try:
            gpu = perf.detect_gpu()
        except Exception:
            gpu = ""
    return {"gpu": gpu, "host": platform.node()}


def _resolve_configurations(spec: str) -> list[dict]:
    """Resolve ``--conf``: explicit path, path under ``src/Tests/regression/``,
    or a built-in group name."""
    direct = Path(spec)
    if direct.is_file():
        return config_loader.load_configurations(direct)
    under_regression = REGRESSION_DIR / spec
    if under_regression.is_file():
        return config_loader.load_configurations(under_regression)
    return config_loader.load_builtin_group(spec)


def _print_entries(configurations: list[dict]) -> None:
    """``--list`` output: one line per loaded configuration."""
    for conf in configurations:
        case_rel = str(conf.get("case", ""))
        target = conf.get("target", f"{Path(case_rel).name}_cuda")
        print(
            f"{conf.get('case-tag', '')}"
            f"  case={case_rel}"
            f"  eval={conf.get('evaluation-function', '')}"
            f"  target={target}"
            f"  init={conf.get('init-script', 'init.py')}"
        )


def main(argv: list[str] | None = None) -> int:
    """CLI entry; see the module docstring for the exit-code contract."""
    parser = argparse.ArgumentParser(
        prog="python3 -m regression.runner",
        description="Regression tier runner of the TNL-SPH unified test suite.",
    )
    parser.add_argument(
        "--conf",
        default="suite_all",
        metavar="MODULE_OR_GROUP",
        help="configuration module path (absolute, CWD-relative, or relative to "
        "src/Tests/regression/) or a built-in group name (default: suite_all)",
    )
    parser.add_argument(
        "--group",
        default=None,
        metavar="SUBSTR",
        help="substring filter on case-tag",
    )
    parser.add_argument(
        "--store-results",
        action="store_true",
        help="rename results/ to results_<sanitized tag> after a passing run",
    )
    parser.add_argument(
        "--no-store",
        action="store_true",
        help="legacy no-op alias (results are stored ONLY with --store-results)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        metavar="N",
        help="parallel build jobs (default: 2, AGENTS hard rule)",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="skip the build phase (the case binary must already exist)",
    )
    parser.add_argument("--verbose", action="store_true", help="stream case output")
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_only",
        help="print the loaded entries (after PREFLIGHT) and exit 0",
    )
    args = parser.parse_args(argv)

    configurations = _resolve_configurations(args.conf)

    if args.list_only:
        _preflight()
        _print_entries(_filter_entries(configurations, args.group))
        return 0

    store_results = bool(args.store_results and not args.no_store)
    started = time.time()
    failed, entries = run(
        configurations,
        store_results=store_results,
        verbose=args.verbose,
        jobs=args.jobs,
        build=not args.no_build,
        group=args.group,
    )
    reporting.print_final_summary({TIER: entries})

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime(_TIMESTAMP_FORMAT)
    session = results.session_log(entries, _machine_info(), started=started)
    json_path = results.write_json_log(session, LOG_DIR / f"log_{stamp}.json")
    md_path = results.write_markdown_log(entries, LOG_DIR / f"log_{stamp}.md")
    print(f"session log (json): {json_path}")
    print(f"session log (md):   {md_path}")
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit as exc:
        if isinstance(exc.code, str):
            # harness errors (config load, PREFLIGHT, empty --group match)
            print(exc.code, file=sys.stderr)
            sys.exit(2)
        raise
    except Exception:
        import traceback

        traceback.print_exc()
        sys.exit(2)

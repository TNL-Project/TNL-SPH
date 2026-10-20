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

Unfrozen-reference gate (todo 10, run mode): after a successful run phase and
BEFORE ``evaluate`` the runner consults
``regression.metrics.references.staging_refs_only(case-tag)``.  When the tag
has a staging reference but NO committed one the entry is reported as SKIP
with reason ``"unfrozen reference (staging only)"`` and evaluation does not
happen (Metis H6: never a silent validate against an unverified staging
file; such a tag previously died in ``load_reference`` with a bare
``FileNotFoundError``).  Immediately thereafter in the same seam,
:func:`_check_committed_reference_frozen` rejects a committed reference whose
metrics carry ``atol: null`` (an unfrozen staging ref promoted without
tolerances) with an informative ``RuntimeError`` naming the file and the
offending metrics -- numpy's ``allclose`` would otherwise die on the
``None`` atol with a cryptic ``TypeError``.  Skipped entries do not store
results even with ``--store-results`` (the evaluate phase did not succeed),
and skips never affect the exit code.

Reference-staging mode (todo 10): ``--generate-references
[--from-stored]`` switches the CLI from the tier pipeline to
"collect-from-existing-artifacts".  For every resolved configuration (same
``--conf``/``--group`` resolution as run mode) the runner locates the
metrics source directory -- ``<case>/results/`` by default, or the
``results_<sanitize_tag(case-tag)>/`` directory produced by
``--store-results`` when ``--from-stored`` is given (via an ephemeral
``results`` symlink alias in a temp dir, because registry collectors append
``results/`` to the case dir themselves) -- and calls the entry's
``evaluation-function`` registry function with ``collect=True``.  Each
written reference carries ``"verified": false``, per-metric ``"atol": null``
plus a ``"tolerance_pending": true`` marker (NO invented heuristics -- the
human sets tolerances at freeze time), provenance and date, and is written
by :func:`regression.metrics.references.save_reference` -- which refuses any
target outside ``references/`` -- so this mode can only ever write inside
``references/_staging/``.  Skips per entry (printed, never failures):

- ``no evaluation to collect`` -- ``evaluation-function`` absent/None/``""``;
- ``unknown evaluation function '<name>'`` -- not in the registry;
- ``no results to collect from (run first)`` -- live ``results/`` missing;
- ``no stored results to collect from: <dir> (run --store-results first)``
  -- ``--from-stored`` source missing;
- ``collector produced no metrics (stub function?)`` -- falsy collect dict.

PREFLIGHT deliberately does not run in staging mode: the mode never rewrites
tracked templates (the H5 net guards the pipeline that does), and a sibling
tier's transient init/build window must not block staging.  ``--from-stored``
without ``--generate-references`` is a harness error (``argparse`` exit 2).
The mode ends with a MANUAL VERIFICATION CHECKLIST per written file
(eyeball the results, set atol for every metric at freeze time, remove the
pending markers, ``git mv`` out of ``_staging/``, commit) -- nothing is
auto-verified -- and a final written/skipped count.

Exit codes: ``0`` = no FAIL entry; ``1`` = at least one FAIL entry; ``2`` =
harness error (configuration load failure, PREFLIGHT abort, empty ``--group``
match, unexpected exception).  ``--list`` prints the loaded entries (after
the PREFLIGHT) and exits 0 without touching anything.  In staging mode the
exit code is ``0`` when at least one staging reference was written, ``1``
when nothing was staged, and ``2`` on harness errors.

CLI::

    PYTHONPATH=src/Tests python3 -m regression.runner \
        --conf <module-path-or-builtin-name> [--group SUBSTR] \
        [--store-results] [--no-store] [--jobs N] [--no-build] \
        [--verbose] [--list] [--generate-references [--from-stored]]

``--conf`` accepts an explicit file path (absolute or CWD-relative), a path
relative to ``src/Tests/regression/`` (e.g. ``configurations/qa_smoke.py``),
or a built-in group name (default ``suite_all``).  ``--group`` is a case-tag
substring filter.

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
import datetime
import platform
import sys
import tempfile
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


# Pinned reason for the unfrozen-reference gate (Metis H6).
_UNFROZEN_SKIP_REASON = "unfrozen reference (staging only)"


def _staging_only_reference(tag: str) -> bool:
    """True when ``tag`` has a staging reference but NO committed one."""
    from regression.metrics import references  # lazy: heavy import chain

    return references.staging_refs_only(tag) if tag else False


def _check_committed_reference_frozen(tag: str, eval_name: str | None) -> None:
    """Reject a committed reference carrying ``atol: null`` metrics.

    Such a file is an unfrozen staging ref that was promoted without setting
    tolerances (``tolerance_pending``).  The comparator never accepts a null
    atol -- ``numpy.allclose`` dies on it with a cryptic ``TypeError`` -- so
    the runner fails AT COMPARE TIME with an informative error that names the
    file and the offending metrics and states the freeze-time remedy.  The
    check is tag-scoped: only the committed reference of the entry's own
    ``case-tag`` is scanned, and only when the entry actually evaluates
    (a non-empty ``evaluation-function``).  Missing references are the
    unfrozen gate's business, not this one's.
    """
    if not tag or not eval_name:
        return
    from regression.metrics import references  # lazy: heavy import chain

    path = references.reference_path(tag)
    if not path.is_file():
        return
    ref = references.load_reference(path)
    pending = sorted(
        name
        for name, metric in ref.get("metrics", {}).items()
        if isinstance(metric, dict) and metric.get("atol") is None
    )
    if pending:
        raise RuntimeError(
            f"reference {path} is unfrozen: metrics {pending} carry atol=null"
            " -- set atol for every metric at freeze time (git mv out of"
            " references/_staging/) before committing the reference"
        )


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
        if _staging_only_reference(tag):
            # Unfrozen gate (Metis H6): SKIP, never a silent staging validate.
            entry.skip_reason = _UNFROZEN_SKIP_REASON
            print(f"SKIP {tag}: {entry.skip_reason}", flush=True)
        else:
            _check_committed_reference_frozen(tag, conf.get("evaluation-function"))
            entry.tests_passed, entry.tests_total = runner_core.evaluate(conf, case_dir)
        entry.comp_time = runner_core.parse_comp_time(case_dir)
        if store_results and not entry.skip_reason:
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


def _print_generation_checklist(
    written: list[Path], skipped: list[tuple[str, str]]
) -> None:
    """Print the skip recap and the manual-verification checklist.

    Counts are the FINAL line (pinned): ``references generated: N, skipped: M``.
    """
    for tag, reason in skipped:
        print(f"  skipped {tag}: {reason}")
    for path in written:
        print(f"  staged  {path}")
    if written:
        print(
            "\nMANUAL VERIFICATION CHECKLIST (H1/H2 freeze flow; NOTHING is"
            " auto-verified):"
        )
        for path in written:
            print(f"- {path}")
            print(
                "    1. eyeball the results_<tag>/ (or results/) outputs the"
            )
            print("       metrics were collected from")
            print(
                "    2. set atol for EVERY metric at freeze time"
                ' ("atol" is null and "tolerance_pending"'
            )
            print(
                '       was written by design -- the comparator never accepts'
                " a null atol;"
            )
            print('       remove the marker and set "verified": true')
            print(
                "    3. git mv the file from references/_staging/ to"
            )
            print("       references/ (same relative path) and commit")
        print(
            "\nUntil then the tier reports every staging-only tag as SKIP"
            f' "{_UNFROZEN_SKIP_REASON}" and never validates against it.'
        )
    print(f"references generated: {len(written)}, skipped: {len(skipped)}", flush=True)


def _record_skip(skipped: list[tuple[str, str]], tag: str, reason: str) -> None:
    """Append ``(tag, reason)`` to the skip list and print the SKIP line."""
    skipped.append((tag, reason))
    print(f"SKIP {tag}: {reason}", flush=True)


def generate_references(
    configurations: list[dict],
    *,
    group: str | None = None,
    from_stored: bool = False,
) -> int:
    """Stage unverified reference JSONs from existing result artifacts.

    Collect-from-existing-artifacts mode (todo 10): for every resolved
    configuration locate the metrics source directory, call the registry
    function in collect mode, and write the reference under
    ``references/_staging/`` (guarded by
    :func:`regression.metrics.references.save_reference`).  Returns the
    number of staging files written; skip entries are printed, never raise.
    """
    matched = _filter_entries(configurations, group)
    from regression import metrics  # lazy: heavy import chain (todo 8)
    from regression.metrics import references  # lazy: heavy import chain

    registry = metrics.REGISTRY
    total = len(matched)
    origin = "stored" if from_stored else "live"
    print(
        f"reference staging: collecting metrics from {origin} results"
        f" ({total} entr{'y' if total == 1 else 'ies'})",
        flush=True,
    )
    written: list[Path] = []
    skipped: list[tuple[str, str]] = []
    with tempfile.TemporaryDirectory(prefix="tnl-sph-refgen-") as alias_root:
        for index, conf in enumerate(matched, 1):
            tag = str(conf.get("case-tag", ""))
            print(f"\n[{index}/{total}] {tag}", flush=True)
            eval_name = conf.get("evaluation-function")
            if not eval_name:
                _record_skip(skipped, tag, "no evaluation to collect")
                continue
            fn = registry.get(eval_name) if isinstance(eval_name, str) else None
            if fn is None:
                _record_skip(skipped, tag, f"unknown evaluation function '{eval_name}'")
                continue
            case_dir = paths.EXAMPLES_DIR / str(conf.get("case", ""))
            if from_stored:
                stored_dir = case_dir / ("results_" + references.sanitize_tag(tag))
                if not stored_dir.is_dir():
                    _record_skip(
                        skipped,
                        tag,
                        f"no stored results to collect from: {stored_dir}"
                        " (run --store-results first)",
                    )
                    continue
                # Collectors read <case_dir>/results/: alias the stored
                # dir as "results" of an ephemeral scratch dir.
                collect_dir = Path(alias_root) / f"case_{index}"
                collect_dir.mkdir()
                (collect_dir / "results").symlink_to(stored_dir)
            elif (case_dir / "results").is_dir():
                collect_dir = case_dir
            else:
                _record_skip(skipped, tag, "no results to collect from (run first)")
                continue
            collected = fn(collect_dir, collect=True)
            if not collected:
                _record_skip(skipped, tag, "collector produced no metrics (stub function?)")
                continue
            ref = {
                "case-tag": tag,
                "source": (
                    f"generated from {origin} results on {platform.node()}"
                    " by --generate-references"
                ),
                "date": datetime.date.today().isoformat(),
                "verified": False,
                "metrics": {
                    name: {
                        "value": value,
                        "atol": None,  # set at freeze time; never invented
                        "rtol": 0.0,
                        "tolerance_pending": True,
                    }
                    for name, value in collected.items()
                },
            }
            path = references.save_reference(
                ref, references.reference_path(tag, staging=True)
            )
            written.append(path)
            print(f"staged {path}", flush=True)
    print("")
    _print_generation_checklist(written, skipped)
    return len(written)


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
    parser.add_argument(
        "--generate-references",
        action="store_true",
        help="reference-staging mode: collect metrics from existing results "
        "and write unverified references under references/_staging/ "
        "(0 = at least one staged, 1 = nothing staged, 2 = harness error)",
    )
    parser.add_argument(
        "--from-stored",
        action="store_true",
        help="with --generate-references: collect from the "
        "results_<sanitized tag>/ directories produced by --store-results "
        "instead of results/",
    )
    args = parser.parse_args(argv)

    if args.from_stored and not args.generate_references:
        parser.error("--from-stored only makes sense with --generate-references")

    configurations = _resolve_configurations(args.conf)

    if args.list_only:
        _preflight()
        _print_entries(_filter_entries(configurations, args.group))
        return 0

    if args.generate_references:
        written = generate_references(
            configurations, group=args.group, from_stored=args.from_stored
        )
        return 0 if written else 1

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

"""Orchestrator CLI of the TNL-SPH unified test suite (todo 27).

Single entry point ``main(argv) -> int`` dispatched from the thin
``src/Tests/runTests.py`` shebang script.  Dispatch is strictly SERIAL: the
selected tiers are executed one after another in one process in the pinned
order ``unit -> infrastructure -> regression`` (D13/G18 -- no concurrency
inside this module, ever).

Dispatch map (lazy imports per tier -- import-cost isolation, and
regression-only runs never touch unit/fixture modules):

    unit            -> unit.runner.run(groups, verbose)
    infrastructure  -> see the todo-21 seam below
    regression      -> suite.config_loader resolution of --conf, then
                       regression.runner.run(...)

Todo-21 infrastructure seam (pinned): the boot/distributed driver lands as
``src/Tests/infrastructure/infra_runner.py`` exposing
``cli_run(groups=<list|None>, verbose=<bool>) -> (fail_count, entries)``.
This module imports it LAZILY inside the infrastructure dispatch: when the
module exists and exposes a callable ``cli_run`` the whole tier is delegated
to it; otherwise the todo-20 fixtures runner alone drives the tier
(``infrastructure.fixture_runner.run``) and a notice is printed.  When even
the fixtures runner cannot be imported, a SystemExit with a clear
"infrastructure boot driver unavailable" message aborts with exit code 2
(harness error) -- but only when the infrastructure tier was actually
requested (the laziness is exactly what keeps unit-only runs working).

Entry normalization: regression entries are native
:class:`suite.results.ResultEntry` objects and pass through untouched;
unit/fixtures entries are the runners' local plain dicts
``{tag, group, tier, status, seconds, error}`` and are wrapped -- never
reshaped -- by :func:`_to_result_entry` (status becomes the ``run`` phase,
SKIP error text becomes ``skip_reason``), keeping the session JSON uniform.

Verdict rule (todo-27 pin, Metis H2): an entry fails the session when
:func:`suite.results.result_failed` is true OR when its test counts mismatch
(``tests_total > 0 and tests_passed != tests_total``).  The mismatch case is
deliberately a CLI-level verdict on top of the raw entry: todo 9 keeps a
tests-mismatched entry phase-PASS (red Tests cell, ``tests_*`` recorded in
the session log) and the gate lives here -- today this makes the pinned
qa_smoke run exit 1 because the one gated DBC reference is stale on this
branch (1/4), until the user freezes references at H0/H2.

Exit codes (pinned):
    0 - every entry passed (fixture/test SKIPs do NOT fail the session),
    1 - at least one FAIL entry or a tests mismatch,
    2 - harness error (bad flag/group/conf, missing build product under
        --no-build, import plumbing, configure failure, unexpected
        exception; argparse own usage errors also exit 2).

``--jobs N`` flows into :func:`suite.build.set_jobs` and NOWHERE else (single
ownership; ``regression.runner.run`` is then called with ``jobs=None`` so it
never overrides the module default).  ``--build`` (default) configures first
via :func:`suite.build.ensure_configured` when a build-requiring tier
(regression) is selected; ``--no-build`` skips both the configure AND target
builds (``build=False`` runner arg) and requires prebuilt case binaries -- a
missing binary in that mode is a harness error (exit 2, never 1), since no
phase can rebuild it.

``regression.runner``'s own ``__main__`` may grow ``--generate-references``
(todo 10) -- it is intentionally NOT an orchestrator flag today; reconciliation
is a todo-28 concern (no flag drift from parallel work).

Console output goes through :mod:`suite.reporting` (the shared rich console);
harness-error diagnostics print to stderr mirroring ``regression.runner``'s
own main, and tracebacks surface only at ``--verbose``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import sys
import time
import traceback
from pathlib import Path
from typing import Callable, cast

# See the module docstring for the pinned todo-21 seam this alias encodes.
_Driver = Callable[..., "tuple[int, list]"]

# Self-locate: put src/Tests on sys.path for direct execution of this file
# (`python3 src/Tests/suite/cli.py`); a no-op when runTests.py already did it.
_TESTS_DIR = Path(__file__).resolve().parent.parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))

import suite.bootstrap  # noqa: E402,F401  (single allowed sys.path hook)

from suite import paths, reporting, results  # noqa: E402
from suite import build as suite_build  # noqa: E402
from suite.results import PhaseResult, ResultEntry  # noqa: E402

SUITE_VERSION = "1.0"
SESSION_SCHEMA = 1  # pinned by suite.results.session_log

# Pinned canonical tier order (also the serial dispatch order).
TIERS = ("unit", "infrastructure", "regression")

DEFAULT_CONF = "suite_all"
DEFAULT_LOG_DIR = paths.TESTS_DIR / "logs"
# Session log name: log_YYYY-MM-DD_HH-MM-SS.{json,md} (todo 27 pin).
_LOG_STAMP_FORMAT = "%Y-%m-%d_%H-%M-%S"


# ── console helpers (all console output goes through suite.reporting) ──────


def _say(message: str) -> None:
    """Print one line through the shared rich console (markup allowed)."""
    reporting._console().print(message)


def _say_safe(text: str) -> str:
    """Escape a payload for embedding inside a markup-carrying ``_say`` line."""
    return reporting._escape(text)


# ── argument parsing ───────────────────────────────────────────────────────


def _parse_args(argv: "list[str] | None") -> argparse.Namespace:
    """Build the parser and parse ``argv`` (argparse itself exits 2 on usage errors)."""
    parser = argparse.ArgumentParser(
        prog="runTests.py",
        description="Unified TNL-SPH test-suite orchestrator: unit, "
        "infrastructure and regression tiers behind one command, one "
        "summary and one session log.  Tiers execute serially.",
        epilog="""\
Per-tier --group semantics (tier-native, by design):
  unit            exact check-group names (unit/ module stems, e.g. decomposition)
  infrastructure  substring match over fixture names (boot/distributed land with todo 21)
  regression      substring match over case-tag within the resolved --conf

examples:
  %(prog)s --tier unit --group decomposition --group vtk_io
  %(prog)s --tier regression --conf configurations/qa_smoke.py
  %(prog)s --tier all --no-build --jobs 2
  %(prog)s --report src/Tests/logs/log_2026-10-20_12-00-00.json
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--tier",
        action="append",
        choices=("all",) + TIERS,
        default=None,
        metavar="TIER",
        help="tier to run: all, unit, infrastructure or regression; may be "
        "given multiple times ('all' expands to the three, default: all)",
    )
    parser.add_argument(
        "--group",
        dest="groups",
        action="append",
        default=None,
        metavar="NAME",
        help="restrict the run to a group within each selected tier "
        "(repeatable; tier-native semantics, see below)",
    )
    parser.add_argument(
        "--conf",
        default=None,
        metavar="PATH_OR_GROUP",
        help="regression configurations: a python file path (absolute, "
        "CWD-relative or relative to src/Tests/regression/) or a built-in "
        f"group name from src/Tests/regression/configurations/ "
        f"(default: {DEFAULT_CONF})",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=2,
        metavar="N",
        help="parallel build jobs (default: 2, AGENTS hard rule; flows into "
        "suite.build.set_jobs only)",
    )
    parser.add_argument(
        "--build",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="build required targets before running (default), configuring "
        "the build tree first when needed; --no-build skips configure AND "
        "target builds and requires prebuilt case binaries",
    )
    parser.add_argument(
        "--build-type",
        choices=("Release", "Debug", "RelWithDebInfo"),
        default="Release",
        help="CMAKE_BUILD_TYPE used only when a fresh configure is needed "
        "(default: Release; an existing cache is never reconfigured)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list per-tier crumbs (unit groups with check counts, fixture "
        "names, regression groups or --conf entries) without executing",
    )
    parser.add_argument(
        "--report",
        default=None,
        metavar="LOG_JSON",
        help="reprint mode: load a session JSON log and render its summary "
        "with its own recorded machine line (ignores tier selection flags)",
    )
    parser.add_argument(
        "--store-results",
        action="store_true",
        help="regression: rename results/ to results_<tag>/ after a passing "
        "run (never overwrites an existing backup)",
    )
    parser.add_argument(
        "--log-dir",
        default=str(DEFAULT_LOG_DIR),
        metavar="DIR",
        help=f"directory for session logs (default: {DEFAULT_LOG_DIR})",
    )
    parser.add_argument(
        "--gpu",
        default=None,
        metavar="NAME",
        help="GPU name override for the informational perf report (default: "
        "autodetect via nvidia-smi)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="stream per-tier detail and print full tracebacks on harness errors",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="print the suite and session-log schema versions and exit 0",
    )
    return parser.parse_args(argv)


def _normalize_tiers(raw: "list[str] | None") -> list[str]:
    """Expand ``all`` and dedupe, preserving the canonical tier order."""
    requested = set(raw or ["all"])
    if "all" in requested:
        return list(TIERS)
    return [tier for tier in TIERS if tier in requested]


# ── entry normalization ────────────────────────────────────────────────────


def _to_result_entry(raw, fallback_tier: str) -> ResultEntry:
    """Wrap one runner entry into :class:`ResultEntry` (wrap-only, no reshape).

    Native ``ResultEntry`` (regression) passes through untouched.  Plain dicts
    (unit: PASS/FAIL, fixtures: PASS/FAIL/SKIP) map their single ``status``
    onto the ``run`` phase; a SKIP entry's ``error`` text is a reason, not a
    failure, so it lands in ``skip_reason`` -- this is exactly what keeps
    fixture SKIPs out of the failure count (exit-code rule).
    """
    if isinstance(raw, ResultEntry):
        return raw
    status = str(raw.get("status", "PASS"))
    seconds = float(raw.get("seconds", 0.0) or 0.0)
    error = str(raw.get("error") or "")
    entry = ResultEntry(
        tag=str(raw.get("tag", "")),
        group=str(raw.get("group", "")),
        tier=str(raw.get("tier", "") or fallback_tier),
        phases={
            "run": PhaseResult(
                status=status,
                returncode=0 if status != "FAIL" else 1,
                seconds=seconds,
            )
        },
        comp_time=seconds,
    )
    if status == "SKIP":
        entry.skip_reason = error
    elif error:
        entry.error = error
    return entry


def _tests_mismatch(entry: ResultEntry) -> bool:
    """Pinned verdict extension (Metis H2): counted as FAIL at the CLI gate."""
    return entry.tests_total > 0 and entry.tests_passed != entry.tests_total


# ── per-tier dispatch (lazy imports, serial execution) ─────────────────────


def _dispatch_unit(groups: "list[str] | None", verbose: bool) -> list[ResultEntry]:
    """Run the unit tier (unit.runner's GroupFilterError propagates -> exit 2)."""
    from unit import runner as unit_runner

    _failed, raw_entries = unit_runner.run(groups=groups, verbose=verbose)
    return [_to_result_entry(raw, "unit") for raw in raw_entries]


def _infrastructure_driver(verbose_notice: bool = True) -> "tuple[_Driver, str]":
    """Resolve the infrastructure tier driver under the pinned todo-21 seam.

    Returns ``(callable, source)`` where ``callable(groups, verbose)`` returns
    ``(fail_count, entries)``.  Prefers ``infrastructure.infra_runner.cli_run``
    when that module exists and exposes it (todo 21's contract); otherwise
    falls back to ``infrastructure.fixture_runner.run`` (todo 20, fixtures
    only).  Raises ``SystemExit`` (harness error, exit 2) when neither is
    importable.
    """
    cli_run = None
    if importlib.util.find_spec("infrastructure.infra_runner") is not None:
        # Dynamic import: the module may be absent today (todo 21 in flight)
        # and a static ``from infrastructure import infra_runner`` would be a
        # lie to every static analyser; importlib keeps the seam honest.
        try:
            infra_runner = importlib.import_module("infrastructure.infra_runner")
        except ImportError as exc:
            _say(
                "[yellow]infrastructure.infra_runner exists but failed to "
                f"import ({_say_safe(str(exc))}); falling back to fixtures[/yellow]"
            )
        else:
            cli_run = getattr(infra_runner, "cli_run", None)
    if callable(cli_run):
        return cast(_Driver, cli_run), "infra_runner.cli_run"

    def fixtures_only(groups=None, verbose=False) -> "tuple[int, list]":
        """Fallback driver: the todo-20 fixtures runner alone."""
        from infrastructure import fixture_runner

        return cast(
            "tuple[int, list]", fixture_runner.run(groups=groups, verbose=verbose)
        )

    try:
        from infrastructure import fixture_runner  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "suite.cli: infrastructure boot driver unavailable - neither "
            "infrastructure.infra_runner.cli_run (todo 21) nor "
            f"infrastructure.fixture_runner (todo 20) could be imported: {exc}"
        ) from exc
    if verbose_notice:
        _say(
            "[yellow]infrastructure boot driver (infra_runner.cli_run) "
            "unavailable - running fixtures only (boot/distributed groups "
            "land with todo 21)[/yellow]"
        )
    return fixtures_only, "fixture_runner (fixtures-only fallback)"


def _dispatch_infrastructure(
    groups: "list[str] | None", verbose: bool
) -> list[ResultEntry]:
    driver, source = _infrastructure_driver()
    _say(f"[dim]infrastructure driver: {_say_safe(source)}[/dim]")
    _failed, raw_entries = driver(groups=groups, verbose=verbose)
    return [_to_result_entry(raw, "infrastructure") for raw in raw_entries]


def _resolve_configurations(spec: str) -> list[dict]:
    """Resolve ``--conf``: explicit path, path under src/Tests/regression/,
    or a built-in group name (mirrors regression.runner's own resolution)."""
    from suite import config_loader

    direct = Path(spec)
    if direct.is_file():
        return config_loader.load_configurations(direct)
    under_regression = paths.TESTS_DIR / "regression" / spec
    if under_regression.is_file():
        return config_loader.load_configurations(under_regression)
    return config_loader.load_builtin_group(spec)


def _filter_regression(
    configurations: list[dict], groups: "list[str] | None"
) -> list[dict]:
    """Any-substring ``--group`` filter on case-tag (regression is filtered
    here, not inside the runner, so every repeat of the flag is honoured).
    An empty match is a selection error (exit 2) naming the available tags."""
    if not groups:
        return configurations
    matched = [
        conf
        for conf in configurations
        if any(sub in str(conf.get("case-tag", "")) for sub in groups)
    ]
    if not matched:
        available = ", ".join(str(conf.get("case-tag", "")) for conf in configurations)
        raise SystemExit(
            f"suite.cli: no regression entries matched --group {groups!r} "
            f"(available: {available or '<none>'})"
        )
    return matched


def _require_prebuilt_binaries(configurations: list[dict]) -> None:
    """``--no-build`` guard: prebuilt binaries must exist up front.

    With the build phase skipped, a missing binary could never be rebuilt and
    the case would die as a plain run FAILURE (exit 1) inside the runner.
    The pinned semantics say a missing required build product is a harness
    error, hence this orchestrator-side check (exit 2, names every missing
    binary).  With ``--build`` the runner's own post-build assert covers it
    (cmake rebuilds whatever is stale).
    """
    missing = []
    for conf in configurations:
        case_rel = str(conf.get("case", ""))
        binary = suite_build.case_binary(case_rel, conf.get("target"))
        if not binary.is_file():
            missing.append((str(conf.get("case-tag", "")), binary))
    if missing:
        lines = "\n".join(f"  {tag}: {path}" for tag, path in missing)
        raise SystemExit(
            "suite.cli: --no-build requires prebuilt case binaries; "
            f"missing {len(missing)}:\n{lines}"
        )


def _dispatch_regression(
    args: argparse.Namespace, groups: "list[str] | None", verbose: bool
) -> list[ResultEntry]:
    """Resolve ``--conf``, filter, guard binaries, delegate to regression.runner."""
    from regression import runner as regression_runner

    configurations = _filter_regression(
        _resolve_configurations(args.conf or DEFAULT_CONF), groups
    )
    if not args.build:
        _require_prebuilt_binaries(configurations)
    _failed, entries = regression_runner.run(
        configurations,
        store_results=args.store_results,
        verbose=verbose,
        jobs=None,  # --jobs flows through suite.build.set_jobs ONLY (single owner)
        build=args.build,
        group=None,  # --group(s) were applied above (any-substring, repeatable)
    )
    return [_to_result_entry(raw, "regression") for raw in entries]


# ── --list crumbs ─────────────────────────────────────────────────────────


def _list_unit_crumbs() -> None:
    """Unit groups with check counts (reuses unit.runner's own listing helpers)."""
    from unit import runner as unit_runner

    counts = unit_runner._available_groups(unit_runner._unit_dir())
    _say(f"[bold]unit[/bold]: {len(counts)} group(s), {sum(counts.values())} check(s)")
    for group in sorted(counts):
        _say(f"  {_say_safe(group)}: {counts[group]} check(s)")


def _fixture_names() -> list[str]:
    """Fixture names from infrastructure/tests.yaml (no comparison executed).

    Reuses the machinery re-exported by the todo-20 fixture_runner module
    (load_test_list returns ``(tolerance, cases)``; unpack by position --
    its upstream annotation mismatch is pre-existing, see fixture_runner).
    """
    from infrastructure import fixture_runner

    loaded = fixture_runner.run_sph_tests.load_test_list(
        fixture_runner.run_sph_tests._default_tests_yaml
    )
    # load_test_list returns (tolerance, cases); its header misdeclares the
    # return as list[dict] only (pre-existing upstream annotation mismatch,
    # the same situation fixture_runner solves with a positional cast).
    cases = cast("list[dict]", loaded[1])
    return [str(case.get("name", "")) for case in cases]


def _list_infrastructure_crumbs() -> None:
    names = _fixture_names()
    _say(f"[bold]infrastructure[/bold]: {len(names)} fixture(s)")
    for name in names:
        _say(f"  {_say_safe(name)}")
    if importlib.util.find_spec("infrastructure.infra_runner") is None:
        _say("  [dim](boot/distributed groups land with todo 21's infra_runner)[/dim]")


def _list_regression_crumbs(conf_spec: "str | None") -> None:
    from suite import config_loader

    if conf_spec is None:
        groups = [g for g in config_loader.list_builtin_groups() if g != "__init__"]
        _say(f"[bold]regression[/bold]: {len(groups)} built-in group(s) available")
        for group in sorted(groups):
            count = len(config_loader.load_builtin_group(group))
            _say(f"  {_say_safe(group)}: {count} entries")
        return
    # --conf given: list that module's entries (regression.runner --list style).
    configurations = _resolve_configurations(conf_spec)
    _say(f"[bold]regression[/bold]: {_say_safe(conf_spec)} -> {len(configurations)} entries")
    for conf in configurations:
        case_rel = str(conf.get("case", ""))
        target = conf.get("target", f"{Path(case_rel).name}_cuda")
        _say(
            f"  {_say_safe(str(conf.get('case-tag', '')))}"
            f"  case={_say_safe(case_rel)}"
            f"  eval={_say_safe(str(conf.get('evaluation-function', '') or ''))}"
            f"  target={_say_safe(str(target))}"
            f"  init={_say_safe(str(conf.get('init-script', 'init.py')))}"
        )


def _list_crumbs(tiers: list[str], conf_spec: "str | None") -> int:
    """Per-tier crumb listing without executing anything; exit 0."""
    for index, tier in enumerate(tiers):
        if index:
            _say("")
        if tier == "unit":
            _list_unit_crumbs()
        elif tier == "infrastructure":
            _list_infrastructure_crumbs()
        else:
            _list_regression_crumbs(conf_spec)
    return 0


# ── --report reprint mode ─────────────────────────────────────────────────


def _report(log_path: str) -> int:
    """Reprint a session JSON log; exit 0 on successful parse, 2 otherwise.

    The recorded machine line is honored verbatim (no live re-detection);
    pass/fail counts are re-derived from the stored entries but never affect
    the exit code -- reprint always exits 0 on a well-formed session log.
    """
    path = Path(log_path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise SystemExit(f"suite.cli: cannot read session log {path}: {exc}") from exc
    try:
        data = json.loads(text)
    except ValueError as exc:
        raise SystemExit(f"suite.cli: cannot parse session log {path}: {exc}") from exc
    if not isinstance(data, dict) or data.get("schema") != SESSION_SCHEMA:
        raise SystemExit(
            f"suite.cli: {path} is not a schema-{SESSION_SCHEMA} session log"
        )
    raw_entries = data.get("entries")
    if not isinstance(raw_entries, list):
        raise SystemExit(f"suite.cli: {path}: session log has no entries list")
    entries = [ResultEntry.from_json(raw) for raw in raw_entries]

    machine = data.get("machine") or {}
    _say(
        f"[bold]session[/bold] {_say_safe(path.name)}"
        f"  started={_say_safe(str(data.get('started', '?')))}"
        f"  duration_s={_say_safe(str(data.get('duration_s', '?')))}"
    )
    _say(
        f"machine: gpu={_say_safe(str(machine.get('gpu') or '-'))}"
        f"  host={_say_safe(str(machine.get('host') or '-'))}"
    )
    by_tier: dict[str, list[ResultEntry]] = {}
    for entry in entries:
        by_tier.setdefault(entry.tier or "?", []).append(entry)
    perf_fn = None
    gpu = str(machine.get("gpu") or "")
    if gpu:
        try:
            from regression import perf as regression_perf
        except Exception:
            regression_perf = None
        if regression_perf is not None:
            ref_times = regression_perf.load_reference_times(gpu)
            perf_fn = lambda entry: regression_perf.perf_columns(entry, ref_times)
    reporting.print_final_summary(by_tier, perf_fn=perf_fn)
    mismatches = [e for e in entries if _tests_mismatch(e) and not results.result_failed(e)]
    if mismatches:
        # Informational only: the session that PRODUCED this log exited 1
        # under the tests-mismatch verdict pin; reprint itself stays exit 0.
        _say(
            f"[bold red]{len(mismatches)} entr{'y' if len(mismatches) == 1 else 'ies'} "
            "recorded with tests mismatch (tests_passed != tests_total); reprint "
            "always exits 0 on successful parse.[/bold red]"
        )
    return 0


# ── session logs ──────────────────────────────────────────────────────────


def _machine_info(gpu_hint: str) -> dict:
    """Session machine line: detected/override GPU + host (gpu lazy/tolerated)."""
    if gpu_hint:
        return {"gpu": gpu_hint, "host": platform.node()}
    try:
        from regression import perf as regression_perf
    except Exception:
        gpu = ""
    else:
        try:
            gpu = regression_perf.detect_gpu()
        except Exception:
            gpu = ""
    return {"gpu": gpu, "host": platform.node()}


def _write_session_logs(
    all_entries: list[ResultEntry],
    machine: dict,
    started: float,
    log_dir: Path,
    ref_times,
) -> None:
    """Write the pinned ``log_<stamp>.{json,md}`` pair into ``log_dir``.

    The cli owns log-dir creation here (results writers deliberately never
    create directories silently).
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime(_LOG_STAMP_FORMAT, time.localtime(started))
    session = results.session_log(all_entries, machine, started=started)
    json_path = results.write_json_log(session, log_dir / f"log_{stamp}.json")
    md_path = results.write_markdown_log(
        all_entries, log_dir / f"log_{stamp}.md", ref_times=ref_times
    )
    _say(f"session log (json): {_say_safe(str(json_path))}")
    _say(f"session log (md):   {_say_safe(str(md_path))}")


# ── main ──────────────────────────────────────────────────────────────────


def _run_tiers(args: argparse.Namespace, tiers: list[str]) -> int:
    """Serial per-tier dispatch, final summary, session logs, verdict."""
    if args.build and "regression" in tiers:
        # Configure FIRST (cache-reuse no-op for existing trees) so configure
        # failures surface here as harness errors, before any tier output.
        suite_build.ensure_configured(build_type=args.build_type)

    started = time.time()
    by_tier: dict[str, list[ResultEntry]] = {}
    all_entries: list[ResultEntry] = []
    for tier in tiers:  # SERIAL: one tier at a time, single process (D13/G18)
        if tier == "unit":
            entries = _dispatch_unit(args.groups, args.verbose)
        elif tier == "infrastructure":
            entries = _dispatch_infrastructure(args.groups, args.verbose)
        else:
            entries = _dispatch_regression(args, args.groups, args.verbose)
        by_tier[tier] = entries
        all_entries.extend(entries)

    # Informational perf report (regression tier only; ref times loaded once).
    gpu = ""
    ref_times = None
    perf_fn = None
    if "regression" in by_tier:
        from regression import perf as regression_perf

        gpu = regression_perf.detect_gpu(args.gpu)
        ref_times = regression_perf.load_reference_times(gpu)
        if ref_times is None:
            _say(
                f"[dim]perf: no reference times for '{_say_safe(gpu)}', "
                "perf report skipped[/dim]"
            )
        perf_fn = lambda entry: regression_perf.perf_columns(entry, ref_times)

    rc = reporting.print_final_summary(by_tier, perf_fn=perf_fn)

    # Pinned verdict extension: a tests mismatch (tests_passed != tests_total,
    # tests_total > 0) counts as FAIL even though the entry stays phase-PASS
    # (todo 9 contract).  The summary counts it via its own result_failed
    # rule, so an unexplained exit-1 is avoided by this explicit line.
    mismatches = [e for e in all_entries if _tests_mismatch(e) and not results.result_failed(e)]
    if mismatches:
        rc = 1
        _say(
            f"[bold red]{len(mismatches)} entr{'y' if len(mismatches) == 1 else 'ies'} "
            "with tests mismatch (tests_passed != tests_total) -- exit 1.[/bold red]"
        )

    machine = _machine_info(gpu)
    _write_session_logs(
        all_entries, machine, started, Path(args.log_dir), ref_times
    )
    return rc


def main(argv: "list[str] | None" = None) -> int:
    """CLI entry; see the module docstring for the exit-code contract."""
    args = _parse_args(argv)
    try:
        if args.version:
            _say(
                f"tnl-sph unified test suite {SUITE_VERSION} "
                f"(session schema {SESSION_SCHEMA})"
            )
            return 0
        if args.report is not None:
            return _report(args.report)
        tiers = _normalize_tiers(args.tier)
        suite_build.set_jobs(args.jobs)
        if args.list:
            return _list_crumbs(tiers, args.conf)
        return _run_tiers(args, tiers)
    except SystemExit as exc:
        if isinstance(exc.code, str):
            # Harness errors (config load, selection errors, preflight aborts)
            # arrive as SystemExit(str) from the landed contracts.
            print(exc.code, file=sys.stderr)
            return 2
        return int(exc.code or 0)
    except Exception as exc:
        if args.verbose:
            traceback.print_exc()
        else:
            print(f"suite.cli: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

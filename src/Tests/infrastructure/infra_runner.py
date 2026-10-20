"""Infrastructure-tier boot driver of the TNL-SPH unified test suite (todo 21).

Orchestrates the three infrastructure subgroups behind the pinned
``cli_run(groups, verbose) -> (fail_count, entries)`` seam consumed by
``suite.cli._dispatch_infrastructure``:

- ``fixtures``    -- delegated verbatim to
  :func:`infrastructure.fixture_runner.run` (todo 20).
- ``boot``        -- for every discovered non-distributed example: preflight
  ``requires`` -> snapshot templates -> :func:`suite.runner_core.init_case`
  (per-case init script, never ``run.py``) ->
  :func:`suite.build.build_target` -> run the BINARY directly with
  ``--config sources/config.jsonc`` under a wall-clock timeout ->
  :func:`infrastructure.validators.validate_case.validate`.  A timeout with
  the init report already on disk is a BOOT PASS (the report is written at
  init -- that IS the check); a timeout without it is a FAIL.
- ``distributed`` -- for the two ``*_distributed`` cases: init with the
  ``subdomains-x/y`` overrides (Metis G5: their default is 3), then
  ``mpirun -np <N> <binary> --config sources/config.jsonc`` under the same
  timeout discipline, then ``validate(..., distributed=N)`` (rank-file
  count + tiling rule from todo 19).
- ``all`` (or no ``--group``) -- fixtures first, then boot, then
  distributed, strictly SERIALLY (one solver at a time, D13/G18).

Boot entries are native :class:`suite.results.ResultEntry` objects with the
``init``/``build``/``run``/``validate`` phases so the suite log keeps the
per-phase breakdown (fixture entries stay the runner's local dicts, wrapped
by the CLI).

Process-safety discipline: every solver is spawned in its own process group
(``start_new_session=True``); on timeout the whole group is SIGKILLed via
:func:`os.killpg` and reaped, which also reaps mpirun's rank processes -- no
orphans survive the tier.  Templates are snapshotted before init and
restored in ``finally`` per case; a crashed previous run is caught by the
``check_templates_clean`` preflight (Metis G18/H5) which aborts the whole
tier with the dirty files named.

CLI::

    PYTHONPATH=src/Tests python3 -m infrastructure.infra_runner \
        [--list] [--group NAME ...] [--verbose]

Exit codes: 0 = every entry PASS or SKIP, 1 = at least one FAIL entry,
2 = selection/usage errors (argparse, unknown fixtures substring,
templates-clean preflight abort).
"""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

# Self-locate: make src/Tests importable before bootstrapping (same pattern
# as fixture_runner.py / runTests.py; a no-op when PYTHONPATH is preset).
_TESTS_DIR = Path(__file__).resolve().parent.parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))

import suite.bootstrap  # noqa: E402,F401  (single allowed sys.path hook)

from suite import paths  # noqa: E402
from suite import build as suite_build  # noqa: E402
from suite import runner_core  # noqa: E402
from suite.results import PhaseResult, ResultEntry, result_failed  # noqa: E402

from infrastructure import fixture_runner  # noqa: E402
from infrastructure.boot_configs import BOOT_CONFIGS  # noqa: E402
from infrastructure.validators.validate_case import validate as _validate_report  # noqa: E402

BOLD_GREEN = "\033[1;32m"
BOLD_RED = "\033[1;31m"
BOLD_YELLOW = "\033[1;33m"
DIM = "\033[2m"
RESET = "\033[0m"

TIER = "infrastructure"
KNOWN_GROUPS = ("fixtures", "boot", "distributed", "all")

GENCASE_BINARY = paths.TOOLS_DIR / "3rdparty" / "GenCase_linux64"
GENCASE_SKIP_REASON = "requires gencase"
MPIRUN_SKIP_REASON = "requires mpirun (not found on PATH)"

# Cap for subprocess output tails recorded into phase logs.
_LOG_TAIL = 2000


# ── discovery ─────────────────────────────────────────────────────────────


def discover_cases() -> list[str]:
    """Scan ``examples/*/`` for dirs carrying both ``case.h`` and
    ``template/`` (Metis G3/G6: NOT ``init.py`` -- the 3 MR/LTS dirs have
    none).  Returns examples-relative paths, sorted and stable."""
    cases: list[str] = []
    for family in sorted(p for p in paths.EXAMPLES_DIR.iterdir() if p.is_dir()):
        for case_dir in sorted(p for p in family.iterdir() if p.is_dir()):
            if (case_dir / "case.h").is_file() and (case_dir / "template").is_dir():
                cases.append(f"{family.name}/{case_dir.name}")
    return cases


# ── preflight ─────────────────────────────────────────────────────────────


def _preflight_templates() -> None:
    """Metis G18/H5 preflight: abort when a previous crashed run left a
    tracked ``template/config.h`` dirty.  Raises ``SystemExit`` (harness
    error, CLI exit 2) naming every dirty file."""
    dirty = runner_core.check_templates_clean()
    if dirty:
        raise SystemExit(
            "infra_runner: refusing to run with dirty template/config.h "
            "(a previous run may have crashed); restore them first:\n  "
            + "\n  ".join(dirty)
        )


def _requires_missing(row: dict) -> str | None:
    """Return the SKIP reason when a row's ``requires`` cannot be met."""
    missing = []
    for requirement in row.get("requires", []):
        if requirement == "gencase":
            if not (GENCASE_BINARY.is_file() and os.access(GENCASE_BINARY, os.X_OK)):
                missing.append(
                    f"{GENCASE_SKIP_REASON} ({GENCASE_BINARY} missing or not executable)"
                )
        else:  # pragma: no cover - no other requirement exists today
            missing.append(f"requires {requirement} (unknown requirement)")
    return "; ".join(missing) if missing else None


# ── phase helpers ─────────────────────────────────────────────────────────


def _phases() -> dict:
    return {}


def _mark(phases: dict, name: str, status: str, seconds: float = 0.0,
          returncode: int = 0, log: str = "") -> None:
    phases[name] = PhaseResult(
        status=status, seconds=round(seconds, 3), returncode=returncode, log=log
    )


def _skip_remaining(phases: dict, names: tuple[str, ...], log: str = "") -> None:
    for name in names:
        if name not in phases:
            _mark(phases, name, "SKIP", log=log)


def _tail(output: str) -> str:
    return output[-_LOG_TAIL:] if output else ""


def _wait_stable(receipts: "list[Path]", settle_s: float = 0.3) -> bool:
    """True when every receipt exists, is non-empty and size-stable across
    ``settle_s`` (guards against racing the report's own file dump)."""
    sizes = []
    for path in receipts:
        try:
            size = path.stat().st_size
        except OSError:
            return False
        if size <= 0:
            return False
        sizes.append(size)
    time.sleep(settle_s)
    for path, size in zip(receipts, sizes):
        try:
            if path.stat().st_size != size:
                return False
        except OSError:
            return False
    return True


def _run_boot(cmd: list[str], cwd: Path, timeout_s: float,
              receipts: "list[Path]") -> "tuple[str, int, float, str]":
    """Run ``cmd`` in ``cwd``; stop as soon as the boot receipt stabilizes.

    Returns ``(output, returncode, seconds, outcome)`` with outcome one of
    ``"completed"`` (process exited on its own), ``"receipt"`` (all receipt
    files present and stable -> process group SIGKILLed, this is the normal
    boot path: the report is written at init and long simulations are never
    awaited), or ``"timeout"`` (deadline hit first -> group SIGKILLed).
    The child runs in its own session, so killpg also reaps mpirun's rank
    processes -- no orphans survive a boot.  SIGKILL receipts surface as
    ``returncode == -signal.SIGKILL``.
    """
    started = time.perf_counter()
    deadline = started + timeout_s
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        start_new_session=True,
    )
    while True:
        try:
            output, _ = proc.communicate(timeout=0.5)
            return (output or "", int(proc.returncode or 0),
                    time.perf_counter() - started, "completed")
        except subprocess.TimeoutExpired:
            pass
        if _wait_stable(receipts):
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            output, _ = proc.communicate()
            return (output or "", int(proc.returncode),
                    time.perf_counter() - started, "receipt")
        if time.perf_counter() >= deadline:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            output, _ = proc.communicate()
            return (output or "", int(proc.returncode),
                    time.perf_counter() - started, "timeout")


# ── boot (single-resolution examples) ─────────────────────────────────────


def _init_conf(rel: str, row: dict) -> dict:
    """Merge the row's init overrides into the conf handed to init_case.

    RESERVED-keys semantics (runner_core): suite fields stay OUT of the
    forwarded argv; every remaining item (the overrides) is forwarded as
    ``--<key> <value>``; ``init-script`` selects the script itself.
    """
    return {
        "case": rel,
        "case-tag": rel,
        "init-script": row["init-script"],
        **row.get("init-overrides", {}),
    }


def _report_receipt(case_dir: Path) -> Path:
    """The boot receipt checked on timeout (the task pins results/)."""
    return case_dir / "results" / "init_report.json"


def _boot_case(rel: str, verbose: bool) -> ResultEntry:
    """Boot one non-distributed example; always returns an entry (never
    aborts the tier)."""
    row = BOOT_CONFIGS[rel]
    case_dir = paths.EXAMPLES_DIR / rel
    entry = ResultEntry(tag=rel, group="boot", tier=TIER)
    phases = _phases()

    # 0) requirements preflight BEFORE anything else.
    missing = _requires_missing(row)
    if missing is not None:
        entry.skip_reason = missing
        _skip_remaining(phases, ("init", "build", "run", "validate"), log=missing)
        entry.phases = phases
        return entry

    snapshot = runner_core.snapshot_tracked_templates(case_dir)
    try:
        # 1) init
        started = time.perf_counter()
        try:
            runner_core.init_case(case_dir, _init_conf(rel, row), verbose)
            _mark(phases, "init", "PASS", time.perf_counter() - started)
        except RuntimeError as exc:
            _mark(phases, "init", "FAIL", time.perf_counter() - started,
                  returncode=1, log=str(exc))
            entry.error = str(exc)
            _skip_remaining(phases, ("build", "run", "validate"))
            entry.phases = phases
            return entry

        # 2) build
        started = time.perf_counter()
        try:
            suite_build.build_target(rel, target=row["target"])
            _mark(phases, "build", "PASS", time.perf_counter() - started)
        except RuntimeError as exc:
            message = str(exc)
            if "No rule to make target" in message:
                # Target not registered in this build tree (e.g. the
                # examples/*/CMakeLists.txt add_subdirectory line is
                # commented out) -- an environment state, not a build
                # defect: SKIP with the exact cause named.
                reason = (
                    f"build target '{row['target']}' not configured in "
                    f"{paths.BUILD_DIR} (case's CMake add_subdirectory is "
                    "commented out)"
                )
                entry.skip_reason = reason
                _mark(phases, "build", "SKIP", time.perf_counter() - started,
                      log=reason)
                _skip_remaining(phases, ("run", "validate"), log=reason)
            else:
                _mark(phases, "build", "FAIL", time.perf_counter() - started,
                      returncode=1, log=_tail(message))
                entry.error = message.splitlines()[0]
                _skip_remaining(phases, ("run", "validate"))
            entry.phases = phases
            return entry

        # 3) run the binary directly (NEVER bare: --config <path> only)
        binary = suite_build.case_binary(rel, row["target"])
        receipt = _report_receipt(case_dir)
        receipt.unlink(missing_ok=True)  # no stale receipt may pass for this boot
        output, returncode, seconds, outcome = _run_boot(
            [str(binary), "--config", "sources/config.jsonc"],
            cwd=case_dir,
            timeout_s=float(row["timeout_s"]),
            receipts=[receipt],
        )
        if outcome == "receipt":
            _mark(
                phases, "run", "PASS", seconds, returncode=returncode,
                log=f"boot receipt {receipt.name} "
                    f"(proc SIGKILLed, rc={returncode})",
            )
        elif outcome == "timeout":
            if receipt.is_file():
                _mark(
                    phases, "run", "PASS", seconds, returncode=returncode,
                    log=f"wall-clock timeout after boot receipt "
                        f"{receipt.name} (proc SIGKILLed, rc={returncode})",
                )
            else:
                entry.error = "timeout without init report"
                _mark(phases, "run", "FAIL", seconds, returncode=returncode,
                      log=_tail(output))
                _skip_remaining(phases, ("validate",))
                entry.phases = phases
                return entry
        elif returncode != 0:
            entry.error = f"binary exited with code {returncode}"
            _mark(phases, "run", "FAIL", seconds, returncode=returncode,
                  log=_tail(output))
            _skip_remaining(phases, ("validate",))
            entry.phases = phases
            return entry
        elif not receipt.is_file():
            entry.error = f"init report missing after run: {receipt}"
            _mark(phases, "run", "FAIL", seconds, returncode=0,
                  log="simulation completed but left no init report")
            _skip_remaining(phases, ("validate",))
            entry.phases = phases
            return entry
        else:
            _mark(phases, "run", "PASS", seconds, returncode=0,
                  log="simulation completed within timeout")

        # 4) validate the init report
        started = time.perf_counter()
        try:
            findings = _validate_report(rel)
        except SystemExit as exc:
            entry.error = f"validate_case: {exc.code}"
            _mark(phases, "validate", "FAIL", time.perf_counter() - started,
                  returncode=1, log=str(exc.code))
            entry.phases = phases
            return entry
        except ValueError as exc:
            entry.error = f"validate_case: unreadable init report ({exc})"
            _mark(phases, "validate", "FAIL", time.perf_counter() - started,
                  returncode=1, log=str(exc))
            entry.phases = phases
            return entry
        failed = [f for f in findings if f.severity == "FAIL"]
        if failed:
            names = ", ".join(f.rule for f in failed)
            entry.error = f"validation failed ({names})"
            first = str(failed[0].message)
            _mark(phases, "validate", "FAIL", time.perf_counter() - started,
                  returncode=len(failed), log=f"{names}: {first}")
        else:
            _mark(phases, "validate", "PASS", time.perf_counter() - started)
        entry.phases = phases
        return entry
    finally:
        runner_core.restore_tracked_templates(case_dir, snapshot)


# ── distributed (MPI boots) ───────────────────────────────────────────────


def _distributed_case(rel: str, verbose: bool) -> ResultEntry:
    """Boot one ``*_distributed`` case under ``mpirun -np <N>``; always
    returns an entry."""
    row = BOOT_CONFIGS[rel]
    dist = row["distributed"]
    np_ranks = int(dist["np"])
    case_dir = paths.EXAMPLES_DIR / rel
    entry = ResultEntry(tag=rel, group="distributed", tier=TIER)
    phases = _phases()

    missing = _requires_missing(row)
    if missing is None and shutil.which("mpirun") is None:
        missing = MPIRUN_SKIP_REASON
    if missing is not None:
        entry.skip_reason = missing
        _skip_remaining(phases, ("init", "build", "run", "validate"), log=missing)
        entry.phases = phases
        return entry

    conf = _init_conf(rel, row)
    conf.update(dist.get("extra-init-overrides", {}))

    snapshot = runner_core.snapshot_tracked_templates(case_dir)
    try:
        started = time.perf_counter()
        try:
            runner_core.init_case(case_dir, conf, verbose)
            _mark(phases, "init", "PASS", time.perf_counter() - started)
        except RuntimeError as exc:
            _mark(phases, "init", "FAIL", time.perf_counter() - started,
                  returncode=1, log=str(exc))
            entry.error = str(exc)
            _skip_remaining(phases, ("build", "run", "validate"))
            entry.phases = phases
            return entry

        started = time.perf_counter()
        try:
            suite_build.build_target(rel, target=row["target"])
            _mark(phases, "build", "PASS", time.perf_counter() - started)
        except RuntimeError as exc:
            message = str(exc)
            if "No rule to make target" in message:
                reason = (
                    f"build target '{row['target']}' not configured in "
                    f"{paths.BUILD_DIR} (case's CMake add_subdirectory is "
                    "commented out)"
                )
                entry.skip_reason = reason
                _mark(phases, "build", "SKIP", time.perf_counter() - started,
                      log=reason)
                _skip_remaining(phases, ("run", "validate"), log=reason)
            else:
                _mark(phases, "build", "FAIL", time.perf_counter() - started,
                      returncode=1, log=_tail(message))
                entry.error = message.splitlines()[0]
                _skip_remaining(phases, ("run", "validate"))
            entry.phases = phases
            return entry

        binary = suite_build.case_binary(rel, row["target"])
        receipt = _report_receipt(case_dir)
        rank_receipts = [
            case_dir / "results" / f"init_report_subdomain{rank}.json"
            for rank in range(np_ranks)
        ]
        for path in [receipt, *rank_receipts]:
            path.unlink(missing_ok=True)  # no stale receipt may pass
        output, returncode, seconds, outcome = _run_boot(
            ["mpirun", "-np", str(np_ranks), str(binary),
             "--config", "sources/config.jsonc"],
            cwd=case_dir,
            timeout_s=float(row["timeout_s"]),
            receipts=rank_receipts,
        )
        if outcome == "receipt":
            _mark(
                phases, "run", "PASS", seconds, returncode=returncode,
                log=f"boot receipt: {np_ranks} rank reports "
                    f"(proc group SIGKILLed, rc={returncode})",
            )
        elif outcome == "timeout":
            if receipt.is_file() or any(p.is_file() for p in rank_receipts):
                _mark(
                    phases, "run", "PASS", seconds, returncode=returncode,
                    log=f"wall-clock timeout after boot receipt "
                        f"(proc group SIGKILLed, rc={returncode})",
                )
            else:
                entry.error = "timeout without init report"
                _mark(phases, "run", "FAIL", seconds, returncode=returncode,
                      log=_tail(output))
                _skip_remaining(phases, ("validate",))
                entry.phases = phases
                return entry
        elif returncode != 0:
            entry.error = f"mpirun exited with code {returncode}"
            _mark(phases, "run", "FAIL", seconds, returncode=returncode,
                  log=_tail(output))
            _skip_remaining(phases, ("validate",))
            entry.phases = phases
            return entry
        elif not (receipt.is_file() or any(p.is_file() for p in rank_receipts)):
            entry.error = "init report missing after run"
            _mark(phases, "run", "FAIL", seconds, returncode=0,
                  log="simulation completed but left no init report")
            _skip_remaining(phases, ("validate",))
            entry.phases = phases
            return entry
        else:
            _mark(phases, "run", "PASS", seconds, returncode=0,
                  log="simulation completed within timeout")

        started = time.perf_counter()
        try:
            findings = _validate_report(rel, distributed=np_ranks)
        except SystemExit as exc:
            entry.error = f"validate_case: {exc.code}"
            _mark(phases, "validate", "FAIL", time.perf_counter() - started,
                  returncode=1, log=str(exc.code))
            entry.phases = phases
            return entry
        except ValueError as exc:
            entry.error = f"validate_case: unreadable init report ({exc})"
            _mark(phases, "validate", "FAIL", time.perf_counter() - started,
                  returncode=1, log=str(exc))
            entry.phases = phases
            return entry
        failed = [f for f in findings if f.severity == "FAIL"]
        if failed:
            names = ", ".join(f.rule for f in failed)
            entry.error = f"validation failed ({names})"
            first = str(failed[0].message)
            _mark(phases, "validate", "FAIL", time.perf_counter() - started,
                  returncode=len(failed), log=f"{names}: {first}")
        else:
            _mark(phases, "validate", "PASS", time.perf_counter() - started)
        entry.phases = phases
        return entry
    finally:
        runner_core.restore_tracked_templates(case_dir, snapshot)


# ── group orchestration ───────────────────────────────────────────────────


def _resolve_groups(groups: "list[str] | None") -> "tuple[list[str] | None, bool, bool]":
    """Map the CLI group selection to ``(fixture_substrings, boot?, distributed?)``.

    ``fixture_substrings`` is ``None`` when the fixtures group is NOT
    selected, ``[]`` when selected unfiltered, and a list of substrings
    otherwise (any unknown group token is forwarded to the fixtures runner
    as a substring filter -- its FixtureFilterError propagates as a
    selection error, exit 2 via the CLI).
    """
    if not groups:
        return [], True, True
    if "all" in groups:
        return [], True, True
    fixtures: "list[str] | None" = [] if "fixtures" in groups else None
    boot = "boot" in groups
    distributed = "distributed" in groups
    substrings = [g for g in groups if g not in KNOWN_GROUPS]
    if substrings:
        fixtures = substrings
    return fixtures, boot, distributed


def _print_entry(entry: ResultEntry, verbose: bool) -> None:
    """One board line per boot/distributed entry."""
    phases = entry.phases
    detail = ", ".join(
        f"{name} {phases[name].seconds:.1f}s"
        for name in ("init", "build", "run", "validate")
        if name in phases and phases[name].status not in ("SKIP",)
    )
    if entry.error:
        label, color = "FAIL", BOLD_RED
    elif entry.skip_reason:
        label, color = "SKIP", BOLD_YELLOW
    else:
        label, color = "PASS", BOLD_GREEN
    mark = {"PASS": "✓", "FAIL": "✗", "SKIP": "⊘"}[label]
    reason = entry.error or entry.skip_reason or ""
    extra = f"  {DIM}{reason}{RESET}" if reason else ""
    if detail:
        extra += f"  {DIM}({detail}){RESET}"
    print(f"{color}  {mark}  {label} {entry.tag}{extra}{RESET}", flush=True)
    if verbose:
        for name, phase in phases.items():
            if phase.log:
                print(f"{DIM}      {name}: {phase.log[:400]}{RESET}", flush=True)


def cli_run(groups: "list[str] | None" = None, verbose: bool = False) -> "tuple[int, list]":
    """The pinned todo-21 seam for ``suite.cli`` (see module docstring).

    Returns ``(fail_count, entries)`` where entries mixes fixture-runner
    dicts and boot/distributed :class:`ResultEntry` objects.  FAIL entries
    count toward ``fail_count``; SKIPs never do.
    """
    _preflight_templates()
    fixture_substrings, want_boot, want_distributed = _resolve_groups(groups)
    discovered = discover_cases()

    fail_count = 0
    entries: list = []

    # 1) fixtures (delegated verbatim; its entries are local dicts).
    if fixture_substrings is not None:
        fixture_fails, fixture_entries = fixture_runner.run(
            groups=fixture_substrings or None, verbose=verbose
        )
        fail_count += fixture_fails
        entries.extend(fixture_entries)
        if not verbose and fixture_entries:
            f_pass = sum(1 for e in fixture_entries if e["status"] == "PASS")
            f_fail = sum(1 for e in fixture_entries if e["status"] == "FAIL")
            f_skip = sum(1 for e in fixture_entries if e["status"] == "SKIP")
            print(
                f"  fixtures: {f_pass} passed, {f_fail} failed, "
                f"{f_skip} skipped, {len(fixture_entries)} total",
                flush=True,
            )

    # 2) boot: discovered non-distributed cases.
    if want_boot:
        for rel in discovered:
            row = BOOT_CONFIGS.get(rel)
            if row is None:
                # Discovery finds exactly the pinned map by construction;
                # keep the mismatch visible instead of guessing fields.
                entry = ResultEntry(tag=rel, group="boot", tier=TIER)
                entry.skip_reason = "no BOOT_CONFIGS row (add it explicitly)"
                _skip_remaining(entry.phases, ("init", "build", "run", "validate"),
                                log=entry.skip_reason)
                entries.append(entry)
                _print_entry(entry, verbose)
                continue
            if row.get("distributed"):
                continue  # handled by the distributed group
            entry = _boot_case(rel, verbose)
            if result_failed(entry):
                fail_count += 1
            entries.append(entry)
            _print_entry(entry, verbose)

    # 3) distributed: the MPI boots.
    if want_distributed:
        if shutil.which("mpirun") is None:
            print(f"{BOLD_YELLOW}  ⊘  SKIP distributed group: "
                  f"{MPIRUN_SKIP_REASON}{RESET}", flush=True)
        for rel in discovered:
            row = BOOT_CONFIGS.get(rel)
            if row is None or not row.get("distributed"):
                continue
            entry = _distributed_case(rel, verbose)
            if result_failed(entry):
                fail_count += 1
            entries.append(entry)
            _print_entry(entry, verbose)

    return fail_count, entries


# ── CLI ───────────────────────────────────────────────────────────────────


def _list_cases() -> None:
    """Print the discovered cases with their per-case boot fields."""
    discovered = discover_cases()
    boot = [r for r in discovered if not BOOT_CONFIGS.get(r, {}).get("distributed")]
    dist = [r for r in discovered if BOOT_CONFIGS.get(r, {}).get("distributed")]
    print(f"discovered {len(discovered)} cases "
          f"({len(boot)} boot + {len(dist)} distributed):")
    for rel in discovered:
        row = BOOT_CONFIGS.get(rel)
        mark = "" if row is not None else "  [NO BOOT_CONFIGS ROW]"
        print(f"  {rel}{mark}")
        if row is None:
            continue
        print(f"      init-script: {row['init-script']}")
        print(f"      target:      {row['target']}")
        print(f"      overrides:   {row.get('init-overrides', {})}")
        print(f"      timeout_s:   {row['timeout_s']}")
        print(f"      requires:    {row.get('requires', [])}")
        if row.get("distributed"):
            dist_row = row["distributed"]
            extra = dist_row.get("extra-init-overrides", {})
            print(f"      distributed: np={dist_row['np']}, "
                  f"extra-init-overrides={extra}")


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(
        prog="infrastructure.infra_runner",
        description=(__doc__ or "").splitlines()[0],
    )
    parser.add_argument(
        "--group", dest="groups", action="append", default=None,
        metavar="NAME",
        help="fixtures | boot | distributed | all, or a fixture-name "
        "substring (repeatable; default: all three subgroups)",
    )
    parser.add_argument("--list", action="store_true",
                        help="print discovered cases with per-case fields")
    parser.add_argument("--verbose", action="store_true",
                        help="stream per-phase detail")
    args = parser.parse_args(argv)

    if args.list:
        _list_cases()
        return 0

    try:
        fail_count, entries = cli_run(groups=args.groups, verbose=args.verbose)
    except fixture_runner.FixtureFilterError as exc:
        print(exc, file=sys.stderr)
        return 2

    passed, failed, skipped = 0, 0, 0
    for entry in entries:
        if isinstance(entry, ResultEntry):
            if result_failed(entry):
                failed += 1
            elif entry.skip_reason:
                skipped += 1
            else:
                passed += 1
        else:
            status = entry["status"]
            if status == "FAIL":
                failed += 1
            elif status == "SKIP":
                skipped += 1
            else:
                passed += 1
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped, "
          f"{len(entries)} total", flush=True)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

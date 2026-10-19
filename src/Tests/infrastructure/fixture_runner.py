"""Infrastructure-tier fixtures runner for the TNL-SPH unified test suite (stdlib only).

Folds the moved multiresolution fixture machinery (todo 2:
``infrastructure/run_sph_tests.py`` + ``tests.yaml`` + ``fixtures/``) into the
infra tier of the suite *as a library* (todo 20), without changing the
comparison logic: :func:`run` replicates
``infrastructure.run_sph_tests.main()`` semantics WITHOUT calling
``sys.exit`` - ``init=True`` calls its ``do_init(config_name)``,
``run=True`` calls its ``do_run(config_name)`` and ``compare=True`` (the
default) calls its ``run_batch_tests(_default_tests_yaml, cli_tol=...)``.
``parse_sph_log.py`` and ``tests.yaml`` (21 entries, global `1.0e-9`
tolerance + per-test overrides, all preserved verbatim) are untouched;
the standalone script keeps working byte-identically:
``python3 src/Tests/infrastructure/run_sph_tests.py --tests src/Tests/infrastructure/tests.yaml``.

Import contract: this module imports the machinery via ``import
suite.bootstrap`` (single allowed sys.path hook; puts ``src/Tests`` +
``src/tools`` on ``sys.path``) followed by the package import
``from infrastructure import run_sph_tests``.

CLI::

    python3 -m infrastructure.fixture_runner [--init] [--run]
                                             [--config-name N]
                                             [--group SUBSTR ...]
                                             [--verbose]

Prints one line per fixture (``✓`` PASS / ``✗`` FAIL / ``⊘`` SKIP,
including a one-line reason for non-PASS entries) and a final
``Results: X passed, Y failed, Z skipped, S total`` summary. Exit codes:
``0`` iff Y == 0, else ``1``; argparse errors AND a ``--group`` filter
that matches no fixture exit ``2`` (mirroring ``unit.runner``'s selection
-error-is-not-a-test-failure convention).

KNOWN STATE (kept visible, not hidden): the pre-existing
``dummy2D-lts-center`` failure documented in the unified-test-suite
notepad ``problems.md`` (``buffers.1.number_of_cells_in_zone``
expected 2272 vs actual 1692, stale ref/log pair after the LTS buffer
rework on this branch) intentionally surfaces as a FAIL entry here, so
the default run exits 1 with ``20 passed, 1 failed, 0 skipped`` until
the user fixes the branch-side reference. The runner must NOT mask it.

Result entries: :func:`run` returns ``(fail_count, entries)`` where
``fail_count`` is the number of FAIL entries (SKIPs do not count) and
each entry is::

    {"tag": <yaml name>, "group": "fixtures", "tier": "infrastructure",
     "status": "PASS"|"FAIL"|"SKIP", "seconds": <float>, "error": <str>}

- Entries follow the ``tests.yaml`` order (stable over the full 21),
  merging batch-compare results back into the SKIP placeholders.
- A fixture whose ``log`` file is missing becomes SKIP with
  ``error == "log missing (init+run needed)"``, ``seconds == 0.0``, and
  is NOT fed to ``run_batch_tests`` (it would be a FAIL there); run
  ``--init --config-name N`` + ``--run --config-name N`` first to create it.
- ``groups`` filters yaml names by SUBSTRING match (any of the given
  substrings; e.g. ``--group lts-center`` selects exactly
  ``dummy2D-lts-center``, ``--group corner`` selects all corners).
- ``seconds`` is the wall time of the single comparison (0.0 for SKIP).
- ``error`` holds the captured (ANSI-stripped) single-test output for
  FAIL entries - first ``✗`` line names the cause (e.g. a missing ref or
  ``FAIL (N difference(s)):`` followed by the diff paths) - else empty.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import re
import sys
import tempfile
import time
from pathlib import Path
from typing import cast

# Self-contained module use: make `src/Tests` importable before bootstrapping
# (same self-locate pattern as `src/Tests/runTests.py`; bootstrap itself only
# inserts missing paths, so this is a no-op when PYTHONPATH is preset).
_TESTS_DIR = Path(__file__).resolve().parent.parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))

import suite.bootstrap  # noqa: E402,F401

# run_sph_tests imports `parse_sph_log` as a bare script-relative module (it
# resolves when run as a script because the script dir lands on sys.path); as a
# package import the infra dir must be importable itself. File-relative shim,
# CWD-independent (same pattern as regression/configurationsToTest.py:1-4).
_INFRA_DIR = Path(__file__).resolve().parent
if str(_INFRA_DIR) not in sys.path:
    sys.path.insert(0, str(_INFRA_DIR))

from infrastructure import run_sph_tests  # noqa: E402

# ── Fancy output (ported from src/UnitTests/multiresolution/run_sph_tests.py:41-52)
BOLD_GREEN = "\033[1;32m"
BOLD_RED = "\033[1;31m"
BOLD_YELLOW = "\033[1;33m"
RESET = "\033[0m"

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

SKIP_REASON_MISSING_LOG = "log missing (init+run needed)"


class FixtureFilterError(Exception):
    """Raised by :func:`run` when a group filter matches no yaml entry."""

    def __init__(self, groups: "list[str]", available: "list[str]") -> None:
        self.groups = list(groups)
        self.available = sorted(available)
        super().__init__(
            f"group filter {self.groups} matched no fixtures; "
            f"available fixtures: {', '.join(self.available) or '(none)'}"
        )


def _strip_ansi(text: str) -> str:
    return ANSI_RE.sub("", text)


def _first_failure_line(output: str) -> str:
    """Concise one-line cause from a captured run_single_test output."""
    for line in output.splitlines():
        if "✗" in line:
            return line.strip()
    for line in output.splitlines():
        line = line.strip()
        if line and not line.startswith("─"):
            return line
    return ""


def run(
    groups: "list[str] | None" = None,
    verbose: bool = False,
    *,
    init: bool = False,
    run: bool = False,
    compare: bool = True,
    config_name: "str | None" = None,
) -> "tuple[int, list[dict]]":
    """Replicate ``run_sph_tests.main()`` semantics; return ``(fail_count, entries)``.

    ``init``/``run`` delegate to ``run_sph_tests.do_init``/``do_run`` with
    ``config_name`` (a failing subprocess surfaces as ``SystemExit`` from
    ``run_sph_tests._subprocess_run``, propagated unchanged). ``compare``
    (default True) delegates the actual comparison to
    ``run_sph_tests.run_batch_tests`` over a filtered temporary test list
    written next to ``tests.yaml`` (log/ref paths are resolved relative to
    the list file's directory, so the temp list must live there) carrying
    the yaml's global tolerance and the selected cases' verbatim per-test
    overrides; ``cli_tol=global_tol`` reproduces ``main()``'s default
    resolution (``tol = cli_tol if cli_tol != 1e-9 else global_tol``).
    Per-case outcomes are captured by temporarily wrapping
    ``run_sph_tests.run_single_test`` (restored in ``finally``) so the
    batch machinery itself is never reimplemented.
    """
    if init:
        run_sph_tests.do_init(config_name)
    if run:
        run_sph_tests.do_run(config_name)

    fail_count = 0
    entries: "list[dict]" = []
    if not compare:
        return 0, entries

    tests_yaml = run_sph_tests._default_tests_yaml
    # load_test_list returns (global_tolerance, cases); its header declares
    # list[dict] only, so unpack by position with casts (upstream annotation
    # mismatch is pre-existing LSP noise, not touched here).
    loaded = run_sph_tests.load_test_list(tests_yaml)
    global_tol = cast(float, loaded[0])
    cases = cast("list[dict]", loaded[1])
    all_names = [case.get("name", "") for case in cases]

    if groups:
        cases = [
            case
            for case in cases
            if any(substring in case.get("name", "") for substring in groups)
        ]
        if not cases:
            raise FixtureFilterError(groups, all_names)

    # SKIP fixtures whose log is missing; compare the rest in one batch.
    runnable: "list[dict]" = []
    for case in cases:
        name = case.get("name", "")
        if (tests_yaml.parent / case["log"]).exists():
            runnable.append(case)
            entries.append({"_run": name})  # placeholder, replaced below
        else:
            entries.append(
                {
                    "tag": name,
                    "group": "fixtures",
                    "tier": "infrastructure",
                    "status": "SKIP",
                    "seconds": 0.0,
                    "error": SKIP_REASON_MISSING_LOG,
                }
            )

    if not runnable:
        return 0, entries

    captured: "dict[str, tuple[bool, float, str]]" = {}
    original_single_test = run_sph_tests.run_single_test

    def _capture_single_test(log_path, ref_path, tol=1e-9, name=""):
        key = name or log_path.name
        buffer = io.StringIO()
        started = time.perf_counter()
        with contextlib.redirect_stdout(buffer):
            ok = original_single_test(log_path, ref_path, tol=tol, name=name)
        seconds = time.perf_counter() - started
        output = buffer.getvalue()
        if verbose:
            sys.stdout.write(output)
        captured[key] = (ok, seconds, _strip_ansi(output).strip())
        return ok

    # JSON is an accepted list format for run_batch_tests and preserves the
    # case dicts verbatim (name/log/ref[/tolerance]).
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=str(tests_yaml.parent),
        prefix=".fixture_runner-",
        suffix=".json",
        delete=False,
    ) as tmp_file:
        json.dump({"tolerance": global_tol, "tests": runnable}, tmp_file)
        tmp_path = Path(tmp_file.name)

    run_sph_tests.run_single_test = _capture_single_test
    try:
        if verbose:
            run_sph_tests.run_batch_tests(tmp_path, cli_tol=global_tol)
        else:
            with contextlib.redirect_stdout(io.StringIO()):
                run_sph_tests.run_batch_tests(tmp_path, cli_tol=global_tol)
    finally:
        run_sph_tests.run_single_test = original_single_test
        tmp_path.unlink(missing_ok=True)

    # Merge batch results back into the placeholders, keeping yaml order.
    for i, entry in enumerate(entries):
        name = entry.pop("_run", None)
        if name is None:
            continue  # SKIP entry already complete
        ok, seconds, output = captured[name]
        status = "PASS" if ok else "FAIL"
        if not ok:
            fail_count += 1
        entries[i] = {
            "tag": name,
            "group": "fixtures",
            "tier": "infrastructure",
            "status": status,
            "seconds": round(seconds, 6),
            "error": "" if ok else output,
        }

    return fail_count, entries


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(
        prog="infrastructure.fixture_runner",
        description="Infra-tier runner over the multiresolution fixtures "
        "(wraps run_sph_tests as a library; see module docstring).",
    )
    parser.add_argument(
        "--init",
        action="store_true",
        help="Initialize fixture configurations (run_sph_tests.do_init) "
        "before comparing.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run fixture configurations (run_sph_tests.do_run) before "
        "comparing.",
    )
    parser.add_argument(
        "--config-name",
        metavar="N",
        choices=run_sph_tests.CONFIG_NAMES,
        help="Target a single configuration (applies to --init/--run).",
    )
    parser.add_argument(
        "--group",
        dest="groups",
        action="append",
        default=None,
        metavar="SUBSTR",
        help="Compare only fixtures whose yaml name contains SUBSTR; "
        "may be given multiple times.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Stream the underlying per-test comparison output as well.",
    )
    args = parser.parse_args(argv)

    try:
        fail_count, entries = run(
            groups=args.groups,
            verbose=args.verbose,
            init=args.init,
            run=args.run,
            compare=True,
            config_name=args.config_name,
        )
    except FixtureFilterError as error:
        # Distinct from a failing fixture: the selection itself was wrong.
        print(error, file=sys.stderr)
        return 2

    passed = failed = skipped = 0
    for entry in entries:
        status = entry["status"]
        if status == "PASS":
            passed += 1
            print(f"{BOLD_GREEN}  ✓  {entry['tag']}{RESET}")
        elif status == "FAIL":
            failed += 1
            print(
                f"{BOLD_RED}  ✗  {entry['tag']} - "
                f"{_first_failure_line(entry['error'])}{RESET}"
            )
        else:  # SKIP
            skipped += 1
            print(f"{BOLD_YELLOW}  ⊘  {entry['tag']} - {entry['error']}{RESET}")
    total = len(entries)
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped, {total} total")
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        sys.exit(main())
    sys.exit(130)

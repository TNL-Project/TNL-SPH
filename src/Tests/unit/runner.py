"""Unit-tier runner for the TNL-SPH unified test suite (stdlib only).

Runs :mod:`unit.collector` over the ``unit/`` directory, executes every
check sequentially in deterministic order (sorted by file path, then by
check name), and reports ``run_sph_tests.py``-style colored lines.

CLI (ctest later invokes ``--here``)::

    python3 -m unit.runner [--group NAME ...] [--list] [--verbose] [--here]
    python3 -m unit.runner --selftest

Exit codes: ``0`` when every executed check passed (including an empty
discovery, which prints ``0 checks``); ``1`` when at least one check
failed; ``2`` when a ``--group`` filter matched no check (the available
groups are printed so the selection failure is never silent).

Result-entry seam: :func:`run` returns ``(fail_count, entries)`` where
each entry is the LOCAL dict shape::

    {"tag": "<group>.<name>", "group": <group>, "tier": "unit",
     "status": "PASS"|"FAIL", "seconds": <float>, "error": None|<str>}

:func:`entry_adapter` is the documented integration seam: today it is the
identity over that local shape; todo 27 (orchestrator) rewires it to the
suite-wide ``suite.results.ResultEntry`` contract from todo 4 without
touching the execution or printing logic. Note this local shape is NOT
the todo-27 session-log schema (which carries phases/timings) -- unit
entries are the per-check rows the orchestrator will enrich there.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import sys
import tempfile
import textwrap
import time
import traceback
from pathlib import Path

from . import collector

# ── Fancy output (ported from src/UnitTests/multiresolution/run_sph_tests.py:41-52)
BOLD_GREEN = "\033[1;32m"
BOLD_RED = "\033[1;31m"
RESET = "\033[0m"


def pass_msg(msg: str) -> None:
    print(f"{BOLD_GREEN}  ✓  {msg}{RESET}")


def fail_msg(msg: str) -> None:
    print(f"{BOLD_RED}  ✗  {msg}{RESET}")


class GroupFilterError(Exception):
    """Raised by :func:`run` when a group filter matches no check."""

    def __init__(self, groups: "list[str]", available: "list[str]") -> None:
        self.groups = list(groups)
        self.available = sorted(available)
        super().__init__(
            f"group filter {self.groups} matched no checks; "
            f"available groups: {', '.join(self.available) or '(none)'}"
        )


def entry_adapter(entry: dict) -> dict:
    """Integration seam: map a LOCAL unit entry to the suite result shape.

    Standalone today (returns the local dict unchanged). Todo 27 rewires
    this to build/attach the ``suite.results.ResultEntry`` from todo 4;
    keep every adaptation here so execution and printing stay stable.
    """
    return dict(entry)


def _unit_dir() -> Path:
    return Path(__file__).resolve().parent


def _collect_all(tests_dir: Path) -> "list[collector.Check]":
    """Discover and collect every check, fabricating named import FAILs."""
    checks: "list[collector.Check]" = []
    for path in collector.discover(tests_dir):
        try:
            checks.extend(collector.collect(path))
        except Exception:
            message = traceback.format_exc()

            def import_failure(message: str = message) -> str:
                return message

            import_failure.__name__ = f"{path.stem} (collect)"
            checks.append(
                collector.Check(
                    name=f"{path.stem} (collect)",
                    group=collector.group_of(path),
                    fn=import_failure,
                )
            )
    return checks


def run(
    groups: "list[str] | None" = None,
    verbose: bool = False,
    tests_dir: "Path | None" = None,
) -> "tuple[int, list[dict]]":
    """Execute checks and return ``(fail_count, entries)``.

    ``tests_dir`` defaults to this package's directory; inner calls (the
    collector selftests) pass a temp directory instead. Raises
    :class:`GroupFilterError` if ``groups`` is given and matches nothing.
    """
    tests_dir = Path(tests_dir) if tests_dir is not None else _unit_dir()
    checks = _collect_all(tests_dir)

    if groups:
        selected = [check for check in checks if check.group in groups]
        if not selected:
            available = {
                collector.group_of(path) for path in collector.discover(tests_dir)
            }
            raise GroupFilterError(groups, sorted(available))
        checks = selected

    fail_count = 0
    entries: "list[dict]" = []
    for check in checks:
        started = time.perf_counter()
        try:
            message = check.fn()
        except Exception:
            message = traceback.format_exc()
        seconds = time.perf_counter() - started

        status = "PASS" if message is None else "FAIL"
        if message is None:
            if verbose:
                pass_msg(check.name)
        else:
            fail_count += 1
            # Non-verbose output keeps the first line; entries keep it all.
            shown = message if verbose else message.splitlines()[0]
            fail_msg(f"{check.name} - {shown}")
        entries.append(
            entry_adapter(
                {
                    "tag": f"{check.group}.{check.name}",
                    "group": check.group,
                    "tier": "unit",
                    "status": status,
                    "seconds": round(seconds, 6),
                    "error": message,
                }
            )
        )

    if not checks:
        print("0 checks")
    else:
        print(f"{len(checks)} checks, {fail_count} failed")
    return fail_count, entries


def _available_groups(tests_dir: Path) -> "dict[str, int]":
    """Map group name -> check count; zero-check modules show up with 0."""
    counts: "dict[str, int]" = {
        collector.group_of(path): 0 for path in collector.discover(tests_dir)
    }
    for check in _collect_all(tests_dir):
        counts[check.group] = counts.get(check.group, 0) + 1
    return counts


_SELFTEST_SOURCE = textwrap.dedent(
    '''
    def check_alpha():
        return None


    def check_beta():
        return None


    def check_gamma():
        return "seeded failure: gamma"
    '''
)


def _selftest(verbose: bool = True) -> int:
    """Plan-acceptance selftest: seeded module with 2 pass + 1 fail check.

    Runs the real discover/collect/run pipeline over a temp file and exits
    1 (the seeded failure must be detected and named). The permanent suite
    member ``test_collector_selftest.py`` asserts the same through the API.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tests_dir = Path(tmp)
        (tests_dir / "test_seeded.py").write_text(_SELFTEST_SOURCE)
        fail_count, _entries = run(groups=None, verbose=verbose, tests_dir=tests_dir)
    if fail_count == 0:
        fail_msg("selftest: seeded failure was NOT detected")
        return 1
    # Seeded failure detected: exit codes are the runner's (1 = a check failed).
    return 0 if fail_count == 0 else 1


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(
        prog="unit.runner", description="TNL-SPH unit-tier collector runner."
    )
    parser.add_argument(
        "--group",
        dest="groups",
        action="append",
        default=None,
        metavar="NAME",
        help="Run only checks of this group; may be given multiple times.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Print discovered groups and check counts only.",
    )
    parser.add_argument("--verbose", action="store_true", help="Print every check line.")
    parser.add_argument(
        "--here",
        action="store_true",
        help="Discovery rooted at the unit/ package directory (the default "
        "today; named explicitly so ctest is immune to future changes in "
        "the default discovery root).",
    )
    parser.add_argument(
        "--selftest",
        action="store_true",
        help="Run a temp module with 2 passing + 1 failing check; exits 1 "
        "naming the seeded failure (harness sanity check).",
    )
    args = parser.parse_args(argv)

    if args.selftest:
        return _selftest(verbose=args.verbose or True)

    tests_dir = _unit_dir()  # --here and the default share this root today.

    if args.list:
        groups = _available_groups(tests_dir)
        if not groups:
            print("0 groups")
            return 0
        for group in sorted(groups):
            print(f"{group}: {groups[group]} check(s)")
        return 0

    try:
        fail_count, _entries = run(
            groups=args.groups, verbose=args.verbose, tests_dir=tests_dir
        )
    except GroupFilterError as error:
        # Distinct from a failing check: the selection itself was wrong.
        print(error, file=sys.stderr)
        return 2
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        sys.exit(main())
    sys.exit(130)

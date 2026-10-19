"""Selftests of the unit collector/runner - permanent member of the unit tier.

Each ``check_*`` function is itself a unit-tier check (this module's group is
``collector_selftest``) that exercises the collector and runner APIs against
temporary modules; nothing here shells out or imports pytest.

Contract asserted inline:
- ``check_*`` returns ``None`` on PASS or a string message on FAIL; raising
  also means FAIL (the runner captures the traceback as the message).
"""

from __future__ import annotations

import contextlib
import io
import tempfile
import textwrap
from pathlib import Path

from unit import collector, runner

_SEEDED_SOURCE = textwrap.dedent(
    '''
    """Seeded module: two passing checks and one failing check."""

    # Imported alias whose name starts with check_ must NOT be collected.
    from textwrap import dedent as check_imported_alias

    CONSTANT = 1


    def check_alpha():
        return None


    def check_beta():
        return None


    def check_gamma():
        return "seeded failure: gamma"
    '''
)

_EMPTY_SOURCE = '"""Empty test module: no checks here."""\n'

_CTX_SOURCE = textwrap.dedent(
    '''
    """Module with prepare() providing a shared context."""


    def prepare():
        return {"k": 42}


    def check_with_ctx(ctx):
        if ctx != {"k": 42}:
            return f"unexpected ctx: {ctx!r}"
        return None


    def check_zero_arg():
        # Zero-arg check must never receive the context (would TypeError).
        return None
    '''
)

_ENTRY_KEYS = {"tag", "group", "tier", "status", "seconds", "error"}


def _write_modules(directory: Path, modules: "dict[str, str]") -> None:
    for filename, source in modules.items():
        (directory / filename).write_text(source)


def check_selftest_pass_and_fail_discovery():
    """Collector finds the seeded checks; runner reports exactly one FAIL."""
    with tempfile.TemporaryDirectory() as tmp:
        tests_dir = Path(tmp)
        _write_modules(
            tests_dir,
            {"test_seeded.py": _SEEDED_SOURCE, "test_emptyish.py": _EMPTY_SOURCE},
        )

        # Discovery: both files found at top level, deterministically ordered.
        discovered = collector.discover(tests_dir)
        assert [path.name for path in discovered] == [
            "test_emptyish.py",
            "test_seeded.py",
        ], f"unexpected discovery order: {discovered}"

        # Collection: 3 checks in the seeded module, group from the filename.
        seeded_path = tests_dir / "test_seeded.py"
        checks = collector.collect(seeded_path)
        assert [check.name for check in checks] == [
            "check_alpha",
            "check_beta",
            "check_gamma",
        ], f"unexpected collected checks: {[check.name for check in checks]}"
        assert {check.group for check in checks} == {"seeded"}, checks

        # Empty module collects zero checks.
        assert collector.collect(tests_dir / "test_emptyish.py") == []

        # Run over the temp dir: only the seeded module contributes entries.
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            fail_count, entries = runner.run(groups=None, verbose=True, tests_dir=tests_dir)
        output = captured.getvalue()

        assert fail_count == 1, f"expected fail_count 1, got {fail_count}"
        assert len(entries) == 3, f"expected 3 entries, got {len(entries)}"
        assert any("check_gamma" in line for line in output.splitlines()), output

        for entry in entries:
            assert set(entry) == _ENTRY_KEYS, entry
            assert entry["tier"] == "unit", entry
            assert entry["group"] == "seeded", entry
        failing = [entry for entry in entries if entry["status"] == "FAIL"]
        assert len(failing) == 1, entries
        assert failing[0]["tag"] == "seeded.check_gamma", failing
        assert "seeded failure: gamma" in failing[0]["error"], failing
        passing = [entry for entry in entries if entry["status"] == "PASS"]
        assert len(passing) == 2, entries

        # Group filter selects, and an unknown group never passes silently.
        with contextlib.redirect_stdout(io.StringIO()):
            filtered_fail, filtered = runner.run(
                groups=["seeded"], tests_dir=tests_dir
            )
        assert filtered_fail == 1 and len(filtered) == 3, filtered
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                runner.run(groups=["no_such_group"], tests_dir=tests_dir)
        except runner.GroupFilterError as error:
            assert error.available == ["emptyish", "seeded"], error.available
        else:
            raise AssertionError("GroupFilterError not raised for unknown group")

    # Empty directory: "0 checks" printed, failure count 0, no entries.
    with tempfile.TemporaryDirectory() as tmp:
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            fail_count, entries = runner.run(tests_dir=Path(tmp))
        assert fail_count == 0 and entries == [], (fail_count, entries)
        assert "0 checks" in captured.getvalue(), captured.getvalue()

    return None


def check_selftest_prepare_ctx():
    """prepare() ctx reaches one-arg checks; zero-arg checks stay zero-arg."""
    with tempfile.TemporaryDirectory() as tmp:
        tests_dir = Path(tmp)
        _write_modules(tests_dir, {"test_ctxmod.py": _CTX_SOURCE})

        checks = collector.collect(tests_dir / "test_ctxmod.py")
        assert [check.name for check in checks] == [
            "check_with_ctx",
            "check_zero_arg",
        ], checks
        assert {check.group for check in checks} == {"ctxmod"}, checks

        with contextlib.redirect_stdout(io.StringIO()):
            fail_count, entries = runner.run(tests_dir=tests_dir)
        assert fail_count == 0, entries
        assert [entry["status"] for entry in entries] == ["PASS", "PASS"], entries

    return None

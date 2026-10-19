"""Unit-tier check collector for the TNL-SPH unified test suite.

Discovery contract (stdlib only, no pytest/rich):

- ``discover(tests_dir)`` lists ``test_*.py`` files **at the top level** of
  ``tests_dir`` (never recurses into ``cpp/`` or any other subdirectory).
- ``collect(path)`` loads one module via
  :func:`importlib.util.spec_from_file_location` (case ``init.py`` files are
  never executed -- discovery is suite-internal) and returns one
  :class:`Check` per top-level function whose name starts with ``check_``.
- Check signature: ``() -> str | None``. Returning ``None`` (or falling off
  the end) means PASS; returning a string means FAIL with that string as the
  failure message. Raising is treated as FAIL by the runner with the full
  traceback captured.
- The *group* of a check is the module filename stem minus the ``test_``
  prefix (``test_configurations.py`` -> group ``configurations``).
- A module may define an optional ``prepare() -> dict`` returning a shared
  context. The context is passed *as-is* to checks that accept a positional
  argument (signature inspected once at collect time); zero-argument checks
  are collected untouched. The bound closure keeps :attr:`Check.fn`
  zero-argument callable in every case.
"""

from __future__ import annotations

import importlib.util
import inspect
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any, Callable

_TEST_FILE_GLOB = "test_*.py"
_TEST_FILE_PREFIX = "test_"
_CHECK_PREFIX = "check_"
_PREPARE_NAME = "prepare"


@dataclass
class Check:
    """One executable unit check.

    ``fn`` is always a zero-argument callable returning ``str | None``;
    when the source check takes a context argument, the closure created in
    :func:`collect` captures the module's ``prepare()`` result.
    """

    name: str
    group: str
    fn: Callable[[], "str | None"]


def group_of(path: Path) -> str:
    """Return the group name derived from a test file path."""
    stem = path.stem
    if stem.startswith(_TEST_FILE_PREFIX):
        return stem[len(_TEST_FILE_PREFIX):]
    return stem


def discover(tests_dir: Path) -> list[Path]:
    """Return sorted ``test_*.py`` paths at the top level of ``tests_dir``.

    Non-recursive on purpose: ``unit/cpp/`` holds the C++ test tree, which
    is driven by CTest, not by this collector. A missing directory yields an
    empty list (treated as "0 checks" by the runner).
    """
    tests_dir = Path(tests_dir)
    if not tests_dir.is_dir():
        return []
    return sorted(p for p in tests_dir.glob(_TEST_FILE_GLOB) if p.is_file())


def _load_module(path: Path) -> ModuleType:
    """Load ``path`` as a fresh module without importing sibling packages."""
    # Registered under the canonical ``unit.<stem>`` name so relative imports
    # (``from . import collector``) inside collected modules resolve against
    # the real ``unit`` package; that package must be importable (via
    # PYTHONPATH=src/Tests or suite.bootstrap) whenever modules are collected.
    module_name = f"unit.{path.stem}"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot build import spec for {path}")
    module = importlib.util.module_from_spec(spec)
    # Register before exec so intra-module ``import unit.collected.<stem>``
    # and traceback rendering behave like a normal import.
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _accepts_positional_arg(fn: Callable[..., Any]) -> bool:
    """Return True if ``fn`` can be called with one positional argument."""
    signature = inspect.signature(fn)
    for parameter in signature.parameters.values():
        if parameter.kind in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.VAR_POSITIONAL,
        ):
            return True
    return False


def _bind_ctx(fn: Callable[..., "str | None"], ctx: Any) -> Callable[[], "str | None"]:
    """Wrap ``fn`` in a zero-argument closure passing ``ctx``."""
    def bound() -> "str | None":
        return fn(ctx)

    bound.__name__ = fn.__name__
    bound.__doc__ = fn.__doc__
    return bound


def collect(path: Path) -> list[Check]:
    """Collect all checks from one test module, sorted by check name.

    Import-time exceptions and ``prepare()`` failures propagate to the
    caller; the runner turns them into a named group-level FAIL entry so a
    broken module never aborts the tier.
    """
    path = Path(path)
    module = _load_module(path)
    group = group_of(path)

    ctx = None
    prepare = getattr(module, _PREPARE_NAME, None)
    if prepare is not None and callable(prepare):
        ctx = prepare()

    checks = []
    for name in dir(module):
        if not name.startswith(_CHECK_PREFIX):
            continue
        fn = getattr(module, name)
        if not inspect.isfunction(fn):
            continue
        # Only collect functions defined in this module (guard against
        # imported helpers that happen to be named check_*).
        if fn.__module__ != module.__name__:
            continue
        if _accepts_positional_arg(fn):
            checks.append(Check(name=name, group=group, fn=_bind_ctx(fn, ctx)))
        else:
            checks.append(Check(name=name, group=group, fn=fn))
    checks.sort(key=lambda check: check.name)
    return checks

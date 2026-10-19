"""Configuration loading and validation of the TNL-SPH unified test suite.

Loader contract (port of ``src/tools/multirunner.py:14-23`` semantics): a
configuration module is ANY python file defining a top-level
``configurations`` list; it is loaded dynamically via
:func:`importlib.util.spec_from_file_location`. Built-in group modules live in
``src/Tests/regression/configurations/<name>.py`` (todo 7).

Entry-shape rules (pinned, consumed by todos 7/9/12/21):

- ``RESERVED_KEYS`` are consumed by the suite itself and are NEVER forwarded
  to a case's init script (every other key is, see
  :func:`suite.runner_core.init_case`).
- ``evaluation-function`` (D5) is either absent/``None``/``""`` or a STRING
  naming a function in ``regression.metrics.REGISTRY`` -- raw callables in
  config data are rejected by :func:`validate_entries` (the legacy
  ``configurationsToTest.py`` keeps callables and is validated only with
  ``allow_callable_eval=True``).
- Unknown keys beyond ``RESERVED_KEYS`` are allowed unconditionally -- they are
  init-script kwargs; todo 12 checks them against the case's argparse.

Selftest (writes only into a temp dir)::

    PYTHONPATH=src/Tests python3 -m suite.config_loader
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any

from . import paths

# Suite-consumed keys, pinned exactly (downstream todos grep this list):
#   case               - example path relative to examples/ (required)
#   case-tag           - unique result tag (required)
#   evaluation-function- registry NAME of the metrics function (D5; string only)
#   init-script        - init script file name, default "init.py" (Metis G3)
#   target             - binary target override, default "<leaf-dir>_cuda" (G12)
#   init-overrides     - reserved for orchestrator-level init kwarg overrides
#   timeout_s          - reserved per-case run timeout
#   requires           - reserved per-case requirements annotation
RESERVED_KEYS = frozenset(
    {
        "case",
        "case-tag",
        "evaluation-function",
        "init-script",
        "target",
        "init-overrides",
        "timeout_s",
        "requires",
    }
)

# Directory of the built-in configuration group modules (todo 7).
BUILTIN_CONFIGURATIONS_DIR = paths.TESTS_DIR / "regression" / "configurations"


def load_configurations(conf_path: Path) -> list[dict]:
    """Load a configuration FILE and return its top-level ``configurations`` list.

    Ported from ``src/tools/multirunner.py:14-23``: the file is imported
    dynamically (its parent directory is prepended to ``sys.path`` so sibling
    data modules can be imported, e.g. an aggregate ``suite_all.py``).
    Raises ``SystemExit`` naming the path -- cleanly, no traceback -- when
    the file does not exist, cannot be loaded, or does not define a
    top-level ``configurations`` list.
    """
    conf_path = Path(conf_path)
    if not conf_path.is_file():
        raise SystemExit(f"config_loader: configuration file not found: {conf_path}")
    spec = importlib.util.spec_from_file_location("configurations_module", conf_path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not load config module from {conf_path}.")
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(conf_path.resolve().parent))
    spec.loader.exec_module(module)
    if not hasattr(module, "configurations") or not isinstance(module.configurations, list):
        raise SystemExit(
            f"Config file {conf_path} must define a top-level 'configurations' list."
        )
    return module.configurations


def list_builtin_groups() -> list[str]:
    """Return the names (file stems) of the built-in configuration groups.

    Returns an empty list when the configurations directory does not exist
    yet (it is created by todo 7) -- callers report the missing-group error.
    """
    conf_dir = BUILTIN_CONFIGURATIONS_DIR
    if not conf_dir.is_dir():
        return []
    return sorted(path.stem for path in conf_dir.glob("*.py"))


def load_builtin_group(name: str) -> list[dict]:
    """Load the built-in configuration group ``<name>.py`` by name.

    Raises ``SystemExit`` -- cleanly, no traceback -- naming the group and
    listing the available module names (``<none>`` when the configurations
    directory is absent) when the configurations directory or the module
    for ``name`` is missing.
    """
    conf_path = BUILTIN_CONFIGURATIONS_DIR / f"{name}.py"
    if not conf_path.is_file():
        available = ", ".join(list_builtin_groups()) or "<none>"
        raise SystemExit(
            f"config_loader: unknown configurations group '{name}' "
            f"(available: {available})"
        )
    return load_configurations(conf_path)


def validate_entries(
    entries: list[dict], *, allow_callable_eval: bool = False
) -> list[str]:
    """Validate loaded configuration entries; returns a list of problem strings.

    Empty list means clean. Pure function, no I/O (unit-testable).

    Per-entry checks:

    - the entry is a dict carrying ``case`` and ``case-tag``;
    - ``evaluation-function`` is MISSING, ``None``, ``""`` or a STRING registry
      name -- raw callables are rejected (D5) unless ``allow_callable_eval``
      is set (legacy data only); the entry's case-tag is named in the message;
    - unknown keys beyond ``RESERVED_KEYS`` need NO validation here: they are
      init-script kwargs checked later by todo 12's argparse comparison.
      Underscore-prefixed keys (``_group``) are internal suite plumbing and
      are likewise ignored.
    """
    problems: list[str] = []
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict):
            problems.append(f"entry[{i}]: not a dict (got {type(entry).__name__})")
            continue
        tag = entry.get("case-tag")
        if "case" not in entry:
            problems.append(
                f"entry[{i}] (case-tag {tag if tag else '<missing>'}): "
                f"missing required key 'case'"
            )
        if not tag:
            problems.append(
                f"entry[{i}] (case {entry.get('case', '<missing>')}): "
                f"missing required key 'case-tag'"
            )
        label = str(tag) if tag else f"entry[{i}]"
        if "evaluation-function" in entry:
            eval_fn: Any = entry["evaluation-function"]
            if eval_fn is None or eval_fn == "":
                continue
            if isinstance(eval_fn, str):
                continue
            if callable(eval_fn):
                if not allow_callable_eval:
                    problems.append(
                        f"entry '{label}': 'evaluation-function' must be a STRING "
                        f"registry name; raw callables are rejected (D5)"
                    )
            else:
                problems.append(
                    f"entry '{label}': 'evaluation-function' must be a string "
                    f"registry name or None; got {type(eval_fn).__name__}"
                )
    return problems


def _selftest() -> None:
    """Selftest: temp config modules through the loader; validate_entries rules."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # 1) happy path: a config module with a top-level `configurations` list.
        good = tmp_path / "good_config.py"
        good.write_text(
            "configurations = [\n"
            "    {\n"
            "        'case': 'WCSPH-DBC/damBreak2D_WCSPH-DBC',\n"
            "        'case-tag': 'WCSPH-DBC/damBreak2D_WCSPH-DBC_default',\n"
            "        'evaluation-function': 'damBreak2D_WCSPH_DBC',\n"
            "        'bc-type': 'DBC',\n"
            "        'h-coef': 1.41421356237,\n"
            "    },\n"
            "    {\n"
            "        'case': 'RSPH/damBreak2D_RSPH',\n"
            "        'case-tag': 'RSPH/damBreak2D_RSPH',\n"
            "        'evaluation-function': None,\n"
            "    },\n"
            "]\n",
            encoding="utf-8",
        )
        entries = load_configurations(good)
        assert len(entries) == 2, entries
        assert entries[0]["case-tag"] == "WCSPH-DBC/damBreak2D_WCSPH-DBC_default"
        assert validate_entries(entries) == []  # string eval + None eval accepted

        # 2) missing `configurations` -> SystemExit naming the path.
        bare = tmp_path / "bare_config.py"
        bare.write_text("X = 1\n", encoding="utf-8")
        try:
            load_configurations(bare)
        except SystemExit as exc:
            assert str(bare) in str(exc), exc
        else:
            raise AssertionError("load_configurations must raise SystemExit when 'configurations' is missing")

        # 3) wrong type of `configurations` -> SystemExit naming the path.
        wrong = tmp_path / "wrong_config.py"
        wrong.write_text("configurations = {'case': 'x'}\n", encoding="utf-8")
        try:
            load_configurations(wrong)
        except SystemExit as exc:
            assert str(wrong) in str(exc), exc
        else:
            raise AssertionError("load_configurations must raise SystemExit when 'configurations' is not a list")

        # 4) nonexistent file -> CLEAN SystemExit naming the path (no traceback).
        try:
            load_configurations(tmp_path / "does_not_exist.py")
        except SystemExit as exc:
            assert "does_not_exist" in str(exc), exc
            assert "Traceback" not in str(exc), exc
            print(f"clean failure (missing file): {exc}")
        else:
            raise AssertionError("load_configurations must raise SystemExit for a missing file")

        # 5) `configurations` that is a list of non-dicts still loads
        #    (loader does not validate) but validate_entries flags every item.
        not_dicts = tmp_path / "not_dicts.py"
        not_dicts.write_text("configurations = [42]\n", encoding="utf-8")
        problems = validate_entries(load_configurations(not_dicts))
        assert len(problems) == 1 and "not a dict" in problems[0], problems

        # 6) validate_entries: raw callable REJECTED with the case-tag named,
        #    string eval ACCEPTED, unknown init-kwargs ignored.
        def raw_callable(case_dir: Path) -> tuple[int, int]:  # noqa: ARG001
            return 0, 0

        raw = {
            "case": "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "case-tag": "WCSPH-DBC/damBreak2D_WCSPH-DBC:raw",
            "evaluation-function": raw_callable,
        }
        problems = validate_entries([raw])
        assert len(problems) == 1, problems
        assert "WCSPH-DBC/damBreak2D_WCSPH-DBC:raw" in problems[0], problems[0]
        # legacy data may still be validated when callables are allowed
        assert validate_entries([raw], allow_callable_eval=True) == []

        string_eval = dict(raw)
        string_eval["evaluation-function"] = "damBreak2D_WCSPH_DBC"
        string_eval["bc-type"] = "DBC"  # unknown key beyond RESERVED_KEYS: fine
        assert validate_entries([string_eval]) == []

        # missing case / case-tag are reported
        assert len(validate_entries([{"case-tag": "x", "evaluation-function": None}])) == 1
        assert len(validate_entries([{"case": "X", "evaluation-function": None}])) == 1
        # non-string, non-callable eval value is reported
        bad_eval = dict(string_eval)
        bad_eval["evaluation-function"] = 17
        assert len(validate_entries([bad_eval])) == 1

        # 7) RESERVED_KEYS pinned exactly (downstream grep contract).
        assert RESERVED_KEYS == frozenset(
            {
                "case",
                "case-tag",
                "evaluation-function",
                "init-script",
                "target",
                "init-overrides",
                "timeout_s",
                "requires",
            }
        ), RESERVED_KEYS

    # 8) builtin group discovery: missing dir -> [] (no crash);
    #    load_builtin_group raises a CLEAN SystemExit naming the group and
    #    listing the available names (<none> when the dir is absent).
    groups_dir = BUILTIN_CONFIGURATIONS_DIR
    groups = list_builtin_groups()
    assert isinstance(groups, list) and all(isinstance(g, str) for g in groups)
    assert groups_dir.is_dir() or groups == []
    try:
        load_builtin_group("definitely_not_a_group")
    except SystemExit as exc:
        msg = str(exc)
        assert "definitely_not_a_group" in msg, msg
        assert "Traceback" not in msg, msg
        for g in groups:
            assert g in msg, msg
        if not groups:
            assert "<none>" in msg, msg
        print(f"clean failure (unknown group): {exc}")
    else:
        raise AssertionError("load_builtin_group must raise SystemExit for an absent group")

    print("config_loader selftest PASS")


if __name__ == "__main__":
    _selftest()

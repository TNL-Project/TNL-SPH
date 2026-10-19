"""Config-rot + schema unit checks of the regression configurations (group ``configurations``).

Suite-side guards against the classic rot classes (todo 12): a configuration
entry whose case dir vanished, whose init kwargs no longer exist in the case's
init script (e.g. MR/LTS ``init_multi_resolution*.py`` variants), whose
``evaluation-function`` name drifted from the metrics REGISTRY, whose binary
target drifted from the case's ``CMakeLists.txt`` (the
``damBreak2D_WCSPH-BI_MR_cuda`` != ``<leaf-dir>_cuda`` class), or whose
reference JSON went stale.  Everything is static: configuration modules are
loaded through :func:`suite.config_loader.load_builtin_group` and init scripts
are AST-parsed for ``add_argument("--x")`` option strings -- NEVER imported or
executed, and NO builds are run.

Context built by :func:`prepare`:

- ``groups``: ``{module_name: entries}`` for every builtin configuration
  module (``list_builtin_groups()`` minus ``__init__``);
- ``records``: one derived record per unique ``case-tag`` across all groups
  (``suite_all`` duplicates its source groups on purpose and is deduped);
- ``registry``: ``regression.metrics.REGISTRY`` (11 names, ``empty`` counts);
- ``references_dir``: root of the committed reference store.

Contract: ``check_*`` returns ``None`` on PASS or a string naming the broken
rule on FAIL.  The ``check_rot_*`` members are the seeded-rot proof: they
feed doctored records/registry/staging trees into the shared workers and
assert the primary checks fail by name (acceptance scenario of todo 12).
"""

from __future__ import annotations

import ast
import re
import tempfile
from pathlib import Path

from regression.metrics import REGISTRY
from regression.metrics import references
from suite import config_loader, paths

# Reserved suite keys are never init-script kwargs (loader contract).
RESERVED_KEYS = config_loader.RESERVED_KEYS

# Cross-module duplicate case-tags are intentional ONLY across these module
# pairs (todo 12 design note): ``suite_all`` aggregates the five source
# groups, and ``qa_smoke`` re-selects two pinned tags from ``defaults`` (same
# tag -> same reference file, filtered quick runs).  Duplicates inside ONE
# module always fail; duplicates across any non-whitelisted pair fail.
_AGGREGATED_GROUPS = ("wcsph_dbc", "wcsph_bi", "rsph", "defaults", "multiresolution")
ALLOWED_DUPLICATE_PAIRS = frozenset(
    [frozenset({"suite_all", group}) for group in _AGGREGATED_GROUPS]
    + [frozenset({"qa_smoke", group}) for group in ("defaults", "suite_all")]
)

# ``references/_staging/`` is gated by design: staged refs are untracked and
# comparators ignore them.  Tolerate the directory being absent or holding
# only a placeholder marker -- anything else is rot committed by accident.
_STAGING_MARKER_NAMES = frozenset({".gitkeep", "MARKER", "README", "README.md"})


def _effective_init_script_name(entry: dict) -> str:
    """Init script file name of a config entry (Metis G3: MR/LTS override)."""
    return entry.get("init-script") or "init.py"


def _effective_target(entry: dict) -> str:
    """Binary target of a config entry (G12: default ``<leaf-dir>_cuda``)."""
    if entry.get("target"):
        return entry["target"]
    case_leaf = str(entry.get("case", "")).rstrip("/").split("/")[-1]
    return f"{case_leaf}_cuda"


def _init_args(script_path: Path) -> "set[str]":
    """Collect ``--option`` strings of every ``add_argument`` call in a script.

    Source text + :func:`ast.parse` ONLY -- the script is never imported or
    executed.  Both ``parser.add_argument(...)`` and
    ``group.add_argument(...)`` calls are matched (``g = parser.
    add_argument_group(...)`` is this repo's idiom).  Raises
    :class:`SyntaxError` for unparseable sources (the caller turns it into a
    named failure).
    """
    source = script_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(script_path))
    options: "set[str]" = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "add_argument"):
            continue
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                if arg.value.startswith("--"):
                    options.add(arg.value)
        # Also honor the keyword form: add_argument(flags=("--x", "--y")).
        for keyword in node.keywords:
            if keyword.arg != "flags":
                continue
            value = keyword.value
            if isinstance(value, (ast.List, ast.Tuple, ast.Set)):
                for element in value.elts:
                    if (
                        isinstance(element, ast.Constant)
                        and isinstance(element.value, str)
                        and element.value.startswith("--")
                    ):
                        options.add(element.value)
    return options


def _make_record(entry, index: int, module: str) -> dict:
    """Derived per-entry data used by the checks (failures are captured)."""
    if not isinstance(entry, dict):
        return {
            "tag": f"<not-a-dict entry[{index}] in {module}>",
            "entry": entry,
            "module": module,
            "case": None,
            "case_dir": None,
            "template_dir": None,
            "cmake_file": None,
            "target": None,
            "kwargs": set(),
            "init_script": None,
            "init_error": "entry is not a dict",
        }
    tag = entry.get("case-tag") or f"<missing case-tag entry[{index}] in {module}>"
    case = entry.get("case")
    case_dir = paths.EXAMPLES_DIR / case if isinstance(case, str) else None
    init_script = None
    init_error = None
    options: "set[str] | None" = None
    if case_dir is not None:
        init_script = case_dir / _effective_init_script_name(entry)
        if init_script.is_file():
            try:
                options = _init_args(init_script)
            except SyntaxError as exc:
                init_error = f"cannot parse {init_script}: {exc}"
    kwargs = {
        key
        for key in entry
        if key not in RESERVED_KEYS and not key.startswith("_")
    }
    return {
        "tag": tag,
        "entry": entry,
        "module": module,
        "case": case,
        "case_dir": case_dir,
        "template_dir": (case_dir / "template") if case_dir is not None else None,
        "cmake_file": (case_dir / "CMakeLists.txt") if case_dir is not None else None,
        "target": _effective_target(entry) if isinstance(entry, dict) else None,
        "kwargs": kwargs,
        "init_script": init_script,
        "options": options,
        "init_error": init_error,
    }


def prepare() -> dict:
    """Load every builtin configuration group and derive the check context."""
    groups = {}
    records = []
    seen_tags = set()
    for module in sorted(config_loader.list_builtin_groups()):
        if module == "__init__":
            continue
        entries = config_loader.load_builtin_group(module)
        groups[module] = entries
        for index, entry in enumerate(entries):
            tag = entry.get("case-tag") if isinstance(entry, dict) else None
            if tag is not None and tag in seen_tags:
                # suite_all duplicates its source groups, qa_smoke re-selects
                # pinned tags: one record per unique tag is enough for the
                # filesystem/registry checks (module-level dups fail in G).
                continue
            if tag is not None:
                seen_tags.add(tag)
            records.append(_make_record(entry, index, module))
    return {
        "groups": groups,
        "records": records,
        "registry": REGISTRY,
        "references_dir": references.REFERENCES_DIR,
    }


# ---------------------------------------------------------------------------
# Shared workers (also driven with doctored data by the rot checks below).


def _required_key_failures(records) -> "list[str]":
    """A: every entry has ``case``/``case-tag``; eval name stays a string."""
    problems = []
    for record in records:
        tag = record["tag"]
        entry = record["entry"]
        if not isinstance(entry, dict):
            problems.append(f"entry '{tag}': not a dict (got {type(entry).__name__})")
            continue
        if not entry.get("case"):
            problems.append(f"entry '{tag}': missing required key 'case'")
        if not entry.get("case-tag"):
            problems.append(f"entry '{tag}': missing required key 'case-tag'")
        if "evaluation-function" in entry:
            eval_fn = entry["evaluation-function"]
            if eval_fn is None or eval_fn == "":
                continue
            if not isinstance(eval_fn, str):
                # Raw callables are rejected upstream by validate_entries;
                # this asserts that assumption still holds in the data.
                problems.append(
                    f"entry '{tag}': 'evaluation-function' must be a STRING "
                    f"registry name or None; got {type(eval_fn).__name__}"
                )
    return problems


def _case_query_failures(records) -> "list[str]":
    """B: case dir + template/ + effective init script all exist."""
    problems = []
    for record in records:
        tag = record["tag"]
        case_dir = record["case_dir"]
        if case_dir is None:
            problems.append(f"entry '{tag}': no valid 'case' path (got {record['case']!r})")
            continue
        if not case_dir.is_dir():
            problems.append(f"entry '{tag}': case dir does not exist: {case_dir}")
            continue
        if not record["template_dir"].is_dir():
            problems.append(f"entry '{tag}': missing template/ dir in {case_dir}")
        if not record["init_script"].is_file():
            problems.append(
                f"entry '{tag}': init script does not exist: {record['init_script']}"
            )
    return problems


def _kwargs_failures(records) -> "list[str]":
    """C: every init kwarg is a real ``--flag`` of the entry's init script."""
    problems = []
    for record in records:
        tag = record["tag"]
        if record["init_error"] is not None:
            problems.append(f"entry '{tag}': {record['init_error']}")
            continue
        if record["options"] is None:
            continue  # missing script is reported by check_case_dirs_and_scripts
        for key in sorted(record["kwargs"]):
            flag = f"--{key}"
            if flag not in record["options"]:
                problems.append(
                    f"entry '{tag}': kwarg '{key}' is not an argparse option of "
                    f"{record['init_script']} (expected flag '{flag}')"
                )
    return problems


def _registry_failures(records, registry) -> "list[str]":
    """D: each ``evaluation-function`` string names a REGISTRY member."""
    problems = []
    for record in records:
        entry = record["entry"]
        if not isinstance(entry, dict):
            continue
        name = entry.get("evaluation-function")
        if name is None or name == "":
            continue  # absent eval == no reference comparison
        if not isinstance(name, str):
            continue  # reported by check_required_keys
        if name not in registry:
            problems.append(
                f"entry '{record['tag']}': evaluation-function '{name}' is not in "
                f"regression.metrics.REGISTRY ({len(registry)} registered names)"
            )
    return problems


def _target_failures(records) -> "list[str]":
    """E: the effective target appears as an ``add_executable`` target."""
    problems = []
    for record in records:
        tag = record["tag"]
        cmake_file = record["cmake_file"]
        target = record["target"]
        if cmake_file is None or not cmake_file.is_file():
            if record["case_dir"] is not None and record["case_dir"].is_dir():
                problems.append(f"entry '{tag}': missing CMakeLists.txt: {cmake_file}")
            continue  # missing case dir is reported by check_case_dirs_and_scripts
        cmake_text = cmake_file.read_text(encoding="utf-8")
        pattern = re.compile(
            r"add_executable\s*\(\s*" + re.escape(target) + r"(?=[\s)])"
        )
        if not pattern.search(cmake_text):
            problems.append(
                f"entry '{tag}': target '{target}' is not an add_executable target "
                f"in {cmake_file}"
            )
    return problems


def _metric_schema_problems(metric_name: str, metric) -> "list[str]":
    """One metric entry must hold numeric ``value`` + ``atol``."""
    problems = []
    if not isinstance(metric, dict):
        return [f"metric '{metric_name}': not a dict"]
    value = metric.get("value", "<missing>")
    atol = metric.get("atol", "<missing>")
    numbers = (int, float)
    valid_value = isinstance(value, numbers) or (
        isinstance(value, list) and all(isinstance(v, numbers) for v in value)
    )
    if isinstance(value, bool) or not valid_value:
        problems.append(f"metric '{metric_name}': 'value' must be a number or list of numbers")
    if isinstance(atol, bool) or not isinstance(atol, numbers):
        problems.append(f"metric '{metric_name}': 'atol' must be a number")
    if "rtol" in metric and (isinstance(metric["rtol"], bool) or not isinstance(metric["rtol"], numbers)):
        problems.append(f"metric '{metric_name}': 'rtol' must be a number when present")
    return problems


def _reference_failures(records, references_dir) -> "list[str]":
    """F: committed reference JSONs parse and hold the pinned schema fields."""
    problems = []
    for record in records:
        ref_path = references_dir / f"{references.sanitize_tag(record['tag'])}.json"
        if not ref_path.is_file():
            continue  # most tags have no frozen reference yet (H1/H2 gate)
        try:
            data = references.load_reference(ref_path)
        except (OSError, ValueError) as exc:
            problems.append(f"reference for '{record['tag']}': {exc}")
            continue
        if not isinstance(data, dict):
            problems.append(f"reference for '{record['tag']}': top level must be a dict")
            continue
        if not isinstance(data.get("case-tag"), str):
            problems.append(f"reference {ref_path}: missing/invalid 'case-tag'")
        elif data["case-tag"] != record["tag"]:
            problems.append(
                f"reference {ref_path}: 'case-tag' {data['case-tag']!r} does not "
                f"match the entry tag {record['tag']!r}"
            )
        if not isinstance(data.get("verified"), bool):
            problems.append(f"reference {ref_path}: missing/invalid 'verified' bool")
        metrics = data.get("metrics")
        if not isinstance(metrics, dict) or not metrics:
            problems.append(f"reference {ref_path}: missing/empty 'metrics' dict")
        else:
            for metric_name, metric in metrics.items():
                for problem in _metric_schema_problems(metric_name, metric):
                    problems.append(f"reference {ref_path}: {problem}")
    return problems


def _staging_failures(references_dir) -> "list[str]":
    """F': ``references/_staging/`` is absent or holds only marker files."""
    staging = references_dir / "_staging"
    if not staging.exists():
        return []
    if not staging.is_dir():
        return [f"{staging} exists and is not a directory"]
    offenders = sorted(
        path
        for path in staging.rglob("*")
        if path.is_file() and path.name not in _STAGING_MARKER_NAMES
    )
    return [
        f"staging gate: unexpected file under {staging}: {offender} "
        f"(staging refs are untracked; never commit them)"
        for offender in offenders
    ]


def _unique_tag_failures(groups) -> "list[str]":
    """G: no duplicate case-tag inside ONE module; cross-module whitelist."""
    problems = []
    tag_modules: "dict[str, list[str]]" = {}
    for module, entries in groups.items():
        seen = set()
        for index, entry in enumerate(entries):
            tag = entry.get("case-tag") if isinstance(entry, dict) else None
            if not tag:
                continue  # reported by check_required_keys
            if tag in seen:
                problems.append(
                    f"module '{module}': duplicate case-tag '{tag}' (entry[{index}]) "
                    f"- duplicates inside one module are never intentional"
                )
            seen.add(tag)
            tag_modules.setdefault(tag, []).append(module)
    for tag, modules in sorted(tag_modules.items()):
        if len(modules) <= 1:
            continue
        for left in modules:
            for right in modules:
                if left >= right:
                    continue
                if frozenset({left, right}) not in ALLOWED_DUPLICATE_PAIRS:
                    problems.append(
                        f"case-tag '{tag}' appears in both '{left}' and '{right}': "
                        f"cross-module duplicates are whitelisted to the suite_all "
                        f"aggregate pairs and qa_smoke re-selections only"
                    )
    return problems


# ---------------------------------------------------------------------------
# Primary checks (todo 12 A-G).


def check_required_keys(ctx):
    problems = _required_key_failures(ctx["records"])
    if problems:
        return f"{len(problems)} required-key problem(s):\n" + "\n".join(problems)
    return None


def check_case_dirs_and_scripts(ctx):
    problems = _case_query_failures(ctx["records"])
    if problems:
        return f"{len(problems)} case-dir/init-script problem(s):\n" + "\n".join(problems)
    return None


def check_kwargs_are_real_argparse_options(ctx):
    problems = _kwargs_failures(ctx["records"])
    if problems:
        return f"{len(problems)} unknown init kwarg(s):\n" + "\n".join(problems)
    return None


def check_registry_names(ctx):
    problems = _registry_failures(ctx["records"], ctx["registry"])
    if problems:
        return f"{len(problems)} unknown registry name(s):\n" + "\n".join(problems)
    return None


def check_targets_exist(ctx):
    problems = _target_failures(ctx["records"])
    if problems:
        return f"{len(problems)} missing build target(s):\n" + "\n".join(problems)
    return None


def check_reference_files(ctx):
    problems = _reference_failures(ctx["records"], ctx["references_dir"])
    problems += _staging_failures(ctx["references_dir"])
    if problems:
        return f"{len(problems)} reference-file problem(s):\n" + "\n".join(problems)
    return None


def check_unique_tags(ctx):
    problems = _unique_tag_failures(ctx["groups"])
    if problems:
        return f"{len(problems)} duplicate case-tag problem(s):\n" + "\n".join(problems)
    return None


# ---------------------------------------------------------------------------
# Seeded-rot proofs: the primary workers must catch doctored data by name.


def _doctored_record(record: dict, **overrides) -> dict:
    doctored = dict(record)
    doctored["entry"] = dict(record["entry"])
    doctored.update(overrides)
    return doctored


def check_rot_broken_kwarg_detected(ctx):
    """Seeded rot: a bogus init kwarg must fail naming key + case-tag."""
    record = next(r for r in ctx["records"] if r["options"] is not None)
    doctored_entry = dict(record["entry"])
    doctored_entry["definitely-not-a-flag"] = 42
    doctored = _doctored_record(record, entry=doctored_entry, kwargs={"definitely-not-a-flag"})
    problems = _kwargs_failures([doctored])
    if (
        len(problems) != 1
        or "definitely-not-a-flag" not in problems[0]
        or record["tag"] not in problems[0]
    ):
        return f"seeded kwarg rot not detected as required: {problems}"
    return None


def check_rot_unregistered_eval_detected(ctx):
    """Seeded rot: an unknown evaluation-function must name the offender."""
    record = ctx["records"][0]
    doctored_entry = dict(record["entry"])
    doctored_entry["evaluation-function"] = "no_such_metric"
    doctored = _doctored_record(record, entry=doctored_entry)
    problems = _registry_failures([doctored], ctx["registry"])
    if (
        len(problems) != 1
        or "no_such_metric" not in problems[0]
        or record["tag"] not in problems[0]
    ):
        return f"seeded registry rot not detected as required: {problems}"
    return None


def check_rot_wrong_target_detected(ctx):
    """Seeded rot: a drifted binary target must be caught via CMakeLists."""
    record = next(r for r in ctx["records"] if r["cmake_file"] and r["cmake_file"].is_file())
    doctored = _doctored_record(record, target="nonsenseTarget42_cuda")
    problems = _target_failures([doctored])
    if (
        len(problems) != 1
        or "nonsenseTarget42_cuda" not in problems[0]
        or record["tag"] not in problems[0]
    ):
        return f"seeded target rot not detected as required: {problems}"
    return None


def check_rot_duplicate_tag_detected(ctx):
    """Seeded rot: duplicate tags inside ONE module must fail by name."""
    module = next(name for name, entries in ctx["groups"].items() if entries)
    first = ctx["groups"][module][0]
    doctored_entry = dict(first)
    doctored_entry["case-tag"] = "ROT/duplicate"
    twin = dict(doctored_entry)
    problems = _unique_tag_failures({module: [doctored_entry, twin]})
    if len(problems) != 1 or module not in problems[0] or "ROT/duplicate" not in problems[0]:
        return f"seeded duplicate-tag rot not detected as required: {problems}"
    return None


def check_rot_staging_gate_detected(ctx):
    """Seeded rot: a staged reference inside _staging/ must fail the gate."""
    record = ctx["records"][0]
    with tempfile.TemporaryDirectory() as tmp:
        fake_root = Path(tmp)
        staging = fake_root / "_staging"
        staging.mkdir()
        staged = staging / f"{references.sanitize_tag(record['tag'])}.json"
        staged.parent.mkdir(parents=True)
        staged.write_text("{}", encoding="utf-8")
        problems = _staging_failures(fake_root)
        if len(problems) != 1 or str(staged) not in problems[0]:
            return f"seeded staging rot not detected as required: {problems}"
        # A pure marker placeholder stays tolerated.
        staged.unlink()
        (staging / ".gitkeep").touch()
        problems = _staging_failures(fake_root)
        if problems:
            return f"staging marker must be tolerated, got: {problems}"
    return None


def check_rot_reference_schema_detected(ctx):
    """Seeded rot: a reference JSON with metric missing 'atol' must fail."""
    record = ctx["records"][0]
    with tempfile.TemporaryDirectory() as tmp:
        fake_root = Path(tmp)
        ref_path = fake_root / f"{references.sanitize_tag(record['tag'])}.json"
        ref_path.parent.mkdir(parents=True)
        ref_path.write_text(
            (
                '{"case-tag": "' + record["tag"] + '", "verified": true,'
                ' "metrics": {"wl_avg": {"value": 0.1}}}'
            ),
            encoding="utf-8",
        )
        problems = _reference_failures([record], fake_root)
        if len(problems) != 1 or "atol" not in problems[0] or str(ref_path) not in problems[0]:
            return f"seeded reference-schema rot not detected as required: {problems}"
    return None

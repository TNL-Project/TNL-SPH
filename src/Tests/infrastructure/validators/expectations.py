"""
expectations.py

Builds the validator's *expectations* dictionary from the case's own
configuration artifacts - everything is derived, nothing is hard-coded.

Sources
-------
* ``sources/config.jsonc`` (JSON-with-comments) - domain geometry,
  resolution, particle counts, boundary loader paths.  Key spellings were
  verified against the generated
  ``examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/sources/config.jsonc``:
  ``particles.numberOfParticles`` / ``particles.numberOfBoundaryParticles``
  (NOT ``numberOfFluidParticles``), ``domain.domainOrigin-x/-y``,
  ``domain.domainSize-x/-y``, ``domain.searchRadius``, ``model.dp``,
  ``model.h``, top-level ``numberOfSubdomains``.
* The open/periodic boundary loader files referenced by the main config
  (``open-boundary-config`` / inline ``open-boundary``,
  ``periodic-boundary-config`` / inline ``periodic-boundary``) with
  ``buffer-N-...`` flat key spellings.
* The subdomains loader (``subdomains-config`` / inline ``subdomains``),
  nested ``{"subdomains": {"subdomain-N": {...}}}`` shape.

``domainSize`` semantics follow ``include/SPH/DecompositionTopology.h:60-67``:
it is the domain's MAX CORNER (end point), so grid dimensions per axis are
``ceil((size[i] - origin[i]) / searchRadius)``.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Optional


# ── JSONC loading ────────────────────────────────────────────────────────────

def strip_jsonc_comments(text: str) -> str:
    """
    Strip ``//`` line comments while respecting string literals.

    Generated TNL-SPH configs use ``//`` comments exclusively; a naive
    line-regex would corrupt URLs or ``//`` sequences inside string values,
    hence a tiny state machine.
    """
    out = []
    in_string = False
    escape = False
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if in_string:
            out.append(ch)
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_string = False
            i += 1
        else:
            if ch == '"':
                in_string = True
                out.append(ch)
                i += 1
            elif ch == "/" and i + 1 < n and text[i + 1] == "/":
                while i < n and text[i] != "\n":
                    i += 1
            else:
                out.append(ch)
                i += 1
    return "".join(out)


def parse_jsonc_text(text: str) -> Any:
    """
    Parse JSON-with-comments text.

    Handles ``//`` line comments (string-aware) and tolerates trailing
    commas before ``}`` / ``]``, which generated templates may leave
    behind after placeholder substitution.
    """
    stripped = strip_jsonc_comments(text)
    out = []
    in_string = False
    escape = False
    i = 0
    n = len(stripped)
    while i < n:
        ch = stripped[i]
        if in_string:
            out.append(ch)
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_string = False
            i += 1
        elif ch == '"':
            in_string = True
            out.append(ch)
            i += 1
        elif ch == ",":
            j = i + 1
            while j < n and stripped[j] in " \t\r\n":
                j += 1
            if j < n and stripped[j] in "}]":
                i += 1  # drop the trailing comma
            else:
                out.append(ch)
                i += 1
        else:
            out.append(ch)
            i += 1
    return json.loads("".join(out))


def load_jsonc(path: Path) -> Any:
    """Load a JSONC file from disk."""
    return parse_jsonc_text(Path(path).read_text(encoding="utf-8"))


# ── Small helpers ────────────────────────────────────────────────────────────

def _f(value: Any, default: float = 0.0) -> float:
    return float(value) if value is not None else default


def _i(value: Any, default: int = 0) -> int:
    return int(value) if value is not None else default


def _vec3(mapping: dict, stem: str, cast=float) -> list:
    """Read a ``<stem>-x/-y(/-z)`` triple; z defaults to 0 (2D configs)."""
    return [
        cast(mapping.get(f"{stem}-x", 0)),
        cast(mapping.get(f"{stem}-y", 0)),
        cast(mapping.get(f"{stem}-z", 0)),
    ]


def expected_grid_dims(origin, size, search_radius) -> list[int]:
    """
    Mirror ``DecompositionTopology.h:68``:
    ``TNL::ceil((size - origin) / radius)`` per axis, where ``size`` is the
    domain's max corner.  A degenerate axis (2D z: size == origin == 0)
    naturally yields 0.
    """
    dims = []
    for k in range(3):
        span = float(size[k]) - float(origin[k])
        if span <= 0.0:
            dims.append(0)
        else:
            dims.append(int(math.ceil(span / float(search_radius))))
    return dims


# ── Boundary-patch loaders ───────────────────────────────────────────────────

def parse_buffer_entries(obj: dict) -> list[dict]:
    """
    Normalize flat ``buffer-N-...`` entries of a boundary loader file into a
    per-patch list ordered by N.  Vector components default to 0 when the
    axis key is absent (2D templates carry only x/y).
    """
    indices = sorted({
        int(key.split("-")[1]) for key in obj
        if key.startswith("buffer-") and len(key.split("-")) > 2
        and key.split("-")[1].isdigit()
    })
    patches = []
    for n in indices:
        prefix = f"buffer-{n}-"
        patch = {
            "index": n,
            "type": obj.get(prefix + "type"),
            "identifier": obj.get(prefix + "identifier"),
            "particles": obj.get(prefix + "particles"),
            "number_of_particles": obj.get(prefix + "numberOfParticles"),
            "orientation": _vec3(obj, prefix + "orientation", float),
            "position_1": _vec3(obj, prefix + "position-1", float),
            "position_2": _vec3(obj, prefix + "position-2", float),
            "buffer_width": _vec3(obj, prefix + "bufferWidth", float),
            "paired_periodic_buffer": _i(obj.get(prefix + "paired-periodic-buffer"), -1),
            "shift_vector": _vec3(obj, prefix + "shiftVector", float),
        }
        patches.append(patch)
    return patches


def _load_loader_config(case_dir: Path, file_value: str,
                        inline_value: str) -> tuple[Optional[dict], bool, Optional[str]]:
    """
    Resolve a ``*-boundary-config`` path / inline JSON pair.

    Returns ``(parsed_object_or_None, declared, error_or_None)`` where
    *declared* marks that the main config announces a loader (non-empty
    path or non-empty inline string).
    """
    if file_value:
        path = Path(file_value)
        if not path.is_absolute():
            path = case_dir / path
        if not path.exists():
            return None, True, f"loader file not found: {path}"
        try:
            return load_jsonc(path), True, None
        except Exception as exc:  # noqa: BLE001 - surfaced to the validator
            return None, True, f"loader parse error in {path}: {exc}"
    if inline_value:
        try:
            return parse_jsonc_text(inline_value), True, None
        except Exception as exc:  # noqa: BLE001
            return None, True, f"inline loader parse error: {exc}"
    return None, False, None


def _subdomain_expectations(subdomains_obj: Optional[dict]) -> Optional[list[dict]]:
    """
    Normalize a subdomains loader object to a list with per-subdomain
    ``origin_global_coordinates`` / ``region_dimensions`` (= refinement
    factor * grid dimensions, the way DecompositionTopology.h reports it).
    Accepts the nested ``{"subdomains": {"subdomain-N": {...}}}`` shape.
    """
    if not subdomains_obj:
        return None
    entries = subdomains_obj.get("subdomains", subdomains_obj)
    result = []
    for key in sorted(entries):
        if not key.startswith("subdomain-"):
            continue
        body = entries[key]
        if not isinstance(body, dict):
            continue
        rf = _f(body.get("refinement-factor"), 1.0)
        dims = _vec3(body, "grid-dimensions", int)
        result.append({
            "index": _i(key.split("-", 1)[1]),
            "refinement_factor": rf,
            "origin_global_coordinates": _vec3(body, "origin-global-coords", int),
            "grid_dimensions": dims,
            "region_dimensions": [rf * dims[0], rf * dims[1], rf * dims[2]],
            "origin": _vec3(body, "origin", float),
            "size": _vec3(body, "size", float),
        })
    return result


# ── Main entry: expectations from a case config ──────────────────────────────

def from_config(config_jsonc: Path, *, init_overrides: dict | None = None) -> dict:
    """
    Parse the case's ``sources/config.jsonc`` and build the expectations
    dictionary consumed by ``infrastructure.validators.invariants``.

    ``init_overrides`` is shallow-merged over the extracted values (used by
    tests and by the CLI to inject rank reports / vtk counts).
    """
    config_jsonc = Path(config_jsonc).resolve()
    cfg = load_jsonc(config_jsonc)
    if not isinstance(cfg, dict):
        raise ValueError(f"config root is not an object: {config_jsonc}")

    case_dir = config_jsonc.parent.parent if config_jsonc.parent.name == "sources" \
        else config_jsonc.parent

    domain = cfg.get("domain", {})
    model = cfg.get("model", {})
    particles = cfg.get("particles", {})
    control = cfg.get("control", {})

    domain_origin = _vec3(domain, "domainOrigin", float)
    domain_size = _vec3(domain, "domainSize", float)
    search_radius = _f(domain.get("searchRadius"))

    n_sub_raw = cfg.get("numberOfSubdomains")
    number_of_subdomains = max(1, _i(n_sub_raw, 1)) if n_sub_raw is not None else 1

    open_obj, open_declared, open_error = _load_loader_config(
        case_dir, cfg.get("open-boundary-config", "") or "",
        cfg.get("open-boundary", "") or "")
    periodic_obj, periodic_declared, periodic_error = _load_loader_config(
        case_dir, cfg.get("periodic-boundary-config", "") or "",
        cfg.get("periodic-boundary", "") or "")
    subdomains_obj, subdomains_declared, subdomains_error = _load_loader_config(
        case_dir, cfg.get("subdomains-config", "") or "",
        cfg.get("subdomains", "") or "")

    open_patches = parse_buffer_entries(open_obj) if open_obj else None
    periodic_patches = parse_buffer_entries(periodic_obj) if periodic_obj else None

    periodic_count_cfg = _i(cfg.get("periodicBoundaryPatches"), 0)
    periodic_expected = max(periodic_count_cfg,
                            len(periodic_patches) if periodic_patches else 0)

    expectations = {
        "config_path": config_jsonc,
        "case_dir": case_dir,
        "domain_origin": domain_origin,
        "domain_size": domain_size,
        "search_radius": search_radius,
        "dp": model.get("dp"),
        "h": model.get("h"),
        "number_of_subdomains": number_of_subdomains,
        "fluid_n": particles.get("numberOfParticles"),
        "boundary_n": particles.get("numberOfBoundaryParticles"),
        "fluid_vtk": control.get("fluid-particles"),
        "boundary_vtk": control.get("boundary-particles"),
        "output_directory": control.get("output-directory", "results"),
        "expected_grid_dims": expected_grid_dims(domain_origin, domain_size,
                                                 search_radius),
        "open_boundary": open_patches,
        "open_boundary_declared": open_declared,
        "open_boundary_error": open_error,
        "periodic_boundary": periodic_patches,
        "periodic_boundary_declared": periodic_declared,
        "periodic_boundary_error": periodic_error,
        "periodic_patches_expected": periodic_expected,
        "subdomains": _subdomain_expectations(subdomains_obj),
        "subdomains_declared": subdomains_declared,
        "subdomains_error": subdomains_error,
        # filled by callers (CLI / tests), not derivable from the config alone
        "vtk_counts": None,
        "rank_reports": None,
        "distributed": False,
    }
    if init_overrides:
        expectations.update(init_overrides)
    return expectations


# ── VTK particle counts ──────────────────────────────────────────────────────

def _import_vtk_reader():
    """
    Import ``src/tools/readVtkParticles.py`` lazily so this module stays
    usable on vtk-less machines (counts then degrade to None + an error
    string, and the affected rules skip count-vs-vtk comparison).
    """
    try:
        import readVtkParticles  # type: ignore
        return readVtkParticles, None
    except ImportError:
        tools_dir = Path(__file__).resolve().parents[3] / "tools"
        if str(tools_dir) not in sys.path:
            sys.path.insert(0, str(tools_dir))
        try:
            import readVtkParticles  # type: ignore
            return readVtkParticles, None
        except ImportError as exc:
            return None, str(exc)


def _resolve_vtk(sources_dir: Path, declared: Optional[str], pattern: str) -> Optional[Path]:
    if declared:
        candidate = Path(declared)
        if not candidate.is_absolute():
            candidate = sources_dir.parent / candidate if declared.startswith("sources/") \
                else sources_dir / candidate
        if candidate.exists():
            return candidate
    matches = sorted(sources_dir.glob(pattern))
    return matches[0] if matches else None


def vtk_counts(sources_dir: Path, *, fluid_vtk: Optional[str] = None,
               boundary_vtk: Optional[str] = None) -> dict:
    """
    Count points in the case's fluid/boundary particle VTKs.

    File resolution order: the ``control.fluid-particles`` /
    ``control.boundary-particles`` paths from the main config, then a
    ``*fluid*.vtk`` / ``*boundary*.vtk`` glob inside ``sources_dir``.
    """
    sources_dir = Path(sources_dir)
    result: dict[str, Any] = {"fluid": None, "boundary": None,
                              "fluid_file": None, "boundary_file": None, "error": None}
    reader_mod, import_error = _import_vtk_reader()
    if reader_mod is None:
        result["error"] = f"readVtkParticles unavailable: {import_error}"
        return result

    for kind, declared, pattern in (
        ("fluid", fluid_vtk, "*fluid*.vtk"),
        ("boundary", boundary_vtk, "*boundary*.vtk"),
    ):
        path = _resolve_vtk(sources_dir, declared, pattern)
        if path is None:
            continue
        try:
            points, _point_data = reader_mod.read_vtk_particles(str(path))
        except Exception as exc:  # noqa: BLE001 - surfaced to the caller
            result["error"] = f"failed reading {path}: {exc}"
            continue
        result[kind] = int(len(points))
        result[f"{kind}_file"] = str(path)
    return result

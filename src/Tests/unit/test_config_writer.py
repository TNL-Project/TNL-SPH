"""Unit checks for ``src/tools/writeInitConfigFile.py`` (group ``config_writer``).

Each check runs the writer inside a ``tempfile.TemporaryDirectory`` whose
``template/`` + ``sources/`` layout is written inline - nothing outside the
temp dir is touched.  Placeholder spellings are copied verbatim from the
replacement tables in ``writeInitConfigFile.py:8-44`` so a rename or a
dropped table entry is caught here by name.

Contract: ``check_*`` returns ``None`` on PASS or a message naming the
broken rule on FAIL (the runner turns raises into FAILs with a traceback).
"""

from __future__ import annotations

import contextlib
import json
import os
import re
import tempfile
from pathlib import Path


def prepare():
    """Import the module under test via the suite bootstrap (adds src/tools)."""
    import suite.bootstrap  # noqa: F401  (sys.path hook for src/tools)
    import numpy as np
    import writeInitConfigFile as cf

    return {"cf": cf, "np": np}


@contextlib.contextmanager
def _cwd(path):
    """chdir into ``path``, always restoring the previous CWD."""
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _strip_jsonc_comments(text: str) -> str:
    """Remove whole-line ``//`` comments so ``json.loads`` accepts JSONC."""
    return "\n".join(
        line for line in text.splitlines() if not line.strip().startswith("//")
    )


def _write_layout(tmp, *, jsonc_template, header_template=None):
    """Create the CWD-relative ``template/`` + ``sources/`` the writer needs."""
    root = Path(tmp)
    (root / "template").mkdir()
    (root / "sources").mkdir()
    (root / "template" / "config_template.jsonc").write_text(jsonc_template)
    if header_template is not None:
        (root / "template" / "config_template.h").write_text(header_template)


# ---------------------------------------------------------------------------
# check A fixtures: placeholder spellings verbatim from `ini_replacements`
# (writeInitConfigFile.py:8-36).  `placeholderTotallyInventedThingy` has no
# setup key on purpose and must be the ONLY residue.
# ---------------------------------------------------------------------------

_INVENTED = "placeholderTotallyInventedThingy"
_PLACEHOLDER_RE = re.compile(r"placeholder[A-Za-z0-9_-]+")

_JSONC_SETUP = {
    "search_radius": 0.02,
    "domain_origin_x": 0.0,
    "domain_origin_y": 0.0,
    "domain_origin_z": 0.0,
    "domain_size_x": 1.6,
    "dp": 0.005,
    "smoothing_length": 0.01,
    "particle_mass": 0.125,
    "boundary_element_size": 0.005,
    "speed_of_sound": 40.0,
    "density": 1000.0,
    "time_step": 0.0001,
    "number_of_subdomains": 2,
    "cfl": 0.2,
    "fluid_n": 200,
    "allocated_fluid_n": 400,
}

_JSONC_TEMPLATE = """\
// synthetic template for the config_writer unit checks
{
    "domain": {
        "origin": [ placeholderDomainOrigin-x, placeholderDomainOrigin-y, placeholderDomainOrigin-z ],
        "size": [ placeholderDomainSize-x, 0.1, 0.05 ]
    },
    "searchRadius": placeholderSearchRadius,
    "dp": placeholderInitParticleDistance,
    "smoothingLength": placeholderSmoothingLength,
    "particleMass": placeholderMass,
    "boundaryElementSize": placeholderBoundaryElementSize,
    "speedOfSound": placeholderSpeedOfSound,
    "density": placeholderDensity,
    "timeStep": placeholderTimeStep,
    "numberOfSubdomains": placeholderNumberOfSubdomains,
    "cfl": placeholderCFL,
    "fluid": {
        "numberOfParticles": placeholderFluidParticles,
        "numberOfAllocatedParticles": placeholderAllocatedFluidParticles
    },
    // deliberately unset placeholder: must remain as the only residue
    "note": "placeholderTotallyInventedThingy"
}
"""

_MINIMAL_JSONC_TEMPLATE = """\
// minimal synthetic template
{
    "density": placeholderDensity
}
"""


def check_no_residues_jsonc(ctx):
    """Every SET `ini_replacements` placeholder is substituted; only the invented one remains; JSON parses."""
    cf = ctx["cf"]
    with tempfile.TemporaryDirectory() as tmp:
        _write_layout(tmp, jsonc_template=_JSONC_TEMPLATE)
        with _cwd(tmp):
            cf.write_simulation_params(dict(_JSONC_SETUP))
            emitted = Path("sources/config.jsonc").read_text()

    residues = sorted(set(_PLACEHOLDER_RE.findall(emitted)))
    if residues != [_INVENTED]:
        return (
            f"unexpected placeholder residues in emitted config.jsonc: {residues} "
            f"(allowed: {_INVENTED!r})"
        )

    parsed = json.loads(_strip_jsonc_comments(emitted))
    expected_values = {
        "searchRadius": 0.02,
        "dp": 0.005,
        "smoothingLength": 0.01,
        "particleMass": 0.125,
        "boundaryElementSize": 0.005,
        "speedOfSound": 40.0,
        "density": 1000.0,
        "timeStep": 0.0001,
        "numberOfSubdomains": 2,
        "cfl": 0.2,
    }
    for key, want in expected_values.items():
        if parsed[key] != want:
            return f"config.jsonc[{key!r}]: expected {want!r}, got {parsed[key]!r}"
    if parsed["domain"]["origin"] != [0.0, 0.0, 0.0]:
        return (
            "config.jsonc['domain']['origin']: expected [0.0, 0.0, 0.0] "
            f"(.7f zeros), got {parsed['domain']['origin']!r}"
        )
    if parsed["fluid"]["numberOfParticles"] != 200:
        return (
            "config.jsonc['fluid']['numberOfParticles']: expected 200, "
            f"got {parsed['fluid']['numberOfParticles']!r}"
        )
    if parsed["fluid"]["numberOfAllocatedParticles"] != 400:
        return (
            "config.jsonc['fluid']['numberOfAllocatedParticles']: expected 400, "
            f"got {parsed['fluid']['numberOfAllocatedParticles']!r}"
        )
    if parsed["note"] != _INVENTED:
        return f"invented placeholder must survive verbatim, got {parsed['note']!r}"
    return None


# ---------------------------------------------------------------------------
# check B fixtures: header_replacements table (writeInitConfigFile.py:38-44).
# ---------------------------------------------------------------------------

_HEADER_PLACEHOLDER_RE = re.compile(r"#placeholder[A-Za-z]+")

_HEADER_SETUP = {
    "density": 1000.0,
    "bc_type": "SyntheticBC",
    "diffusive_term": "SyntheticDiffusion",
    "viscous_term": "SyntheticViscosity",
    "time_integration": "SyntheticVerlet",
}

_HEADER_TEMPLATE = """\
// synthetic header template for the config_writer unit checks
#pragma once

using BoundaryConditionsType = #placeholderBoundaryConditionsType;
using DiffusiveTerm = #placeholderDiffusiveTerm;
using ViscousTerm = #placeholderViscosTerm;
using TimeIntegration = #placeholderTimeIntegration;
"""


def check_no_residues_header(ctx):
    """Emitted config.h: zero `#placeholder` residues with all 4 keys set; no header is emitted without a template."""
    cf = ctx["cf"]
    with tempfile.TemporaryDirectory() as tmp:
        _write_layout(tmp, jsonc_template=_MINIMAL_JSONC_TEMPLATE)
        with _cwd(tmp):
            # Without template/config_template.h no config.h may be emitted
            # (guarded by os.path.exists in write_simulation_params).
            cf.write_simulation_params(dict(_HEADER_SETUP))
            if Path("template/config.h").exists():
                return "template/config.h emitted although config_template.h is absent"
            if not Path("sources/config.jsonc").exists():
                return "sources/config.jsonc missing after write_simulation_params"

            (Path("template") / "config_template.h").write_text(_HEADER_TEMPLATE)
            cf.write_simulation_params(dict(_HEADER_SETUP))
            header_path = Path("template") / "config.h"
            if not header_path.exists():
                return "template/config.h not emitted although config_template.h exists"
            emitted = header_path.read_text()

    residues = sorted(set(_HEADER_PLACEHOLDER_RE.findall(emitted)))
    if residues:
        return f"#placeholder residues in emitted template/config.h: {residues}"
    for token in (
        "SyntheticBC",
        "SyntheticDiffusion",
        "SyntheticViscosity",
        "SyntheticVerlet",
    ):
        if token not in emitted:
            return f"substituted value {token!r} missing in emitted config.h"
    return None


def check_defaults_and_aliases(ctx):
    """`fluid_points_count` aliases `fluid_n`; allocated_* default to 1x the counts; boundary_element_size defaults to 0."""
    cf = ctx["cf"]
    template = """\
// defaults/aliases synthetic template
{
    "fluid": placeholderFluidParticles,
    "allocatedFluid": placeholderAllocatedFluidParticles,
    "allocatedBoundary": placeholderAllocatedBoundaryParticles,
    "boundaryElementSize": placeholderBoundaryElementSize
}
"""
    setup = {"fluid_points_count": 1234}
    with tempfile.TemporaryDirectory() as tmp:
        _write_layout(tmp, jsonc_template=template)
        with _cwd(tmp):
            cf.write_simulation_params(setup)
            emitted = Path("sources/config.jsonc").read_text()

    # The writer documents that the caller's dict is never polluted.
    if "fluid_n" in setup or "allocated_fluid_n" in setup or "boundary_element_size" in setup:
        return f"caller's setup dict was polluted by alias/defaults: {sorted(setup)}"

    parsed = json.loads(_strip_jsonc_comments(emitted))
    expected = {
        "fluid": 1234,
        "allocatedFluid": 1234,
        "allocatedBoundary": 0,
        "boundaryElementSize": 0,
    }
    for key, want in expected.items():
        if parsed[key] != want:
            return f"config.jsonc[{key!r}]: expected {want!r}, got {parsed[key]!r}"
    return None


# ---------------------------------------------------------------------------
# check D fixtures: write_open_boundary_params_inline (writeInitConfigFile.py:164-229).
# Both patches are x-dominant so bufferWidth may be nonzero on x only.
# ---------------------------------------------------------------------------

_OB_SETUP = {
    "open_boundary_patches": ["inlet", "outlet"],
    "inlet_n": 12,
    "outlet_n": 34,
    "inlet_orientation_x": 1.0,
    "inlet_orientation_y": 0.0,
    "inlet_orientation_z": 0.0,
    "outlet_orientation_x": 1.0,
    "outlet_orientation_y": 0.0,
    "outlet_orientation_z": 0.0,
    "inlet_position_x": 0.0,
    "inlet_position_y": 0.1,
    "inlet_position_z": -0.2,
    "inlet_size_x": 0.05,
    "inlet_size_y": 0.3,
    "inlet_size_z": 0.4,
    "outlet_position_x": 1.0,
    "outlet_position_y": 0.11,
    "outlet_position_z": -0.21,
    "outlet_size_x": 0.06,
    "outlet_size_y": 0.31,
    "outlet_size_z": 0.41,
    "inlet_width": 0.0125,
    "outlet_width": 0.025,
    "density": 1000.0,
}


def check_open_boundary_inline_math(ctx):
    """Inline OB writer: valid JSONC, position-2 == position+size per axis, bufferWidth nonzero only on the dominant axis."""
    cf = ctx["cf"]
    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp) / "ob.jsonc"
        cf.write_open_boundary_params_inline(dict(_OB_SETUP), output_path=str(output))
        emitted = output.read_text()

    parsed = json.loads(_strip_jsonc_comments(emitted))

    if parsed["buffer-1-type"] != "inlet":
        return f"buffer-1-type: expected 'inlet', got {parsed['buffer-1-type']!r}"
    if parsed["buffer-2-type"] != "outlet":
        return f"buffer-2-type: expected 'outlet', got {parsed['buffer-2-type']!r}"

    for i, patch in enumerate(_OB_SETUP["open_boundary_patches"], start=1):
        for ax in "xyz":
            want = _OB_SETUP[f"{patch}_position_{ax}"] + _OB_SETUP[f"{patch}_size_{ax}"]
            got = parsed[f"buffer-{i}-position-2-{ax}"]
            if abs(got - want) > 1e-6:
                return (
                    f"buffer-{i}-position-2-{ax}: expected position+size = {want}, "
                    f"got {got}"
                )
            want1 = _OB_SETUP[f"{patch}_position_{ax}"]
            got1 = parsed[f"buffer-{i}-position-1-{ax}"]
            if abs(got1 - want1) > 1e-6:
                return (
                    f"buffer-{i}-position-1-{ax}: expected position = {want1}, "
                    f"got {got1}"
                )

        for ax in "xyz":
            want = _OB_SETUP[f"{patch}_width"] if ax == "x" else 0.0
            got = parsed[f"buffer-{i}-bufferWidth-{ax}"]
            if abs(got - want) > 1e-6:
                return (
                    f"buffer-{i}-bufferWidth-{ax}: dominant orientation axis is x; "
                    f"expected {want}, got {got}"
                )

        got_density = parsed[f"buffer-{i}-density"]
        if abs(got_density - _OB_SETUP["density"]) > 1e-6:
            return (
                f"buffer-{i}-density: expected setup['density'] = "
                f"{_OB_SETUP['density']}, got {got_density}"
            )
    return None


def check_save_params_numpy(ctx):
    """save_params_to_json converts numpy scalars/arrays into JSON-native values."""
    cf, np = ctx["cf"], ctx["np"]
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "params.json"
        cf.save_params_to_json(
            {"a": np.int64(3), "b": np.float32(1.5), "c": np.array([1, 2])},
            str(target),
        )
        parsed = json.loads(target.read_text())

    if parsed["a"] != 3 or not isinstance(parsed["a"], int):
        return f"'a' must round-trip as int 3, got {parsed['a']!r}"
    if abs(parsed["b"] - 1.5) > 1e-6:
        return f"'b' must round-trip as ~1.5, got {parsed['b']!r}"
    if parsed["c"] != [1, 2]:
        return f"'c' must round-trip as [1, 2], got {parsed['c']!r}"
    return None

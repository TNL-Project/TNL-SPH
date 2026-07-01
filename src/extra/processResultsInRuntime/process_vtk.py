#!/usr/bin/env python3
"""
Offscreen ParaView rendering script.

Usage (new — multi-source):
    pvpython process_vtk.py <state_file> <output_png> [--template TEMPLATE]
                            [--project-dir DIR] <vtk_file> [<vtk_file> ...]

Usage (old — single VTK, backward compatible via watch_and_render_simple.py):
    pvpython process_vtk.py <state_file> <output_png> <vtk_file>

The --template flag controls how particle-set filenames are parsed and matched
against paths embedded in the .pvsm state file.  Placeholders:

    {id}   — particle-set identifier (e.g. '0' in 'fluid_subdomain0_...')
    {time} — timestep (e.g. '0.001000')

Default template:  fluid_{time}_particles.vtk   (single-resolution case)

For multiresolution cases with multiple particle sets per state, use e.g.:
    --template fluid_subdomain{id}_{time}_particles.vtk

When {id} is present, each incoming VTK file is matched to the corresponding
.pvsm reference by its id.  When absent, all fluid references are replaced
with the single incoming file (legacy 2D behaviour).
"""
from __future__ import annotations

import sys
import os
import re
import shutil
import tempfile
import argparse
from pathlib import Path


# ---------------------------------------------------------------------------
# Template → regex conversion
# ---------------------------------------------------------------------------

def template_to_regex(template: str) -> re.Pattern:
    """Convert a filename template to a compiled regex with named groups.

    Template placeholders:
        {time} → (?P<time>[\\d.]+)        — digits and dots
        {id}   → (?P<id>[^_]+)            — anything except underscore
        {name} → (?P<name>[^/]+)          — generic placeholder
    """
    parts = re.split(r"(\{\w+\})", template)
    regex_parts: list[str] = []
    for part in parts:
        if not part:
            continue
        if part.startswith("{") and part.endswith("}"):
            name = part[1:-1]
            if name == "time":
                regex_parts.append(r"(?P<time>[\d.]+)")
            elif name == "id":
                regex_parts.append(r"(?P<id>[^_]+)")
            else:
                regex_parts.append(f"(?P<{name}>[^/]+)")
        else:
            regex_parts.append(re.escape(part))
    return re.compile("".join(regex_parts))


# ---------------------------------------------------------------------------
# .pvsm patching
# ---------------------------------------------------------------------------

def patch_pvsm(
    state_file: str,
    vtk_paths: list[str],
    template: str,
    tmp_dir: str,
    project_dir: str | None = None,
) -> str:
    """Replace particle-set paths in the .pvsm XML with the real VTK files.

    Returns the path to a patched copy in *tmp_dir* (original is untouched).
    """
    with open(state_file, "r") as f:
        content = f.read()

    tmpl_regex = template_to_regex(template)

    # For .pvsm matching we allow an arbitrary path prefix that does not
    # cross XML attribute boundaries (quotes, angle brackets).
    pvsm_regex = re.compile(r"[^\"'<>]*" + tmpl_regex.pattern)

    has_id = "id" in tmpl_regex.groupindex

    if has_id and len(vtk_paths) > 1:
        # --- Multi-source: match each .pvsm reference to its VTK by id ---
        id_to_path: dict[str, str] = {}
        for vp in vtk_paths:
            m = tmpl_regex.search(os.path.basename(vp))
            if m:
                fid = m.group("id")
                id_to_path[fid] = os.path.abspath(vp)
            else:
                print(f"[WARN] VTK file does not match template: {vp}")

        def _replacer(m: re.Match) -> str:
            fname_match = tmpl_regex.search(m.group(0))
            if fname_match:
                fid = fname_match.group("id")
                if fid in id_to_path:
                    return id_to_path[fid]
            return m.group(0)  # no matching VTK — leave unchanged

        content = pvsm_regex.sub(_replacer, content)

    else:
        # --- Single-source (2D legacy): replace all fluid refs with one file ---
        vtk_abspath = os.path.abspath(vtk_paths[0])
        content = pvsm_regex.sub(lambda m: vtk_abspath, content)

    # Resolve {project_dir} placeholder so static (e.g. boundary) paths work
    if project_dir:
        content = content.replace("{project_dir}", str(project_dir))

    tmp_state = os.path.join(tmp_dir, os.path.basename(state_file))
    with open(tmp_state, "w") as f:
        f.write(content)
    return tmp_state


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render(
    vtk_files: list[str],
    state_file: str,
    output_png: str,
    template: str,
    project_dir: str | None = None,
):
    from paraview.simple import (  # type: ignore[import-not-found]
        LoadState,
        GetActiveView,
        Render,
        SaveScreenshot,
        GetRenderView,
    )

    tmp_dir = tempfile.mkdtemp()
    try:
        patched_state = patch_pvsm(
            state_file, vtk_files, template, tmp_dir, project_dir
        )

        LoadState(patched_state)

        view = GetRenderView()
        Render(view)

        SaveScreenshot(
            output_png,
            view,
            ImageResolution=[1920, 1080],
            OverrideColorPalette="WhiteBackground",
        )
        print(f"[OK] Saved: {output_png}")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Offscreen ParaView rendering with multi-source support.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("state_file", help="Path to .pvsm state file")
    parser.add_argument("output_png", help="Path to output PNG screenshot")
    parser.add_argument(
        "vtk_files",
        nargs="+",
        help="One or more VTK files to render (order does not matter; "
        "matching is done via the template)",
    )
    parser.add_argument(
        "--template",
        default="fluid_{time}_particles.vtk",
        help="Filename template with {id} and {time} placeholders "
        '(default: "fluid_{time}_particles.vtk")',
    )
    parser.add_argument(
        "--project-dir",
        default=None,
        help="Project root to resolve {project_dir} placeholders in .pvsm",
    )
    args = parser.parse_args()

    if not os.path.exists(args.state_file):
        print(f"[ERROR] State file not found: {args.state_file}")
        sys.exit(1)
    for vf in args.vtk_files:
        if not os.path.exists(vf):
            print(f"[ERROR] VTK file not found: {vf}")
            sys.exit(1)

    render(
        args.vtk_files,
        args.state_file,
        args.output_png,
        args.template,
        args.project_dir,
    )


if __name__ == "__main__":
    main()

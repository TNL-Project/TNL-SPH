#!/usr/bin/env python3
"""
Offscreen ParaView rendering script.
Usage: pvpython process_vtk.py <vtk_file> <state_file> <output_png>
"""
import sys
import os
import re
import shutil
import tempfile

def patch_pvsm(state_file, vtk_path, tmp_dir):
    """Replace the hardcoded VTK path in the .pvsm XML with the real one."""
    with open(state_file, 'r') as f:
        content = f.read()

    # Replace any absolute path ending in fluid_*_particles.vtk
    patched = re.sub(
        r'[^"\'<>]+fluid_[\d.]+_particles\.vtk',
        vtk_path,
        content
    )

    tmp_state = os.path.join(tmp_dir, os.path.basename(state_file))
    with open(tmp_state, 'w') as f:
        f.write(patched)
    return tmp_state


def render(vtk_file, state_file, output_png):
    from paraview.simple import (
        LoadState, GetActiveView, Render,
        SaveScreenshot, GetRenderView
    )
    import paraview.simple as pvs

    # Work with a patched copy of the state so we don't modify the original
    tmp_dir = tempfile.mkdtemp()
    try:
        patched_state = patch_pvsm(state_file, os.path.abspath(vtk_file), tmp_dir)

        # Load the patched state
        LoadState(patched_state)

        view = GetRenderView()
        Render(view)

        SaveScreenshot(output_png, view,
                       ImageResolution=[1920, 1080],
                       OverrideColorPalette='WhiteBackground')
        print(f"[OK] Saved: {output_png}")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: pvpython process_vtk.py <vtk_file> <state_file> <output_png>")
        sys.exit(1)

    vtk_file, state_file, output_png = sys.argv[1], sys.argv[2], sys.argv[3]

    if not os.path.exists(vtk_file):
        print(f"[ERROR] VTK file not found: {vtk_file}")
        sys.exit(1)
    if not os.path.exists(state_file):
        print(f"[ERROR] State file not found: {state_file}")
        sys.exit(1)

    render(vtk_file, state_file, output_png)

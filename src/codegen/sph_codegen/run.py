#!/usr/bin/env python3
"""
run.py  –  generate C++ from a scheme file.
Usage:  python run.py examples/scheme_clean.py --output ./out
"""
import sys, argparse, importlib.util
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

def load(path):
    spec = importlib.util.spec_from_file_location("scheme", path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def main():
    p = argparse.ArgumentParser()
    p.add_argument("scheme")
    p.add_argument("--output", "-o", default="./generated")
    args = p.parse_args()

    mod = load(args.scheme)
    required = ["Fluid","Boundary","Constants","Interactions","Integration"]
    for r in required:
        if not hasattr(mod, r):
            print(f"ERROR: scheme file must define class {r}")
            sys.exit(1)

    from sph_codegen.gen import generate_all
    print(f"Generating from '{args.scheme}' → '{args.output}/'")
    generate_all(
        fluid_cls        = mod.Fluid,
        boundary_cls     = mod.Boundary,
        constants_cls    = mod.Constants,
        interactions_cls = mod.Interactions,
        integration_cls  = mod.Integration,
        output_dir       = args.output,
        namespace        = "WCSPH_DBC",
    )
    print("Done.")

if __name__ == "__main__":
    main()

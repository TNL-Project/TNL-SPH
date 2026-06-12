#! /usr/bin/env python3
"""
parse_sph_log.py

Parses SPH simulation log files and extracts the multiresolution section
(Fluid/Boundary object info + Multiresolution boundary buffers) into JSON.

Usage:
    python parse_sph_log.py <logfile> [--output <outfile.json>]
"""

import re
import json
import argparse
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Low-level helpers
# ─────────────────────────────────────────────────────────────────────────────

def _parse_value(raw: str):
    """Convert a raw string value to int, float, or keep as string."""
    raw = raw.strip()
    # Vector: [ x, y ] or [ x, y, z ]
    vec_match = re.match(r'^\[\s*([-\d.e+]+)\s*,\s*([-\d.e+]+)(?:\s*,\s*([-\d.e+]+))?\s*\]$', raw)
    if vec_match:
        parts = [g for g in vec_match.groups() if g is not None]
        return [_scalar(p) for p in parts]
    return _scalar(raw)


def _scalar(s: str):
    s = s.strip()
    try:
        i = int(s)
        return i
    except ValueError:
        pass
    try:
        f = float(s)
        return f
    except ValueError:
        pass
    return s


def _parse_table_line(line: str):
    """
    Parse a log table row like:
      | Key name:   value |
    Returns (key_snake_case, value) or None if not a data row.
    """
    # Strip outer pipes and whitespace
    inner = line.strip()
    if not (inner.startswith('|') and inner.endswith('|')):
        return None
    inner = inner[1:-1]

    # Must contain a colon (key: value pattern)
    colon = inner.find(':')
    if colon == -1:
        return None

    key_raw = inner[:colon].strip()
    val_raw = inner[colon + 1:].strip()

    if not key_raw or not val_raw:
        return None

    # Skip pure separator/title rows (no value after colon or all dashes)
    if re.match(r'^[-+]+$', val_raw):
        return None

    key = _to_snake(key_raw)
    value = _parse_value(val_raw)
    return key, value


def _to_snake(s: str) -> str:
    """Convert 'Some Key (with parens)' → 'some_key_with_parens'."""
    s = re.sub(r'[()\/]', ' ', s)
    s = re.sub(r'[^a-zA-Z0-9 _]', '', s)
    s = re.sub(r'\s+', '_', s.strip())
    return s.lower()


# ─────────────────────────────────────────────────────────────────────────────
# Section detection
# ─────────────────────────────────────────────────────────────────────────────

SECTION_PATTERNS = [
    (re.compile(r'Fluid\s+(\d+)\s+object\s+information', re.I),      'fluid'),
    (re.compile(r'Boundary\s+(\d+)\s+object\s+information', re.I),   'boundary'),
    (re.compile(r'Multiresolu[t]?in\s+boundary\s+buffer(\d+)', re.I),'buffer'),
]


def _detect_section(line: str):
    """Return (type, index) or None."""
    inner = re.sub(r'[|+\-]', '', line).strip()
    for pattern, kind in SECTION_PATTERNS:
        m = pattern.search(inner)
        if m:
            return kind, int(m.group(1))
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Main parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_log(log_text: str) -> dict:
    """
    Parse the relevant section of an SPH log and return a structured dict:
    {
      "fluids":    { "0": {...}, "1": {...} },
      "boundaries": { "0": {...}, "1": {...} },
      "buffers":   { "1": {...}, "2": {...} }
    }
    """
    result = {
        "fluids": {},
        "boundaries": {},
        "buffers": {},
    }

    lines = log_text.splitlines()
    current_section = None   # (kind, index)
    current_data = {}
    current_subsection = None   # tracks sub-tables inside buffers

    def flush():
        if current_section is None or not current_data:
            return
        kind, idx = current_section
        key = str(idx)
        container = {
            'fluid': result['fluids'],
            'boundary': result['boundaries'],
            'buffer': result['buffers'],
        }[kind]
        if key not in container:
            container[key] = {}
        container[key].update(current_data)

    for line in lines:
        # ── Check for new section header ──────────────────────────────────
        detected = _detect_section(line)
        if detected:
            flush()
            current_section = detected
            current_data = {}
            current_subsection = None
            continue

        if current_section is None:
            continue

        # ── Separator line: reset subsection grouping ─────────────────────
        if re.match(r'^\s*\+[-]+\+\s*$', line):
            current_subsection = None
            continue

        # ── Empty title row (all spaces inside pipes) ─────────────────────
        if re.match(r'^\s*\|\s+\|\s*$', line):
            continue

        # ── Possible sub-section title (no colon, bold text in pipes) ─────
        title_match = re.match(r'^\s*\|\s{1,4}([A-Za-z][^|:]{3,}?)\s*\|\s*$', line)
        if title_match and ':' not in title_match.group(1):
            current_subsection = _to_snake(title_match.group(1).strip())
            if current_subsection not in current_data:
                current_data[current_subsection] = {}
            continue

        # ── Data row ──────────────────────────────────────────────────────
        parsed = _parse_table_line(line)
        if parsed is None:
            continue
        k, v = parsed

        if current_subsection and current_subsection in current_data and isinstance(current_data[current_subsection], dict):
            current_data[current_subsection][k] = v
        else:
            current_data[k] = v

    flush()
    return result


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description='Parse SPH log → JSON')
    ap.add_argument('logfile', help='Path to the .log file')
    ap.add_argument('--output', '-o', help='Output JSON path (default: <logfile>.json)')
    args = ap.parse_args()

    log_path = Path(args.logfile)
    out_path = Path(args.output) if args.output else log_path.with_suffix('.json')

    text = log_path.read_text(encoding='utf-8')
    data = parse_log(text)

    out_path.write_text(json.dumps(data, indent=2), encoding='utf-8')
    print(f'Parsed → {out_path}')
    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    main()

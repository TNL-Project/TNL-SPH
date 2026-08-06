#!/usr/bin/env python3
"""Render TNL-SPH benchmark results from a results folder as tables or charts.

Reads ``summary.json`` produced by ``benchmark_particles.py`` and renders
comparison tables (terminal / Markdown / LaTeX) or grouped bar charts
(matplotlib PNG).

Examples
--------
  # Render table in terminal
  ./plot_results.py benchmark_results/20250806_143022

  # Export as Markdown
  ./plot_results.py benchmark_results/20250806_143022 --format markdown -o table.md

  # Plot grouped bar charts for search + interaction stages
  ./plot_results.py benchmark_results/20250806_143022 --plot

  # Plot specific metrics and save as PNG
  ./plot_results.py benchmark_results/20250806_143022 --plot \
      --plot-metrics search,interact -o chart.png

  # Pivot table, compare search time
  ./plot_results.py benchmark_results/20250806_143022 --view pivot --metric search
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

import numpy as np

METRICS = ["total_time", "search", "interact", "integrate", "run_seconds"]
METRIC_LABELS = {
    "total_time": "Total [s]",
    "search": "Search [s]",
    "interact": "Interact [s]",
    "integrate": "Integrate [s]",
    "run_seconds": "Wall [s]",
}


def load_results(folder: Path) -> list[dict]:
    summary = folder / "summary.json"
    if not summary.is_file():
        sys.exit(f"error: {summary} not found. Pass the results timestamp folder.")
    return json.loads(summary.read_text())


def _fmt(v) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.3f}"
    return str(v)


def _speedup_str(base: Optional[float], val: Optional[float]) -> str:
    if base is None or val is None or val == 0:
        return "—"
    s = base / val
    return f"{s:.2f}x"


def compute_baselines(results: list[dict]) -> dict[float, float]:
    """For each resolution find the baseline total_time.

    Baseline = first variant in the results list at that resolution.
    """
    baselines: dict[float, float] = {}
    for r in results:
        res = r["resolution"]
        if res not in baselines and r.get("total_time") is not None:
            baselines[res] = r["total_time"]
    return baselines


# --------------------------------------------------------------------------- #
# Flat table (one row per run)
# --------------------------------------------------------------------------- #

FLAT_HEADERS = ["Variant", "dp", "Particles", "Total [s]", "Search [s]",
                "Interact [s]", "Integrate [s]", "Steps", "Speedup", "Status"]


def flat_rows(results: list[dict]) -> list[list[str]]:
    baselines = compute_baselines(results)
    rows = []
    for r in results:
        res = r["resolution"]
        base = baselines.get(res)
        tot = r.get("total_time")
        rows.append([
            r["variant"],
            f"{res:g}",
            f"{r.get('fluid_particles', 0):,}" if r.get("fluid_particles") else "—",
            _fmt(tot),
            _fmt(r.get("search")),
            _fmt(r.get("interact")),
            _fmt(r.get("integrate")),
            str(r.get("simulation_steps", "")) or "—",
            _speedup_str(base, tot),
            "ok" if r.get("ok") else f"FAIL",
        ])
    return rows


# --------------------------------------------------------------------------- #
# Pivot table (resolutions as rows, variants as columns)
# --------------------------------------------------------------------------- #

def pivot_data(results: list[dict], metric: str):
    """Organise results as {resolution: {variant: value}}.

    Returns (resolutions_sorted, variants_ordered, data_dict).
    """
    variants: list[str] = []
    data: dict[float, dict[str, Optional[float]]] = {}
    for r in results:
        v = r["variant"]
        if v not in variants:
            variants.append(v)
        res = r["resolution"]
        if res not in data:
            data[res] = {}
        data[res][v] = r.get(metric) if r.get("ok") else None
    resolutions = sorted(data.keys())
    return resolutions, variants, data


def pivot_rows(results: list[dict], metric: str) -> tuple[list[str], list[list[str]]]:
    resolutions, variants, data = pivot_data(results, metric)
    label = METRIC_LABELS.get(metric, metric)

    headers = ["dp", "Particles"] + [f"{v}\n{label}" for v in variants] + \
              [f"{v}\nSpeedup" for v in variants]

    # Find baseline variant (first one) per resolution
    baseline_variant: Optional[str] = variants[0] if variants else None

    rows = []
    for res in resolutions:
        row_data = data[res]
        # particle count from first available variant at this resolution
        particles = "—"
        for r in results:
            if r["resolution"] == res and r.get("fluid_particles"):
                particles = f"{r['fluid_particles']:,}"
                break

        row = [f"{res:g}", particles]
        # metric values
        base_val = row_data.get(baseline_variant) if baseline_variant else None
        for v in variants:
            row.append(_fmt(row_data.get(v)))
        # speedups
        for v in variants:
            val = row_data.get(v)
            row.append(_speedup_str(base_val, val) if v != baseline_variant else "1.00x")
        rows.append(row)

    return headers, rows


# --------------------------------------------------------------------------- #
# Renderers
# --------------------------------------------------------------------------- #

def render_terminal_rich(headers, rows, title=""):
    from rich.console import Console
    from rich.table import Table
    table = Table(title=title or "TNL-SPH benchmark results", show_lines=True)
    for h in headers:
        table.add_column(h)
    for row in rows:
        table.add_row(*[str(c) for c in row])
    Console().print(table)


def render_terminal_plain(headers, rows, title=""):
    print(title or "TNL-SPH benchmark results")
    print()
    str_rows = [[str(c) for c in row] for row in rows]
    widths = [max([len(h)] + [len(r[i]) for r in str_rows])
              for i, h in enumerate(headers)]
    sep = "+".join("-" * (w + 2) for w in widths)
    print(sep)
    print(" | ".join(h.ljust(w) for h, w in zip(headers, widths)))
    print(sep.replace("-", "="))
    for row in str_rows:
        print(" | ".join(c.ljust(w) for c, w in zip(row, widths)))
    print(sep)


def render_terminal(headers, rows, title=""):
    try:
        render_terminal_rich(headers, rows, title)
    except ImportError:
        render_terminal_plain(headers, rows, title)


def render_markdown(headers, rows) -> str:
    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        # newlines in headers (pivot) → replace with <br> for Markdown
        cells = [str(c).replace("\n", "<br>") for c in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def render_latex(headers, rows, title="") -> str:
    # Sanitise headers for LaTeX
    def tex(s):
        return str(s).replace("_", r"\_").replace("\n", r"\\").replace("#", r"\#")

    ncols = len(headers)
    spec = "l" * ncols
    lines = []
    if title:
        lines.append(r"\begin{table}[htbp]")
        lines.append(r"\centering")
        lines.append(r"\caption{" + tex(title) + "}")
    lines.append(r"\begin{tabular}{" + spec + "}")
    lines.append(r"\hline")
    lines.append(" & ".join(tex(h) for h in headers) + r" \\")
    lines.append(r"\hline")
    for row in rows:
        lines.append(" & ".join(tex(c) for c in row) + r" \\")
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    if title:
        lines.append(r"\end{table}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Bar chart renderer (matplotlib)
# --------------------------------------------------------------------------- #

PLOT_METRIC_LABELS = {
    "search": "Neighbor search",
    "interact": "Interaction",
    "integrate": "Integration",
    "total_time": "Total",
    "run_seconds": "Wall time",
}


def render_plot(results: list[dict], metrics: list[str], output: Optional[str],
                title_prefix: str = "TNL-SPH benchmark") -> None:
    """Generate grouped bar charts comparing stages across variants.

    One subplot per metric.  Within each subplot, bars are grouped by
    resolution with one bar per variant.  Speedup vs baseline is annotated
    above each non-baseline bar.
    """
    import matplotlib
    if output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    resolutions, variants, _ = pivot_data(results, metrics[0])
    n_metrics = len(metrics)
    n_res = len(resolutions)
    n_var = len(variants)

    fig, axes = plt.subplots(1, n_metrics, figsize=(6 * n_metrics, 5),
                             constrained_layout=True)
    if n_metrics == 1:
        axes = [axes]

    colors = plt.cm.tab10(np.linspace(0, 0.9, max(n_var, 1)))
    bar_width = 0.8 / max(n_var, 1)
    x = np.arange(n_res)

    for ax, metric in zip(axes, metrics):
        _, _, data = pivot_data(results, metric)
        baseline_variant = variants[0]

        for vi, variant in enumerate(variants):
            vals = []
            base_vals = []
            for res in resolutions:
                v = data[res].get(variant)
                vals.append(v if v is not None else 0)
                bv = data[res].get(baseline_variant)
                base_vals.append(bv if bv is not None else 0)

            offset = (vi - (n_var - 1) / 2) * bar_width
            bars = ax.bar(x + offset, vals, bar_width * 0.9,
                          label=variant, color=colors[vi],
                          edgecolor="white", linewidth=0.5)

            for bi, (bar, val, base_val) in enumerate(zip(bars, vals, base_vals)):
                if val == 0:
                    continue
                if variant != baseline_variant and base_val > 0:
                    speedup = base_val / val
                    ax.text(bar.get_x() + bar.get_width() / 2, val,
                            f"{speedup:.2f}x", ha="center", va="bottom",
                            fontsize=7, fontweight="bold",
                            color=colors[vi])

        ax.set_xticks(x)
        ax.set_xticklabels([f"dp={r:g}" for r in resolutions])
        ax.set_ylabel("Time [s]")
        ax.set_title(PLOT_METRIC_LABELS.get(metric, metric))
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

    fig.suptitle(f"{title_prefix} — {', '.join(variants)}",
                 fontsize=12, fontweight="bold")

    if output:
        fig.savefig(output, dpi=150)
        print(f"Chart saved to {output}")
    else:
        plt.show()


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Render TNL-SPH benchmark results as tables.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("folder", help="path to the results folder (contains summary.json)")
    ap.add_argument("--format", choices=["terminal", "markdown", "latex"],
                    default="terminal", help="output format (default: terminal)")
    ap.add_argument("--view", choices=["flat", "pivot"], default="pivot",
                    help="table layout: flat (one row per run) or pivot "
                         "(resolutions×variants). Default: pivot")
    ap.add_argument("--metric", choices=METRICS, default="search,interact,total_time",
                    help="metric to compare in pivot view (default: total_time)")
    ap.add_argument("--plot", action="store_true",
                    help="generate grouped bar charts instead of a table")
    ap.add_argument("--plot-metrics", default="search,interact",
                    help="comma-separated metrics to plot (default: search,interact). "
                         f"Choices: {', '.join(METRICS)}")
    ap.add_argument("-o", "--output", default=None,
                    help="write to file instead of stdout/terminal "
                         "(tables: .md/.tex; charts: .png)")
    args = ap.parse_args()

    folder = Path(args.folder)
    if not folder.is_dir():
        sys.exit(f"error: not a directory: {folder}")

    results = load_results(folder)
    if not results:
        sys.exit("error: no results found in summary.json")

    # Chart mode
    if args.plot:
        plot_metrics = [m.strip() for m in args.plot_metrics.split(",") if m.strip()]
        for m in plot_metrics:
            if m not in METRICS:
                sys.exit(f"error: unknown metric '{m}'. Choices: {', '.join(METRICS)}")
        render_plot(results, plot_metrics, args.output)
        return 0

    # Table mode
    if args.view == "flat":
        headers, rows = FLAT_HEADERS, flat_rows(results)
        title = ""
    else:
        headers, rows = pivot_rows(results, args.metric)
        metric_label = METRIC_LABELS.get(args.metric, args.metric)
        title = f"TNL-SPH benchmark — {metric_label} (speedup vs {results[0]['variant']})"

    # Render
    if args.format == "terminal":
        render_terminal(headers, rows, title)
    else:
        if args.format == "markdown":
            output = render_markdown(headers, rows)
        else:
            output = render_latex(headers, rows, title)

        if args.output:
            Path(args.output).write_text(output)
            print(f"Written to {args.output}")
        else:
            sys.stdout.write(output)

    return 0


if __name__ == "__main__":
    sys.exit(main())

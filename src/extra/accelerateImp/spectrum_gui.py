"""Interactive matplotlib explorer for the WCSPH stability spectrum."""
from __future__ import annotations

import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Button, RadioButtons, Slider

from spectrum import (
    KERNELS,
    PlotConfig,
    Parameters,
    compute_all_direct,
    compute_all_general,
    plot_eigvals_panel,
    plot_eigvals_panel_random,
    plot_mu_panel,
    plot_mu_panel_random,
)

# --- shared state ----------------------------------------------------------
cfg = PlotConfig()
# Smaller patch / fewer realisations than the static script so each slider
# update stays sub-second. h is set explicitly (independent of H).
params = Parameters(n_cell=8, n_real=3, h=0.006)

# --- figure + plot axes ----------------------------------------------------
fig = plt.figure(figsize=(15, 9))
fig.suptitle("Spectrum explorer", fontsize=cfg.title_fontsize)

gs = GridSpec(
    2, 3,
    top=0.90, bottom=0.30, left=0.06, right=0.97,
    hspace=0.40, wspace=0.22,
)
ax_top = [fig.add_subplot(gs[0, i]) for i in range(3)]  # general spectrum
ax_bot = [fig.add_subplot(gs[1, i]) for i in range(3)]  # direct spectrum

# --- sliders (2 rows x 3 cols) --------------------------------------------
slider_y = [0.22, 0.14]
slider_x = [0.10, 0.41, 0.72]
slider_w = 0.24
slider_h = 0.025


def _add_slider(name, vmin, vmax, vinit, col, row):
    ax_s = fig.add_axes(
        (slider_x[col], slider_y[row], slider_w, slider_h),
        facecolor="lightgoldenrodyellow",
    )
    return Slider(ax_s, name, vmin, vmax, valinit=vinit)


sliders = {
    "c0":    _add_slider("c0",    1.0,   50.0,    10.0,    0, 0),
    "rho0":  _add_slider("rho0",  100.0, 3000.0,  1000.0,  1, 0),
    "delta": _add_slider("delta", 0.0,   1.0,     0.1,     2, 0),
    "H":     _add_slider("H",     0.001, 0.05,    0.012,   0, 1),
    "h":     _add_slider("h",     0.001, 0.05,    0.006,   1, 1),
    "theta": _add_slider("theta", 0.0,   np.pi,   0.0,     2, 1),
}

# --- radio buttons (view mode + kernel) + save button --------------------
rax = fig.add_axes((0.02, 0.02, 0.12, 0.11), facecolor="lightgoldenrodyellow")
rax.set_title("view", fontsize=8)
radio = RadioButtons(rax, ("General 1x3", "Direct 1x3", "All 2x3"), active=2)

krax = fig.add_axes((0.16, 0.02, 0.15, 0.11), facecolor="lightgoldenrodyellow")
krax.set_title("kernel", fontsize=8)
kernel_radio = RadioButtons(krax, tuple(KERNELS.keys()), active=0)

bax = fig.add_axes((0.88, 0.03, 0.10, 0.045))
save_btn = Button(bax, "Save")


# --- single update path ----------------------------------------------------
def update(_val=None) -> None:
    # Pull current slider values. H and h are independent -- changing H never
    # touches h.
    params.c0 = float(sliders["c0"].val)
    params.rho0 = float(sliders["rho0"].val)
    params.delta = float(sliders["delta"].val)
    params.H = float(sliders["H"].val)
    params.h = float(sliders["h"].val)
    params.theta = float(sliders["theta"].val)
    params.kernel = kernel_radio.value_selected

    mode = radio.value_selected

    gen = compute_all_general(params)
    dr = compute_all_direct(params)

    for ax in ax_top + ax_bot:
        ax.clear()

    if mode == "General 1x3":
        for ax in ax_bot:
            ax.set_visible(False)
        for ax in ax_top:
            ax.set_visible(True)
        _, _, mu_c = gen["cartesian"]
        _, _, mu_h = gen["hex"]
        plot_mu_panel(ax_top[0], mu_c, cfg,
                      cfg.color_cartesian_plus, cfg.color_cartesian_minus,
                      "Cartesian (general)", cfg.scatter_size)
        plot_mu_panel(ax_top[1], mu_h, cfg,
                      cfg.color_hex_plus, cfg.color_hex_minus,
                      "Hexagonal (general)", cfg.scatter_size)
        plot_mu_panel_random(ax_top[2], gen["random"], cfg, "Random (general)")
        fig.suptitle(
            f"General spectrum  |  kernel={params.kernel}  c0={params.c0:.1f}  delta={params.delta:.2f}  "
            f"H={params.H:.3f}  h={params.h:.3f}  theta={params.theta:.2f}",
            fontsize=cfg.title_fontsize,
        )

    elif mode == "Direct 1x3":
        for ax in ax_bot:
            ax.set_visible(False)
        for ax in ax_top:
            ax.set_visible(True)
        plot_eigvals_panel(ax_top[0], dr["cartesian"], cfg,
                           cfg.color_cartesian_direct, "Cartesian (direct)",
                           cfg.scatter_size)
        plot_eigvals_panel(ax_top[1], dr["hex"], cfg,
                           cfg.color_hex_direct, "Hexagonal (direct)",
                           cfg.scatter_size)
        plot_eigvals_panel_random(ax_top[2], dr["random"], cfg, "Random (direct)")
        fig.suptitle(
            f"Direct spectrum  |  kernel={params.kernel}  c0={params.c0:.1f}  rho0={params.rho0:.0f}  "
            f"delta={params.delta:.2f}  H={params.H:.3f}  h={params.h:.3f}",
            fontsize=cfg.title_fontsize,
        )

    else:  # "All 2x3"
        for ax in ax_top + ax_bot:
            ax.set_visible(True)
        _, _, mu_c = gen["cartesian"]
        _, _, mu_h = gen["hex"]
        plot_mu_panel(ax_top[0], mu_c, cfg,
                      cfg.color_cartesian_plus, cfg.color_cartesian_minus,
                      "Cartesian (general)", cfg.scatter_size)
        plot_mu_panel(ax_top[1], mu_h, cfg,
                      cfg.color_hex_plus, cfg.color_hex_minus,
                      "Hexagonal (general)", cfg.scatter_size)
        plot_mu_panel_random(ax_top[2], gen["random"], cfg, "Random (general)")
        plot_eigvals_panel(ax_bot[0], dr["cartesian"], cfg,
                           cfg.color_cartesian_direct, "Cartesian (direct)",
                           cfg.scatter_size)
        plot_eigvals_panel(ax_bot[1], dr["hex"], cfg,
                           cfg.color_hex_direct, "Hexagonal (direct)",
                           cfg.scatter_size)
        plot_eigvals_panel_random(ax_bot[2], dr["random"], cfg, "Random (direct)")
        fig.suptitle(
            f"General (top) + Direct (bottom)  |  kernel={params.kernel}  c0={params.c0:.1f}  "
            f"rho0={params.rho0:.0f}  delta={params.delta:.2f}  "
            f"H={params.H:.3f}  h={params.h:.3f}  theta={params.theta:.2f}",
            fontsize=cfg.title_fontsize,
        )

    fig.canvas.draw_idle()


def on_save(_event) -> None:
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    path = f"spectrum_gui_{ts}.png"
    fig.savefig(path, dpi=cfg.dpi)
    print(f"Saved {path}")


for s in sliders.values():
    s.on_changed(update)
radio.on_clicked(update)
kernel_radio.on_clicked(update)
save_btn.on_clicked(on_save)

# initial render
update()
plt.show()

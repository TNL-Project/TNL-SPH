"""WCSPH stability-spectrum analysis library (refactored for importability)."""
# noqa: SIZE_OK - monolithic SPH analysis script preserved as a single
# importable module by explicit user constraint: splitting into more than the
# 2 deliverables is forbidden and ALL math must be copied verbatim.
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kernels import (W_wendlandC2, dWdr_wendlandC2,
                     W_wendlandC4, dWdr_wendlandC4,
                     W_cubic, dWdr_cubic,
                     W_quartic, dWdr_quartic,
                     W_quintic, dWdr_quintic,
                     W_gaussian, dWdr_gaussian)


# =====================================================================
# 1. PHYSICAL / NUMERICAL PARAMETERS
# =====================================================================
@dataclass
class Parameters:
    c0: float = 10.0
    rho0: float = 1000.0
    delta: float = 0.1
    H: float = 1.2e-2
    h: float = 0.0
    theta: float = 0.0
    kernel: str = "wendlandC2"
    n_real: int = 10
    n_cell: int = 12
    n_images: int = 2
    xi_vals: np.ndarray | None = None
    nu: float = 0.0
    dim: int = 2
    eps_visc: float = 0.01

    def __post_init__(self) -> None:
        # h is INDEPENDENTLY adjustable (e.g. the GUI slider). The H/2 default
        # only applies when the user leaves h at its 0.0 sentinel; any explicit
        # value (including a GUI-driven one) is respected as-is.
        if self.h == 0.0:
            self.h = self.H / 2
        if self.xi_vals is None:
            self.xi_vals = np.linspace(0, np.pi, 300)


# =====================================================================
# 2. KERNEL REGISTRY
#    Each entry maps a name to its vector   ized W, dWdr, and support
#    radius factor (in units of H).  The    default ("wendlandC2") matches
#    the original monolithic script exactly.
# =====================================================================
W    = W_wendlandC2
dWdr = dWdr_wendlandC2
W_vec    = np.vectorize(W)
dWdr_vec = np.vectorize(dWdr)

KERNELS = {
    "wendlandC2": {"W_vec": W_vec,                      "dWdr_vec": dWdr_vec,                      "support": 2.0},
    "wendlandC4": {"W_vec": np.vectorize(W_wendlandC4), "dWdr_vec": np.vectorize(dWdr_wendlandC4), "support": 2.0},
    "cubic":      {"W_vec": np.vectorize(W_cubic),      "dWdr_vec": np.vectorize(dWdr_cubic),      "support": 2.0},
    "quartic":    {"W_vec": np.vectorize(W_quartic),    "dWdr_vec": np.vectorize(dWdr_quartic),    "support": 2.5},
    "quintic":    {"W_vec": np.vectorize(W_quintic),    "dWdr_vec": np.vectorize(dWdr_quintic),    "support": 3.0},
    "gaussian":   {"W_vec": np.vectorize(W_gaussian),   "dWdr_vec": np.vectorize(dWdr_gaussian),   "support": 3.0},
}


# =====================================================================
# 2b. VISCOSITY HELPER
#     Convert a target kinematic viscosity nu [m^2/s] into the dimensionless
#     Monaghan-Gingold coefficient alpha, using the standard equivalence
#     between the always-on (switch-dropped) alpha-viscosity and an
#     effective Laplacian/physical viscosity:
#
#         nu_eff = alpha*c0*H / (2*(dim+2))
#         =>  alpha = 2*(dim+2)*nu / (c0*H)
#
#     dim : int
#         Spatial dimension (2 for the 2D lattices used here, 3 for 3D SPH).
# =====================================================================
def alpha_from_nu(nu, c0, H, dim=2):
    return 2*(dim + 2)*nu / (c0*H)


# =====================================================================
# 3. NEIGHBOUR-LIST GENERATORS FOR DIFFERENT PARTICLE ARRANGEMENTS
#    all calibrated to the SAME mean number density n0 = 1/h^2
#    (i.e. Vj = h^2 for every particle in every arrangement) so the
#    comparison isolates the effect of geometry, not of density.
# =====================================================================
def neighbours_cartesian(h, H, margin=1, support_radius=None):
    R = support_radius if support_radius is not None else 2*H
    nx = int(np.ceil(R/h)) + margin
    ix, iy = np.meshgrid(np.arange(-nx, nx+1), np.arange(-nx, nx+1))
    rx = (ix*h).ravel(); ry = (iy*h).ravel()
    r  = np.sqrt(rx**2 + ry**2)
    mask = (r > 1e-12) & (r < R)
    Vj = np.full(mask.sum(), h**2)
    return rx[mask], ry[mask], r[mask], Vj

def neighbours_hex(h, H, margin=2, support_radius=None):
    # triangular lattice; nearest-neighbour distance a chosen so that the
    # area per particle (sqrt(3)/2 * a^2) equals h^2 -> same density as cartesian
    R = support_radius if support_radius is not None else 2*H
    a  = h * np.sqrt(2/np.sqrt(3))
    a1 = np.array([a, 0.0])
    a2 = np.array([a/2, a*np.sqrt(3)/2])
    n  = int(np.ceil(R/a)) + margin
    ii, jj = np.meshgrid(np.arange(-n, n+1), np.arange(-n, n+1))
    pts = ii.ravel()[:, None]*a1 + jj.ravel()[:, None]*a2
    r  = np.linalg.norm(pts, axis=1)
    mask = (r > 1e-12) & (r < R)
    Vj = np.full(mask.sum(), h**2)
    return pts[mask, 0], pts[mask, 1], r[mask], Vj

def neighbours_random(h, H, jitter_frac=0.3, n_cell=16, tile_range=1, seed=0, support_radius=None):
    # start from a cartesian grid (same density), jitter every particle,
    # then tile that SAME disordered patch periodically (so the medium
    # really is periodic, just with a disordered unit cell)
    R = support_radius if support_radius is not None else 2*H
    rng = np.random.default_rng(seed)
    idx = np.arange(-(n_cell//2), n_cell - n_cell//2)
    ix, iy = np.meshgrid(idx, idx)
    base_x = ix.ravel()*h; base_y = iy.ravel()*h
    px = base_x + rng.uniform(-jitter_frac*h, jitter_frac*h, size=base_x.shape)
    py = base_y + rng.uniform(-jitter_frac*h, jitter_frac*h, size=base_y.shape)
    L = n_cell*h
    ref = np.argmin(base_x**2 + base_y**2)   # reference particle nearest cell centre
    refx, refy = px[ref], py[ref]
    rx_list, ry_list = [], []
    for ti in range(-tile_range, tile_range+1):
        for tj in range(-tile_range, tile_range+1):
            rx_list.append(px + ti*L - refx)
            ry_list.append(py + tj*L - refy)
    rx = np.concatenate(rx_list); ry = np.concatenate(ry_list)
    r  = np.sqrt(rx**2 + ry**2)
    mask = (r > 1e-12) & (r < R)
    Vj = np.full(mask.sum(), h**2)
    return rx[mask], ry[mask], r[mask], Vj


# =====================================================================
# 4. DIRECT (NON-FOURIER) PATCH BUILDERS
#    Build a real periodic particle patch and (later) diagonalise the full
#    3N x 3N operator directly -- no plane-wave ansatz, no xi_vals loop.
# =====================================================================
def build_patch_cartesian(n_cell, h):
    idx = np.arange(n_cell)
    ix, iy = np.meshgrid(idx, idx)
    x = (ix.ravel()*h).astype(float)
    y = (iy.ravel()*h).astype(float)
    L = n_cell*h
    Vj = np.full(len(x), h**2)
    return x, y, L, Vj

def build_patch_hex(n_cell, h):
    a  = h * np.sqrt(2/np.sqrt(3))
    a1 = np.array([a, 0.0])
    a2 = np.array([a/2, a*np.sqrt(3)/2])
    idx = np.arange(n_cell)
    ii, jj = np.meshgrid(idx, idx)
    pts = ii.ravel()[:, None]*a1 + jj.ravel()[:, None]*a2
    x, y = pts[:, 0], pts[:, 1]
    L = n_cell*a   # patch is only approximately periodic in (x,y) for hex;
                   # fine for n_images>=2 neighbour search, just not exactly
                   # tileable like the cartesian case
    Vj = np.full(len(x), h**2)
    return x, y, L, Vj

def build_patch_random(n_cell, h, jitter_frac=0.3, seed=0):
    rng = np.random.default_rng(seed)
    idx = np.arange(n_cell)
    ix, iy = np.meshgrid(idx, idx)
    base_x = ix.ravel()*h; base_y = iy.ravel()*h
    x = base_x + rng.uniform(-jitter_frac*h, jitter_frac*h, size=base_x.shape)
    y = base_y + rng.uniform(-jitter_frac*h, jitter_frac*h, size=base_y.shape)
    L = n_cell*h
    Vj = np.full(len(x), h**2)
    return x, y, L, Vj


# =====================================================================
# 5. FOURIER SYMBOLS G_hat(xi), D_hat(xi), AND EIGENVALUES mu_pm(xi)
# =====================================================================
def compute_spectrum(rx, ry, r, Vj, H, h, xi_vals, c0, delta,
                      nu=0.0, dim=2, eps_visc=0.01,
                      W_fn=None, dWdr_fn=None):
    """
    Von Neumann (Fourier symbol) stability analysis of delta-WCSPH,
    including the standard Monaghan-Gingold artificial viscosity.

    nu : float
        Kinematic viscosity [m^2/s]. Converted to the dimensionless
        Monaghan-Gingold coefficient via
            alpha = alpha_from_nu(nu, c0, H, dim) = 2*(dim+2)*nu/(c0*H)
        nu=0 recovers the original inviscid spectrum exactly.
    dim : int
        Spatial dimension used in the nu -> alpha conversion (default 2).
    eps_visc : float
        Regularisation constant in mu_ij's denominator (standard 0.01).

    Linearisation: the true Pi_ij is switched off unless v_ij.r_ij < 0
    (approaching particles) and carries a quadratic beta*mu_ij**2 term.
    Both are nonlinear in the velocity perturbation about the quiescent
    background (u=0), so they cannot enter a constant-coefficient linear
    operator. We keep only the always-on, linear-in-mu_ij piece -- the
    known equivalence between "always-on" alpha-viscosity and an effective
    Laplacian/physical viscosity (Monaghan 1992, 2005).

    Returns G_hat, D_hat, PI_hat, mu_pm (4 values; PI_hat added).
    """
    alpha = alpha_from_nu(nu, c0, H, dim)

    _W  = W_fn if W_fn is not None else W_vec
    _dW = dWdr_fn if dWdr_fn is not None else dWdr_vec
    Wr  = _W(r, H)
    dWr = _dW(r, H)
    dWx = dWr * rx / r

    # Linearised viscous acceleration on particle i (x-momentum, background u=0):
    #   a_i,x^visc = alpha*c0*H**2 * sum_j Vj * dWr*rx**2/(r*(r**2+eps_visc*H**2)) * (ux_i - ux_j)
    # Difference (Laplacian-type) quantity like D_hat -- needs (1-phase)
    # treatment so it vanishes for a uniform velocity field.
    visc_w = dWr * rx**2 / (r * (r**2 + eps_visc*H**2))

    G_hat  = np.zeros(len(xi_vals), dtype=complex)
    D_hat  = np.zeros(len(xi_vals), dtype=complex)
    PI_hat = np.zeros(len(xi_vals), dtype=complex)
    for idx, xi in enumerate(xi_vals):
        k = xi / h
        phase = np.exp(1j*k*rx)
        G_hat[idx]  = np.sum(Vj*dWx*phase)
        D_hat[idx]  = np.sum(Vj*(Wr/r**2)*(1 - phase))
        PI_hat[idx] = np.sum(Vj*visc_w*(1 - phase))
    PI_hat *= alpha*c0*H**2

    mu0_rho  = -delta*H*c0*D_hat
    mu0_visc = PI_hat
    #mu0_visc = -nu*D_hat

    trace = mu0_rho + mu0_visc
    disc  = ((mu0_rho - mu0_visc)/2)**2 - c0**2*np.abs(G_hat)**2 + 0j
    sq    = np.sqrt(disc)
    mu_pm = trace/2 + np.array([1, -1])[:, None]*sq
    return G_hat, D_hat, PI_hat, mu_pm


#FIXME: Again, the general vector doesnt work
def compute_spectrum_general(rx, ry, r, Vj, H, h, xi_vals, theta, c0, delta, W_fn=None, dWdr_fn=None):
    """
    Same as compute_spectrum, but for a wave vector k = (k cos(theta), k sin(theta))
    instead of being restricted to k = (k, 0).

    theta : propagation angle of the wave, in radians (theta=0 reproduces
            the original x-only compute_spectrum exactly).

    NOTE: Gx_hat and Gy_hat are computed as independent Cartesian components
    of the kernel-gradient sum (each with the FULL 2D phase factor), then
    combined as Gx_hat**2 + Gy_hat**2 in the discriminant. This is NOT the
    same as projecting the kernel gradient onto k_hat first and squaring
    that scalar -- the two only coincide for a perfectly isotropic kernel
    sum, which is generally false on a finite lattice. Using the projected
    scalar would artificially erase the very anisotropy this function is
    meant to reveal.
    """
    _W  = W_fn if W_fn is not None else W_vec
    _dW = dWdr_fn if dWdr_fn is not None else dWdr_vec
    Wr  = _W(r, H)
    dWr = _dW(r, H)
    fx = dWr * rx / r        # x-component of kernel gradient
    fy = dWr * ry / r        # y-component of kernel gradient

    kdir_x, kdir_y = np.cos(theta), np.sin(theta)

    Gx_hat = np.zeros(len(xi_vals), dtype=complex)
    Gy_hat = np.zeros(len(xi_vals), dtype=complex)
    D_hat  = np.zeros(len(xi_vals), dtype=complex)
    for idx, xi in enumerate(xi_vals):
        k = xi / h
        kx, ky = k*kdir_x, k*kdir_y
        phase = np.exp(1j*(kx*rx + ky*ry))
        Gx_hat[idx] = np.sum(Vj*fx*phase)
        Gy_hat[idx] = np.sum(Vj*fy*phase)
        D_hat[idx]  = np.sum(Vj*(Wr/r**2)*(1 - phase))

    # mu0  = -delta*H*c0*D_hat
    # disc = (mu0/2)**2 + c0**2*(Gx_hat**2 + Gy_hat**2) + 0j
    # #disc = (mu0/2)**2 + c0**2 * (np.abs(Gx_hat)**2 + np.abs(Gy_hat)**2) + 0j
    # sq   = np.sqrt(disc)
    # mu_pm = mu0/2 + np.array([1, -1])[:, None]*sq

    # # G_hat returned here only for diagnostics/back-compat with compute_spectrum's
    # # 3-value signature -- it's the rotation-invariant combined magnitude,
    # # not a physically meaningful "scalar projection".
    # G_hat = np.sqrt(Gx_hat**2 + Gy_hat**2)
    # return G_hat, D_hat, mu_pm

    G_hat = kdir_x*Gx_hat + kdir_y*Gy_hat      # projection onto propagation direction
    mu0   = -delta*H*c0*D_hat
    disc  = (mu0/2)**2 - c0**2*np.abs(G_hat)**2 + 0j
    sq    = np.sqrt(disc)
    mu_pm = mu0/2 + np.array([1, -1])[:, None]*sq
    return G_hat, D_hat, mu_pm

def get_spectrum_direct(x, y, L, Vj, H, h, c0, delta, rho0=1000.0, n_images=2,
                        nu=0.0, dim=2, eps_visc=0.01,
                        W_fn=None, dWdr_fn=None, support_radius=None):
    """
    Directly assemble the linear stability operator from particle-particle
    interactions (periodic box L x L, n_images periodic shells for neighbour
    search) and return eigenvalues via np.linalg.eigvals.

    nu : float
        Kinematic viscosity [m^2/s]. Converted to alpha via alpha_from_nu.
        nu=0 recovers the original inviscid operator exactly.
    dim : int
        Spatial dimension used in the nu -> alpha conversion (default 2).
    eps_visc : float
        Regularisation constant in mu_ij's denominator (standard 0.01).

    Viscosity here is NOT restricted to the longitudinal (1D) case: since
    v_ij . r_ij mixes both velocity components, the linearised Monaghan
    viscosity couples ux and uy through the full off-diagonal blocks
    L_ux_uy / L_uy_ux, not just the diagonal L_ux_ux / L_uy_uy blocks.
    """
    alpha = alpha_from_nu(nu, c0, H, dim)

    _W  = W_fn if W_fn is not None else W_vec
    _dW = dWdr_fn if dWdr_fn is not None else dWdr_vec
    _R  = support_radius if support_radius is not None else 2*H
    N = len(x)
    Gx = np.zeros((N, N))
    Gy = np.zeros((N, N))
    D  = np.zeros((N, N))
    Pxx = np.zeros((N, N))
    Pxy = np.zeros((N, N))
    Pyy = np.zeros((N, N))

    for i in range(N):
        dxs, dys = [], []
        for sx in range(-n_images, n_images+1):
            for sy in range(-n_images, n_images+1):
                dxs.append(x - x[i] + sx*L)
                dys.append(y - y[i] + sy*L)
        dx = np.concatenate(dxs)
        dy = np.concatenate(dys)
        jidx = np.tile(np.arange(N), (2*n_images+1)**2)
        Vjr  = np.tile(Vj, (2*n_images+1)**2)

        r = np.sqrt(dx**2 + dy**2)
        mask = (r > 1e-12) & (r < _R)
        dx, dy, r, jidx, Vjr = dx[mask], dy[mask], r[mask], jidx[mask], Vjr[mask]

        Wr  = _W(r, H)
        dWr = _dW(r, H)
        rx = -dx; ry = -dy
        wx = dWr * rx / r
        wy = dWr * ry / r
        wD = Wr / r**2

        np.add.at(Gx[i], jidx, Vjr*wx)
        np.add.at(Gy[i], jidx, Vjr*wy)
        np.add.at(D[i],  jidx, Vjr*wD)

        # Monaghan-Gingold viscosity, always-on linearisation:
        # w_ij = alpha*c0*H**2*Vj*dWr/(r*(r**2+eps_visc*H**2))
        # coupling matrices Pxx=w*rx**2, Pxy=w*rx*ry, Pyy=w*ry**2
        visc_k = alpha*c0*H**2 * dWr / (r * (r**2 + eps_visc*H**2))
        np.add.at(Pxx[i], jidx, Vjr*visc_k*rx**2)
        np.add.at(Pxy[i], jidx, Vjr*visc_k*rx*ry)
        np.add.at(Pyy[i], jidx, Vjr*visc_k*ry**2)

    diagD = np.diag(D.sum(axis=1))
    L_rho_rho = delta*H*c0*(D - diagD)
    L_rho_ux  = -rho0*Gx
    L_rho_uy  = -rho0*Gy
    L_ux_rho  = -(c0**2/rho0)*Gx
    L_uy_rho  = -(c0**2/rho0)*Gy

    # Viscous velocity blocks: acceleration is a_i = w_ij*(u_i - u_j),
    # so row i picks up +rowsum on the diagonal and -w_ij off-diagonal --
    # the OPPOSITE sign convention to (D - diagD). Hence diag(rowsum) - P.
    L_ux_ux = np.diag(Pxx.sum(axis=1)) - Pxx
    L_ux_uy = np.diag(Pxy.sum(axis=1)) - Pxy
    L_uy_ux = np.diag(Pxy.sum(axis=1)) - Pxy   # Pxy symmetric (rx*ry) by construction
    L_uy_uy = np.diag(Pyy.sum(axis=1)) - Pyy

    Ltot = np.block([
        [L_rho_rho, L_rho_ux, L_rho_uy],
        [L_ux_rho,  L_ux_ux,  L_ux_uy ],
        [L_uy_rho,  L_uy_ux,  L_uy_uy ],
    ])
    eigvals = np.linalg.eigvals(Ltot)
    return eigvals, Ltot


# =====================================================================
# 6. SHARED PLOT CONFIGURATION
# =====================================================================
@dataclass
class PlotConfig:
    figsize_1x3: tuple = (15, 5.2)
    figsize_2x3: tuple = (15, 10.4)
    figsize_single: tuple = (7.5, 6.5)
    dpi: int = 150
    scatter_size: int = 10
    scatter_size_random: int = 5
    alpha_random: float = 0.4
    color_cartesian_plus: str = 'darkorange'
    color_cartesian_minus: str = 'orangered'
    color_hex_plus: str = 'teal'
    color_hex_minus: str = 'darkslategray'
    color_random: str = 'lightcoral'
    color_cartesian_direct: str = 'darkorange'
    color_hex_direct: str = 'teal'
    xlabel: str = 'Re(mu)  [1/s]'
    ylabel: str = 'Im(mu)  [1/s]'
    zero_line_color: str = 'gray'
    zero_line_lw: float = 0.6
    legend_fontsize: int = 8
    title_fontsize: int = 11


# =====================================================================
# 7. CONVENIENCE COMPUTE-ALL HELPERS
# =====================================================================
def compute_all(params: Parameters) -> dict:
    """Fourier spectra (with viscosity) for cartesian / hex / random."""
    kinfo = KERNELS[params.kernel]
    kW, kdw, kR = kinfo["W_vec"], kinfo["dWdr_vec"], kinfo["support"] * params.H
    out: dict = {}
    rxC, ryC, rC, VjC = neighbours_cartesian(params.h, params.H, support_radius=kR)
    G, D, _PI, mu = compute_spectrum(
        rxC, ryC, rC, VjC, params.H, params.h, params.xi_vals,
        params.c0, params.delta, params.nu, params.dim, params.eps_visc, kW, kdw)
    out['cartesian'] = (G, D, mu)
    rxH, ryH, rH, VjH = neighbours_hex(params.h, params.H, support_radius=kR)
    G, D, _PI, mu = compute_spectrum(
        rxH, ryH, rH, VjH, params.H, params.h, params.xi_vals,
        params.c0, params.delta, params.nu, params.dim, params.eps_visc, kW, kdw)
    out['hex'] = (G, D, mu)
    rnd = []
    for s in range(params.n_real):
        rxR, ryR, rR, VjR = neighbours_random(params.h, params.H, seed=s, support_radius=kR)
        G, D, _PI, mu = compute_spectrum(
            rxR, ryR, rR, VjR, params.H, params.h, params.xi_vals,
            params.c0, params.delta, params.nu, params.dim, params.eps_visc, kW, kdw)
        rnd.append((G, D, mu))
    out['random'] = rnd
    return out


def compute_all_general(params: Parameters) -> dict:
    """Generalised Fourier spectra for cartesian / hex / random at params.theta."""
    kinfo = KERNELS[params.kernel]
    kW, kdw, kR = kinfo["W_vec"], kinfo["dWdr_vec"], kinfo["support"] * params.H
    out: dict = {}
    rxC, ryC, rC, VjC = neighbours_cartesian(params.h, params.H, support_radius=kR)
    out['cartesian'] = compute_spectrum_general(
        rxC, ryC, rC, VjC, params.H, params.h, params.xi_vals,
        params.theta, params.c0, params.delta, kW, kdw)
    rxH, ryH, rH, VjH = neighbours_hex(params.h, params.H, support_radius=kR)
    out['hex'] = compute_spectrum_general(
        rxH, ryH, rH, VjH, params.H, params.h, params.xi_vals,
        params.theta, params.c0, params.delta, kW, kdw)
    rnd = []
    for s in range(params.n_real):
        rxR, ryR, rR, VjR = neighbours_random(params.h, params.H, seed=s, support_radius=kR)
        rnd.append(compute_spectrum_general(
            rxR, ryR, rR, VjR, params.H, params.h, params.xi_vals,
            params.theta, params.c0, params.delta, kW, kdw))
    out['random'] = rnd
    return out


def compute_all_direct(params: Parameters) -> dict:
    """Direct (non-Fourier) eigenvalue spectra for cartesian / hex / random."""
    kinfo = KERNELS[params.kernel]
    kW, kdw, kR = kinfo["W_vec"], kinfo["dWdr_vec"], kinfo["support"] * params.H
    out: dict = {}
    xC, yC, LC, VjC = build_patch_cartesian(params.n_cell, params.h)
    out['cartesian'] = get_spectrum_direct(
        xC, yC, LC, VjC, params.H, params.h, params.c0, params.delta,
        params.rho0, params.n_images, params.nu, params.dim, params.eps_visc,
        kW, kdw, kR)[0]
    xH, yH, LH, VjH = build_patch_hex(params.n_cell, params.h)
    out['hex'] = get_spectrum_direct(
        xH, yH, LH, VjH, params.H, params.h, params.c0, params.delta,
        params.rho0, params.n_images, params.nu, params.dim, params.eps_visc,
        kW, kdw, kR)[0]
    rnd = []
    for s in range(params.n_real):
        xR, yR, LR, VjR = build_patch_random(params.n_cell, params.h, seed=s)
        rnd.append(get_spectrum_direct(
            xR, yR, LR, VjR, params.H, params.h, params.c0, params.delta,
            params.rho0, params.n_images, params.nu, params.dim, params.eps_visc,
            kW, kdw, kR)[0])
    out['random'] = np.concatenate(rnd)
    return out


# =====================================================================
# 8. PLOT HELPERS (operate on a given Axes, styled via PlotConfig)
# =====================================================================
def plot_mu_panel(ax, mu_pm, config: PlotConfig, color_plus, color_minus,
                  title, scatter_size=None):
    s = config.scatter_size if scatter_size is None else scatter_size
    ax.scatter(mu_pm[0].real, mu_pm[0].imag, s=s, c=color_plus)
    ax.scatter(mu_pm[1].real, mu_pm[1].imag, s=s, c=color_minus)
    ax.axhline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.axvline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.set_xlabel(config.xlabel)
    ax.set_ylabel(config.ylabel)
    ax.set_title(title)


def plot_mu_panel_random(ax, mu_list, config: PlotConfig, title):
    # mu_list: list of (G_hat, D_hat, mu_pm) tuples (one per random realisation)
    for _G, _D, mu in mu_list:
        ax.scatter(mu[0].real, mu[0].imag, s=config.scatter_size_random,
                   c=config.color_random, alpha=config.alpha_random)
        ax.scatter(mu[1].real, mu[1].imag, s=config.scatter_size_random,
                   c=config.color_random, alpha=config.alpha_random)
    ax.axhline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.axvline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.set_xlabel(config.xlabel)
    ax.set_ylabel(config.ylabel)
    ax.set_title(title)


def plot_eigvals_panel(ax, eigvals, config: PlotConfig, color, title,
                       scatter_size=None):
    s = config.scatter_size if scatter_size is None else scatter_size
    ax.scatter(eigvals.real, eigvals.imag, s=s, c=color)
    ax.axhline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.axvline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.set_xlabel(config.xlabel)
    ax.set_ylabel(config.ylabel)
    ax.set_title(title)


def plot_eigvals_panel_random(ax, eigvals, config: PlotConfig, title):
    # eigvals: already-concatenated array across all random seeds
    ax.scatter(eigvals.real, eigvals.imag, s=config.scatter_size_random,
               c=config.color_random, alpha=config.alpha_random)
    ax.axhline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.axvline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.set_xlabel(config.xlabel)
    ax.set_ylabel(config.ylabel)
    ax.set_title(title)


# =====================================================================
# 9. STATIC PLOT GENERATION (original script behaviour, import-safe)
# =====================================================================
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    p = Parameters()
    cfg = PlotConfig()
    kinfo = KERNELS[p.kernel]
    kW, kdw, kR = kinfo["W_vec"], kinfo["dWdr_vec"], kinfo["support"] * p.H

    # --- General (Fourier) spectra; original wiring preserved exactly:
    #     cartesian uses compute_spectrum_general with theta=pi/2,
    #     hex/random use plain compute_spectrum. ---
    rxC, ryC, rC, VjC = neighbours_cartesian(p.h, p.H, support_radius=kR)
    G_c, D_c, _PI_c, mu_c = compute_spectrum(rxC, ryC, rC, VjC, p.H, p.h, p.xi_vals, p.c0, p.delta, kW, kdw)
    G_c, D_c, mu_c = compute_spectrum_general(rxC, ryC, rC, VjC, p.H, p.h, p.xi_vals, np.pi/2, p.c0, p.delta, kW, kdw)
    #plt.scatter( rxH, ryH )
    #plt.show()
    print(f"Cartesian : N_neighbours={len(rC):4d}  max Re(mu)={mu_c.real.max():.3e}  max|Im(mu)|={np.abs(mu_c.imag).max():.3e}")

    rxH, ryH, rH, VjH = neighbours_hex(p.h, p.H, support_radius=kR)
    G_h, D_h, _PI_h, mu_h = compute_spectrum(rxH, ryH, rH, VjH, p.H, p.h, p.xi_vals, p.c0, p.delta, kW, kdw)
    #G_h, D_h, mu_h = compute_spectrum_general(rxH, ryH, rH, VjH, p.H, p.h, p.xi_vals, np.pi/4, c0, delta)
    #plt.scatter( rxH, ryH )
    #plt.show()
    print(f"Hexagonal : N_neighbours={len(rH):4d}  max Re(mu)={mu_h.real.max():.3e}  max|Im(mu)|={np.abs(mu_h.imag).max():.3e}")

    mu_rand_list = []
    mu_rand_arr = []
    for s in range(p.n_real):
        rxR, ryR, rR, VjR = neighbours_random(p.h, p.H, seed=s, support_radius=kR)
        #plt.scatter( rxR, ryR )
        #plt.show()
        G_r, D_r, _PI_r, mu_r = compute_spectrum(rxR, ryR, rR, VjR, p.H, p.h, p.xi_vals, p.c0, p.delta, kW, kdw)
        #_, _, mu_r = compute_spectrum_general(rxR, ryR, rR, VjR, p.H, p.h, p.xi_vals, np.pi/4, c0, delta)
        mu_rand_list.append((G_r, D_r, mu_r))
        mu_rand_arr.append(mu_r)
    mu_rand_all = np.array(mu_rand_arr)   # shape (n_real, 2, Nxi)
    print(f"Random    : N_neighbours~{len(rR):4d} (last seed)  "
          f"max Re(mu)={mu_rand_all.real.max():.3e}  max|Im(mu)|={np.abs(mu_rand_all.imag).max():.3e}")

    # --- fig 1: raw spectrum comparison across particle arrangements ---
    fig, ax = plt.subplots(figsize=cfg.figsize_single, dpi=cfg.dpi)

    for s in range(p.n_real):
        mu_r = mu_rand_all[s]
        ax.scatter(mu_r[0].real, mu_r[0].imag, s=cfg.scatter_size_random, c=cfg.color_random, alpha=cfg.alpha_random,
                   label=f'random ({p.n_real} seeds)' if s == 0 else None)
        ax.scatter(mu_r[1].real, mu_r[1].imag, s=cfg.scatter_size_random, c=cfg.color_random, alpha=cfg.alpha_random)

    ax.scatter(mu_c[0].real, mu_c[0].imag, s=cfg.scatter_size, c=cfg.color_cartesian_plus,  label='cartesian  mu_+')
    ax.scatter(mu_c[1].real, mu_c[1].imag, s=cfg.scatter_size, c=cfg.color_cartesian_minus,   label='cartesian  mu_-')
    ax.scatter(mu_h[0].real, mu_h[0].imag, s=cfg.scatter_size, c=cfg.color_hex_plus,        label='hexagonal  mu_+')
    ax.scatter(mu_h[1].real, mu_h[1].imag, s=cfg.scatter_size, c=cfg.color_hex_minus, label='hexagonal  mu_-')

    ax.axhline(0, color=cfg.zero_line_color, lw=cfg.zero_line_lw); ax.axvline(0, color=cfg.zero_line_color, lw=cfg.zero_line_lw)
    ax.set_xlabel(cfg.xlabel); ax.set_ylabel(cfg.ylabel)
    ax.set_title('Raw spectrum of L(k) for different particle arrangements\n'
                 '(same number density n0 = 1/h^2 in all cases, k along x)')

    mu_rand_mean = mu_rand_all.mean(axis=0)   # not plotted; kept for reference only

    print(f"\nWorst-case damping  max|Re(mu)|:  cartesian={-mu_c.real.min():.1f}  "
          f"hexagonal={-mu_h.real.min():.1f}  random={-mu_rand_all.real.min():.1f}  [1/s]")
    print(f"Worst-case oscillation max|Im(mu)|: cartesian={np.abs(mu_c.imag).max():.1f}  "
          f"hexagonal={np.abs(mu_h.imag).max():.1f}  random={np.abs(mu_rand_all.imag).max():.1f}  [1/s]")

    ax.legend(fontsize=cfg.legend_fontsize, loc='upper left')
    plt.tight_layout()
    plt.savefig('raw_spectrum_comparison.png', dpi=cfg.dpi)
    plt.close(fig)
    print("Saved raw_spectrum_comparison.png")

    # --- fig 2: side-by-side panels (one geometry per axis) ---
    fig2, axs = plt.subplots(1, 3, figsize=cfg.figsize_1x3, dpi=cfg.dpi, sharex=True, sharey=True)

    plot_mu_panel(axs[0], mu_c, cfg, cfg.color_cartesian_plus, cfg.color_cartesian_minus,
                  'Cartesian (square) lattice', cfg.scatter_size)
    plot_mu_panel(axs[1], mu_h, cfg, cfg.color_hex_plus, cfg.color_hex_minus,
                  'Hexagonal (triangular) lattice', cfg.scatter_size)
    plot_mu_panel_random(axs[2], mu_rand_list, cfg,
                         f'Random, periodically repeating\n({p.n_real} disordered realisations)')

    fig2.suptitle('Effect of particle arrangement on the semi-discrete spectrum L(k)  (k along x)',
                  fontsize=cfg.title_fontsize)
    plt.tight_layout()
    plt.savefig('raw_spectrum_panels.png', dpi=cfg.dpi)
    plt.close(fig2)
    print("Saved raw_spectrum_panels.png")

    # --- Direct (non-Fourier) spectra ---
    xC, yC, LC, VjC_patch = build_patch_cartesian(p.n_cell, p.h)
    eig_c_direct, _ = get_spectrum_direct(xC, yC, LC, VjC_patch, p.H, p.h, p.c0, p.delta, p.rho0, p.n_images, kW, kdw, kR)
    print(f"Cartesian (direct): N={len(xC):4d}  max Re(mu)={eig_c_direct.real.max():.3e}  "
          f"max|Im(mu)|={np.abs(eig_c_direct.imag).max():.3e}")

    xH, yH, LH, VjH_patch = build_patch_hex(p.n_cell, p.h)
    eig_h_direct, _ = get_spectrum_direct(xH, yH, LH, VjH_patch, p.H, p.h, p.c0, p.delta, p.rho0, p.n_images, kW, kdw, kR)
    print(f"Hexagonal (direct): N={len(xH):4d}  max Re(mu)={eig_h_direct.real.max():.3e}  "
          f"max|Im(mu)|={np.abs(eig_h_direct.imag).max():.3e}")

    eig_rand_direct = []
    for s in range(p.n_real):
        xR, yR, LR, VjR_patch = build_patch_random(p.n_cell, p.h, seed=s)
        eig_r, _ = get_spectrum_direct(xR, yR, LR, VjR_patch, p.H, p.h, p.c0, p.delta, p.rho0, p.n_images, kW, kdw, kR)
        eig_rand_direct.append(eig_r)
    eig_rand_direct = np.concatenate(eig_rand_direct)
    print(f"Random (direct)   : max Re(mu)={eig_rand_direct.real.max():.3e}  "
          f"max|Im(mu)|={np.abs(eig_rand_direct.imag).max():.3e}")

    #def filter_nonzero(eig, tol=1e-6):
    #    return eig[np.abs(eig) > tol]
    #
    #eig_c_acoustic    = filter_nonzero(eig_c_direct)
    #eig_h_acoustic    = filter_nonzero(eig_h_direct)
    #eig_rand_acoustic = filter_nonzero(eig_rand_direct)

    # --- fig 3: direct-matrix spectrum comparison ---
    fig3, ax3 = plt.subplots(figsize=cfg.figsize_single, dpi=cfg.dpi)

    ax3.scatter(eig_rand_direct.real, eig_rand_direct.imag, s=cfg.scatter_size_random, c=cfg.color_random, alpha=cfg.alpha_random,
                label=f'random (direct, {p.n_real} seeds)')
    ax3.scatter(eig_c_direct.real, eig_c_direct.imag, s=cfg.scatter_size, c=cfg.color_cartesian_direct, label='cartesian (direct)')
    ax3.scatter(eig_h_direct.real, eig_h_direct.imag, s=cfg.scatter_size, c=cfg.color_hex_direct, label='hexagonal (direct)')

    ax3.axhline(0, color=cfg.zero_line_color, lw=cfg.zero_line_lw); ax3.axvline(0, color=cfg.zero_line_color, lw=cfg.zero_line_lw)
    ax3.set_xlabel(cfg.xlabel); ax3.set_ylabel(cfg.ylabel)
    ax3.set_title('Direct-matrix eigenvalue spectrum\n(no Fourier ansatz, full 3N x 3N operator, periodic patches)')
    ax3.legend(fontsize=cfg.legend_fontsize, loc='upper left')
    plt.tight_layout()
    plt.savefig('direct_spectrum_comparison.png', dpi=cfg.dpi)
    plt.close(fig3)
    print("Saved direct_spectrum_comparison.png")

    # --- fig 4: direct-matrix spectrum panels ---
    fig4, axs4 = plt.subplots(1, 3, figsize=cfg.figsize_1x3, dpi=cfg.dpi, sharex=True, sharey=True)

    plot_eigvals_panel(axs4[0], eig_c_direct, cfg, cfg.color_cartesian_direct,
                       'Cartesian (direct)', cfg.scatter_size)
    plot_eigvals_panel(axs4[1], eig_h_direct, cfg, cfg.color_hex_direct,
                       'Hexagonal (direct)', cfg.scatter_size)
    plot_eigvals_panel_random(axs4[2], eig_rand_direct, cfg,
                              f'Random (direct, {p.n_real} realisations)')

    fig4.suptitle('Direct (non-Fourier) eigenvalue spectrum of the assembled 3N x 3N operator',
                  fontsize=cfg.title_fontsize)
    plt.tight_layout()
    plt.savefig('direct_spectrum_panels.png', dpi=cfg.dpi)
    plt.close(fig4)
    print("Saved direct_spectrum_panels.png")

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
    """
    Convert a target kinematic viscosity nu [m^2/s] into the dimensionless
    Monaghan-Gingold coefficient alpha, using the standard equivalence
    between the always-on (switch-dropped) alpha-viscosity and an
    effective Laplacian/physical viscosity:

        nu_eff = alpha*c0*H / (2*(dim+2))
        =>  alpha = 2*(dim+2)*nu / (c0*H)

    dim : int
        Spatial dimension (2 for the 2D lattices used here, 3 for 3D SPH).
    """
    return 2*(dim + 2)*nu / (c0*H)


def solve_cubic_vectorized(p2, p1, p0):
    """
    Closed-form (Cardano) roots of x**3 + p2*x**2 + p1*x + p0 = 0, fully
    vectorised over arrays of (in general complex) coefficients p2, p1, p0.
    Returns an array of shape (*p2.shape, 3) with the three roots.

    Standard depressed-cubic substitution x = t - p2/3 gives
        t**3 + p*t + q = 0,   p = p1 - p2**2/3,   q = 2*p2**3/27 - p2*p1/3 + p0
    with roots t_k = omega**k*u + omega**(-k)*v  (k=0,1,2), omega=exp(2j*pi/3),
    where u is *any* cube root of (-q/2 + sqrt((q/2)**2+(p/3)**3)) and v is
    then fixed via v = -p/(3*u) -- NOT independently cube-rooted. This is
    the standard trick that avoids the classic "9 candidate pairs, only 3
    valid" branch-matching problem of naively cube-rooting both terms.

    A tiny-|u| fallback handles the (rare, measure-zero) degenerate case
    where u ~ 0; np.linalg.eigvals on the corresponding 3x3 matrix remains
    a safe way to double check any individual mode if ever in doubt.
    """
    p2 = np.asarray(p2, dtype=complex)
    p1 = np.asarray(p1, dtype=complex)
    p0 = np.asarray(p0, dtype=complex)

    p = p1 - p2**2/3
    q = 2*p2**3/27 - p2*p1/3 + p0
    disc = (q/2)**2 + (p/3)**3
    sqrt_disc = np.sqrt(disc)

    u3 = -q/2 + sqrt_disc
    u = u3**(1/3)   # principal complex cube root

    tiny = np.abs(u) < 1e-12
    v_direct = (-q/2 - sqrt_disc)**(1/3)
    u_safe = np.where(tiny, 1.0, u)
    v = np.where(tiny, v_direct, -p/(3*u_safe))

    omega = np.exp(2j*np.pi/3)
    t0 = u + v
    t1 = omega*u + np.conj(omega)*v
    t2 = np.conj(omega)*u + omega*v

    roots = np.stack([t0, t1, t2], axis=-1) - p2[..., None]/3
    return roots


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
    including a viscous term. Two standard, physically distinct viscosity
    models are provided -- choose ONE by commenting/uncommenting the
    corresponding block marked "OPTION A" / "OPTION B" below.

    OPTION A -- Monaghan-Gingold artificial viscosity (default, active)
        Pi_ij = -alpha*c0*H*mu_ij/rho_bar_ij,
        mu_ij =  H*(v_ij . r_ij)/(r_ij**2 + eps_visc*H**2)
        linearised (switch and quadratic beta*mu_ij**2 term dropped, see
        note below) into
            a_i,x^visc = alpha*c0*H**2 * sum_j Vj*dWr*rx**2/(r*(r**2+eps_visc*H**2)) * (ux_i-ux_j)
        i.e. weighted by rx**2 -- it is a projection onto the propagation
        direction, so this model is anisotropic even in its linearised form.
        alpha is reached from nu via the standard equivalence
            alpha = alpha_from_nu(nu, c0, H, dim) = 2*(dim+2)*nu/(c0*H)
        (Monaghan & Gingold 1983; Monaghan 1992, 2005.)

    OPTION B -- Morris, Fox & Zhu (1997) physical-viscosity Laplacian
        Directly discretises nu*Laplacian(u) as
            a_i,x^visc = 2*nu * sum_j Vj*dWr*r/(r**2+eps_visc*H**2) * (ux_i-ux_j)
        i.e. weighted only by the distance r (isotropic -- no directional
        projection), and uses nu directly with no alpha/c0/H conversion.
        This is the more standard choice when nu is meant as a genuine
        physical kinematic viscosity rather than a shock-capturing device.

    New parameters
    ----------
    nu : float
        Kinematic viscosity [m^2/s]. nu=0 recovers the original inviscid
        spectrum exactly under either option.
    dim : int
        Spatial dimension, used only by OPTION A's nu -> alpha conversion
        (default 2, matching these 2D lattices).
    eps_visc : float
        Regularisation constant in the mu_ij / Laplacian denominator
        (standard value 0.01).

    Linearisation note (applies to OPTION A)
    ----------
    The true Monaghan-Gingold Pi_ij is switched off unless v_ij.r_ij < 0
    (approaching particles) and also carries a quadratic beta*mu_ij**2
    term. Both are nonlinear in the velocity perturbation about the
    quiescent background (u=0), so they cannot enter a constant-coefficient
    linear operator. Following standard practice in the delta-SPH
    linear-stability literature, we keep only the always-on,
    linear-in-mu_ij piece -- the known equivalence between "always-on"
    alpha-viscosity and an effective Laplacian/physical viscosity. OPTION B
    has no such issue: the Morris form is already a linear Laplacian
    discretisation, with no switch or quadratic term to drop.

    Returns
    ----------
    G_hat, D_hat, PI_hat, mu_pm
        NOTE: return signature grew from 3 to 4 values (PI_hat added).
        Existing call sites unpacking `G_hat, D_hat, mu_pm = compute_spectrum(...)`
        must be updated to `G_hat, D_hat, PI_hat, mu_pm = compute_spectrum(...)`.
    """
    _W  = W_fn if W_fn is not None else W_vec
    _dW = dWdr_fn if dWdr_fn is not None else dWdr_vec
    Wr  = _W(r, H)
    dWr = _dW(r, H)
    dWx = dWr * rx / r

    # ==================================================================
    # VISCOSITY MODEL SELECTION -- uncomment exactly ONE of the two blocks
    # ==================================================================

    # --- OPTION A: Monaghan-Gingold artificial viscosity (anisotropic) ---
    alpha = alpha_from_nu(nu, c0, H, dim)
    visc_w         = dWr * rx**2 / (r * (r**2 + eps_visc*H**2))
    visc_prefactor = alpha*c0*H**2

    # --- OPTION B: Morris (1997) physical-viscosity Laplacian (isotropic) ---
    # visc_w         = dWr * r / (r**2 + eps_visc*H**2)
    # visc_prefactor = 2*nu

    # This is a difference (Laplacian-type) quantity, just like D_hat, so it
    # needs the same (1-phase) treatment (a raw undifferenced sum would not
    # vanish for a uniform velocity field, which is unphysical for a viscous term).
    G_hat  = np.zeros(len(xi_vals), dtype=complex)
    D_hat  = np.zeros(len(xi_vals), dtype=complex)
    PI_hat = np.zeros(len(xi_vals), dtype=complex)
    for idx, xi in enumerate(xi_vals):
        k = xi / h
        phase = np.exp(1j*k*rx)
        G_hat[idx]  = np.sum(Vj*dWx*phase)
        D_hat[idx]  = np.sum(Vj*(Wr/r**2)*(1 - phase))
        PI_hat[idx] = np.sum(Vj*visc_w*(1 - phase))
    PI_hat *= visc_prefactor

    mu0_rho  = -delta*H*c0*D_hat      # density-diffusion trace term (unchanged)
    mu0_visc = PI_hat                  # viscous trace term added to the velocity eq.

    trace = mu0_rho + mu0_visc
    disc  = ((mu0_rho - mu0_visc)/2)**2 - c0**2*np.abs(G_hat)**2 + 0j
    sq    = np.sqrt(disc)
    mu_pm = trace/2 + np.array([1, -1])[:, None]*sq
    return G_hat, D_hat, PI_hat, mu_pm

def compute_spectrum_general(rx, ry, r, Vj, H, h, xi_vals, theta, c0, delta,
                              nu=0.0, dim=2, eps_visc=0.01, W_fn=None, dWdr_fn=None):
    """
    Same as compute_spectrum, but for a wave vector k = (k cos(theta), k sin(theta))
    instead of being restricted to k = (k, 0). Now includes the standard
    Monaghan-Gingold artificial viscosity (see compute_spectrum's docstring
    for the exact linearisation used -- switch and beta*mu**2 term dropped).

    theta : propagation angle of the wave, in radians (theta=0 reproduces
            the original x-only compute_spectrum's inviscid part exactly).
    nu : float
        Kinematic viscosity [m^2/s]; converted to alpha via
        alpha_from_nu(nu, c0, H, dim) for OPTION A only (see below).
        nu=0 disables viscosity under either option.
    dim : int
        Spatial dimension used in the nu -> alpha conversion (OPTION A only,
        default 2).
    eps_visc : float
        Regularisation constant in the Monaghan-Gingold mu_ij denominator.

    Two viscosity models are available -- choose ONE by commenting/
    uncommenting the corresponding block below (see compute_spectrum's
    docstring for the full derivation of each):
      OPTION A: Monaghan-Gingold (anisotropic tensor PIxx/PIxy/PIyy, alpha
                reached from nu via alpha_from_nu). Active by default.
      OPTION B: Morris (1997) physical-viscosity Laplacian (isotropic,
                PIxy=0, uses nu directly).

    WHY THE OLD 2-EIGENVALUE FORMULA WAS WRONG (and couldn't just take
    viscosity as an extra term)
    ----------
    The previous version collapsed the 2D velocity field to a single scalar
    by projecting the kernel-gradient sum onto the propagation direction
    (G_hat = kdir_x*Gx_hat + kdir_y*Gy_hat) before forming the (rho, u)
    quadratic -- this only reduces to the true 3-DOF (rho, ux, uy) physics
    when the kernel sum is isotropic, which is generally false on a finite
    lattice for theta != 0. Viscosity in particular acts differently on the
    longitudinal and transverse (shear) velocity components, so there is no
    slot for it in a formula that has already thrown the transverse
    component away -- this is exactly why it was flagged broken instead of
    silently extended.

    Fix: assemble the genuine 3x3 linear operator for the state
    (rho', ux', uy') at each wavevector k = (xi/h)*(cos theta, sin theta),

        d/dt [rho'; ux'; uy'] = A(k) [rho'; ux'; uy']

        A(k) = [[ -delta*H*c0*D_hat,     -rho0*Gx_hat,        -rho0*Gy_hat      ],
                [ -(c0**2/rho0)*Gx_hat,   PIxx_hat,            PIxy_hat         ],
                [ -(c0**2/rho0)*Gy_hat,   PIxy_hat,            PIyy_hat         ]]

    and diagonalise it -- exactly the same physics as get_spectrum_direct,
    just evaluated in Fourier space instead of built from a finite periodic
    particle set. Since A(k) is only 3x3, its eigenvalues are the closed-form
    (Cardano) roots of a cubic, computed via solve_cubic_vectorized rather
    than a per-k np.linalg.eigvals call -- fully vectorised over xi_vals.
    rho0 only rescales the eigenvectors, not the eigenvalues, so it is fixed
    at 1.0 here (as in a standard von Neumann analysis); pass a different
    rho0 if needed.

    Returns
    ----------
    Gx_hat, Gy_hat, D_hat, PIxx_hat, PIxy_hat, PIyy_hat, mu_eigs
        mu_eigs has shape (len(xi_vals), 3): the three eigenvalues of A(k)
        at each xi. NOTE: this is a signature change from the previous
        (broken) 3-value return (G_hat, D_hat, mu_pm) -- update call sites.
    """
    alpha = alpha_from_nu(nu, c0, H, dim)
    rho0 = 1.0  # eigenvalues of A(k) are independent of rho0 (see docstring)

    _W  = W_fn if W_fn is not None else W_vec
    _dW = dWdr_fn if dWdr_fn is not None else dWdr_vec
    Wr  = _W(r, H)
    dWr = _dW(r, H)
    fx = dWr * rx / r        # x-component of kernel gradient
    fy = dWr * ry / r        # y-component of kernel gradient

    # ==================================================================
    # VISCOSITY MODEL SELECTION -- uncomment exactly ONE of the two blocks
    # (see compute_spectrum's docstring for the full derivation of each)
    # ==================================================================

    # --- OPTION A: Monaghan-Gingold artificial viscosity (anisotropic tensor) ---
    visc_base      = dWr / (r * (r**2 + eps_visc*H**2))
    wxx, wxy, wyy  = visc_base*rx**2, visc_base*rx*ry, visc_base*ry**2
    visc_prefactor = alpha*c0*H**2

    # --- OPTION B: Morris (1997) physical-viscosity Laplacian (isotropic) ---
    # visc_base      = dWr * r / (r**2 + eps_visc*H**2)
    # wxx, wxy, wyy  = visc_base, 0.0*visc_base, visc_base   # no ux-uy shear coupling
    # visc_prefactor = 2*nu

    kdir_x, kdir_y = np.cos(theta), np.sin(theta)

    Gx_hat   = np.zeros(len(xi_vals), dtype=complex)
    Gy_hat   = np.zeros(len(xi_vals), dtype=complex)
    D_hat    = np.zeros(len(xi_vals), dtype=complex)
    PIxx_hat = np.zeros(len(xi_vals), dtype=complex)
    PIxy_hat = np.zeros(len(xi_vals), dtype=complex)
    PIyy_hat = np.zeros(len(xi_vals), dtype=complex)
    for idx, xi in enumerate(xi_vals):
        k = xi / h
        kx, ky = k*kdir_x, k*kdir_y
        phase = np.exp(1j*(kx*rx + ky*ry))
        Gx_hat[idx]   = np.sum(Vj*fx*phase)
        Gy_hat[idx]   = np.sum(Vj*fy*phase)
        D_hat[idx]    = np.sum(Vj*(Wr/r**2)*(1 - phase))
        PIxx_hat[idx] = np.sum(Vj*wxx*(1 - phase))
        PIxy_hat[idx] = np.sum(Vj*wxy*(1 - phase))
        PIyy_hat[idx] = np.sum(Vj*wyy*(1 - phase))
    PIxx_hat *= visc_prefactor
    PIxy_hat *= visc_prefactor
    PIyy_hat *= visc_prefactor

    mu0_rho = -delta*H*c0*D_hat

    # Characteristic polynomial of A(k) = lambda**3 - tr(A)*lambda**2 + M*lambda - det(A) = 0,
    # i.e. lambda**3 + p2*lambda**2 + p1*lambda + p0 = 0 with:
    a, b, c = mu0_rho, -rho0*Gx_hat, -rho0*Gy_hat
    d, e, f = -(c0**2/rho0)*Gx_hat, PIxx_hat, PIxy_hat
    g, h    = -(c0**2/rho0)*Gy_hat, PIyy_hat
    tr_A  = a + e + h
    M_A   = (a*e - b*d) + (a*h - c*g) + (e*h - f*f)
    det_A = a*(e*h - f*f) - b*(d*h - f*g) + c*(d*f - e*g)

    mu_eigs = solve_cubic_vectorized(-tr_A, M_A, -det_A)

    return Gx_hat, Gy_hat, D_hat, PIxx_hat, PIxy_hat, PIyy_hat, mu_eigs


def get_spectrum_direct(x, y, L, Vj, H, h, c0, delta, rho0=1000.0, n_images=2,
                        nu=0.0, dim=2, eps_visc=0.01,
                        W_fn=None, dWdr_fn=None, support_radius=None):

    """
    Directly assemble the linear stability operator from particle-particle
    interactions (periodic box L x L, n_images periodic shells for neighbour
    search) and return eigenvalues via np.linalg.eigvals. rho0 is taken as a
    keyword argument (default 1000.0) so callers -- e.g. the interactive GUI
    -- can control it explicitly.

    New parameters
    ----------
    nu : float
        Kinematic viscosity [m^2/s]. nu=0 recovers the original inviscid
        operator exactly under either viscosity model below.
    dim : int
        Spatial dimension, used only by OPTION A's nu -> alpha conversion
        (default 2).
    eps_visc : float
        Regularisation constant in the viscous denominator (standard 0.01).

    Two viscosity models are available -- choose ONE by commenting/
    uncommenting the corresponding block below (see compute_spectrum's
    docstring for the full derivation of each):
      OPTION A: Monaghan-Gingold artificial viscosity (active by default).
                Genuinely couples ux and uy: since v_ij.r_ij mixes both
                velocity components, the linearised term fills the full
                off-diagonal blocks L_ux_uy / L_uy_ux, not just the
                diagonal L_ux_ux / L_uy_uy. alpha reached from nu via
                alpha_from_nu(nu, c0, H, dim).
      OPTION B: Morris (1997) physical-viscosity Laplacian. Isotropic --
                acts on ux and uy independently, so L_ux_uy = L_uy_ux = 0.
                Uses nu directly, no alpha/c0/H conversion.
    """
    alpha = alpha_from_nu(nu, c0, H, dim)

    _W  = W_fn if W_fn is not None else W_vec
    _dW = dWdr_fn if dWdr_fn is not None else dWdr_vec
    _R  = support_radius if support_radius is not None else 2*H
    N = len(x)
    Gx = np.zeros((N, N))
    Gy = np.zeros((N, N))
    D  = np.zeros((N, N))
    # Viscous coupling matrices (built alongside Gx,Gy,D). Pxx/Pxy/Pyy feed
    # OPTION A (Monaghan-Gingold); Piso feeds OPTION B (Morris).
    Pxx  = np.zeros((N, N))
    Pxy  = np.zeros((N, N))
    Pyy  = np.zeros((N, N))
    Piso = np.zeros((N, N))

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

        # --- viscosity, always-on linearisation (both models computed; pick one post-loop) ---
        # OPTION A (Monaghan-Gingold) pairwise weight: alpha*c0*H**2*Vj*dWr/(r*(r**2+eps_visc*H**2)),
        # coupling matrices Pxx=w*rx**2, Pxy=w*rx*ry, Pyy=w*ry**2 (anisotropic tensor).
        visc_k = alpha*c0*H**2 * dWr / (r * (r**2 + eps_visc*H**2))
        np.add.at(Pxx[i], jidx, Vjr*visc_k*rx**2)
        np.add.at(Pxy[i], jidx, Vjr*visc_k*rx*ry)
        np.add.at(Pyy[i], jidx, Vjr*visc_k*ry**2)

        # OPTION B (Morris 1997) pairwise weight: 2*nu*Vj*dWr*r/(r**2+eps_visc*H**2)
        # (isotropic -- no rx/ry projection, so no ux-uy shear coupling).
        visc_k_iso = 2*nu * dWr * r / (r**2 + eps_visc*H**2)
        np.add.at(Piso[i], jidx, Vjr*visc_k_iso)

    diagD = np.diag(D.sum(axis=1))
    L_rho_rho = delta*H*c0*(D - diagD)
    L_rho_ux  = -rho0*Gx
    L_rho_uy  = -rho0*Gy
    L_ux_rho  = -(c0**2/rho0)*Gx
    L_uy_rho  = -(c0**2/rho0)*Gy

    # ==================================================================
    # VISCOSITY MODEL SELECTION -- uncomment exactly ONE of the two blocks
    # ==================================================================
    # Note the acceleration is a_i = w_ij*(u_i - u_j), i.e. row i picks up
    # +rowsum on the diagonal and -w_ij off-diagonal -- the OPPOSITE sign
    # convention to (D - diagD) above (which represents sum_j D_ij*(rho_j-rho_i)).
    # Hence diag(rowsum) - P here, not P - diag(rowsum).

    # --- OPTION A: Monaghan-Gingold (anisotropic, ux-uy coupled) ---
    L_ux_ux = np.diag(Pxx.sum(axis=1)) - Pxx
    L_ux_uy = np.diag(Pxy.sum(axis=1)) - Pxy
    L_uy_ux = np.diag(Pxy.sum(axis=1)) - Pxy   # Pxy is symmetric (rx*ry) by construction
    L_uy_uy = np.diag(Pyy.sum(axis=1)) - Pyy

    # --- OPTION B: Morris (1997) physical viscosity (isotropic, no ux-uy coupling) ---
    # L_ux_ux = np.diag(Piso.sum(axis=1)) - Piso
    # L_ux_uy = np.zeros((N, N))
    # L_uy_ux = np.zeros((N, N))
    # L_uy_uy = np.diag(Piso.sum(axis=1)) - Piso

    #diagD = np.diag(D.sum(axis=1))
    #L_rho_rho = delta * H * c0 * (D - diagD)

    #Sx = np.diag(Gx.sum(axis=1))
    #Sy = np.diag(Gy.sum(axis=1))

    #L_rho_ux = -rho0 * (Gx - Sx)           # was: -rho0*Gx
    #L_rho_uy = -rho0 * (Gy - Sy)           # was: -rho0*Gy
    #L_ux_rho = -(c0**2/rho0) * (Gx + Sx)   # note: PLUS Sx here, not minus
    #L_uy_rho = -(c0**2/rho0) * (Gy + Sy)

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
    color_shear: str = 'mediumslateblue'
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
    """Generalised Fourier spectra (with viscosity) for cartesian / hex / random at params.theta."""
    kinfo = KERNELS[params.kernel]
    kW, kdw, kR = kinfo["W_vec"], kinfo["dWdr_vec"], kinfo["support"] * params.H
    out: dict = {}
    rxC, ryC, rC, VjC = neighbours_cartesian(params.h, params.H, support_radius=kR)
    Gx, _Gy, D, _PIxx, _PIxy, _PIyy, mu = compute_spectrum_general(
        rxC, ryC, rC, VjC, params.H, params.h, params.xi_vals,
        params.theta, params.c0, params.delta, params.nu, params.dim, params.eps_visc, kW, kdw)
    out['cartesian'] = (Gx, D, mu.T)
    rxH, ryH, rH, VjH = neighbours_hex(params.h, params.H, support_radius=kR)
    Gx, _Gy, D, _PIxx, _PIxy, _PIyy, mu = compute_spectrum_general(
        rxH, ryH, rH, VjH, params.H, params.h, params.xi_vals,
        params.theta, params.c0, params.delta, params.nu, params.dim, params.eps_visc, kW, kdw)
    out['hex'] = (Gx, D, mu.T)
    rnd = []
    for s in range(params.n_real):
        rxR, ryR, rR, VjR = neighbours_random(params.h, params.H, seed=s, support_radius=kR)
        Gx, _Gy, D, _PIxx, _PIxy, _PIyy, mu = compute_spectrum_general(
            rxR, ryR, rR, VjR, params.H, params.h, params.xi_vals,
            params.theta, params.c0, params.delta, params.nu, params.dim, params.eps_visc, kW, kdw)
        rnd.append((Gx, D, mu.T))
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
    colors = [color_plus, color_minus, config.color_shear]
    for i in range(len(mu_pm)):
        ax.scatter(mu_pm[i].real, mu_pm[i].imag, s=s, c=colors[i % len(colors)])
    ax.axhline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.axvline(0, color=config.zero_line_color, lw=config.zero_line_lw)
    ax.set_xlabel(config.xlabel)
    ax.set_ylabel(config.ylabel)
    ax.set_title(title)


def plot_mu_panel_random(ax, mu_list, config: PlotConfig, title):
    for _G, _D, mu in mu_list:
        for i in range(len(mu)):
            ax.scatter(mu[i].real, mu[i].imag, s=config.scatter_size_random,
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
    _Gx_c, _Gy_c, _D_c, _PIxx_c, _PIxy_c, _PIyy_c, mu_eigs_c = compute_spectrum_general(
        rxC, ryC, rC, VjC, p.H, p.h, p.xi_vals, np.pi/2, p.c0, p.delta, kW, kdw)
    mu_c = mu_eigs_c.T
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

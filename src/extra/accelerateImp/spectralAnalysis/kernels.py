import numpy as np

#Wendland C2 (Dehnen & Aly 2012)
def W_wendlandC2(r, H):
    q = r / H
    if q >= 2: return 0.0
    return (7/(4*np.pi*H**2)) * (1 - q/2)**4 * (2*q + 1)

def dWdr_wendlandC2(r, H):
    q = r / H
    if q >= 2: return 0.0
    return (7/(4*np.pi*H**2)) * (-2*(1-q/2)**3*(2*q+1) + 2*(1-q/2)**4) / H

#Wendland C4
def W_wendlandC4(r, H):
    q = r / H
    if q >= 2: return 0.0
    sigma = 9.0 / (4.0*np.pi*H**2)
    return sigma * (1 - q/2)**6 * (35/12*q**2 + 3*q + 1)

def dWdr_wendlandC4(r, H):
    q = r / H
    if q >= 2: return 0.0
    sigma = 9.0 / (4.0*np.pi*H**2)
    poly  = 35/12*q**2 + 3*q + 1
    dpoly = 35/6*q + 3
    return sigma * ( -3*(1 - q/2)**5*poly + (1 - q/2)**6*dpoly ) / H

# Cubic spline (Monaghan 1992, "M4")
def W_cubic(r, H):
    q = r / H
    sigma = 10.0 / (7.0*np.pi*H**2)
    if q >= 2: return 0.0
    if q < 1:
        return sigma * (1 - 1.5*q**2 + 0.75*q**3)
    else:
        return sigma * 0.25 * (2 - q)**3

def dWdr_cubic(r, H):
    q = r / H
    sigma = 10.0 / (7.0*np.pi*H**2)
    if q >= 2: return 0.0
    if q < 1:
        dWdq = -3*q + 2.25*q**2
    else:
        dWdq = -0.75 * (2 - q)**2
    return sigma * dWdq / H

# Quartic spline ("M5")
def W_quartic(r, H):
    q = r / H
    sigma = 96.0 / (1199.0*np.pi*H**2)
    if q >= 2.5: return 0.0
    val = (2.5 - q)**4
    if q < 1.5: val -= 5*(1.5 - q)**4
    if q < 0.5: val += 10*(0.5 - q)**4
    return sigma * val

def dWdr_quartic(r, H):
    q = r / H
    sigma = 96.0 / (1199.0*np.pi*H**2)
    if q >= 2.5: return 0.0
    dval = -4*(2.5 - q)**3
    if q < 1.5: dval += 20*(1.5 - q)**3
    if q < 0.5: dval -= 40*(0.5 - q)**3
    return sigma * dval / H

# Quintic spline ("M6", Morris 1996)
def W_quintic(r, H):
    q = r / H
    sigma = 7.0 / (478.0*np.pi*H**2)
    if q >= 3: return 0.0
    val = (3 - q)**5
    if q < 2: val -= 6*(2 - q)**5
    if q < 1: val += 15*(1 - q)**5
    return sigma * val

def dWdr_quintic(r, H):
    q = r / H
    sigma = 7.0 / (478.0*np.pi*H**2)
    if q >= 3: return 0.0
    dval = -5*(3 - q)**4
    if q < 2: dval += 30*(2 - q)**4
    if q < 1: dval -= 75*(1 - q)**4
    return sigma * dval / H

# Gaussian (truncated)
def W_gaussian(r, H, q_cut=3.0):
    q = r / H
    if q >= q_cut: return 0.0
    # normalization chosen so that 2*pi*int_0^{q_cut*H} W r dr = 1
    alpha = 1.0 / (np.pi*H**2 * (1 - np.exp(-q_cut**2)))
    return alpha * np.exp(-q**2)

def dWdr_gaussian(r, H, q_cut=3.0):
    q = r / H
    if q >= q_cut: return 0.0
    alpha = 1.0 / (np.pi*H**2 * (1 - np.exp(-q_cut**2)))
    return alpha * (-2*q) * np.exp(-q**2) / H

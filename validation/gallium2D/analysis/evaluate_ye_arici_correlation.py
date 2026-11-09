#!/usr/bin/env python3
"""
Evaluates the Ye & Arıcı (2023) benchmark correlation (Eq. 6) for 2D pure gallium melting:
  phi_n = 2.708 * tau^{0.712} * Ra^{0.0243} * Ar^{-0.0515}
from:
  W.-B. Ye and M. Arıcı, "3D validation, 2D feasibility, corrected and developed correlations
  for pure solid-gallium phase change modeling by enthalpy-porosity methodology",
  International Communications in Heat and Mass Transfer, Vol. 144, 106780 (2023).
  DOI: 10.1016/j.icheatmasstransfer.2023.106780

where:
  tau = Fo * Ste
  Fo = alpha * t / W^2
  Ste = Cp * (T_hot - T_m) / L
  Ra = g * beta * (T_hot - T_m) * H^3 / (nu * alpha)
  Ar = H / W
Also evaluates sensitivity to the Rayleigh number exponent (0.0243).
While direct sensitivity to further Ra variation is weak (Ra^0.0243),
the exponent tau^0.712 and prefactor 2.708 embody substantial natural
convection enhancement over pure conduction (~0.57 vs ~0.37).
"""

import numpy as np

def evaluate_ye_arici(t_seconds=1020.0):
    # Geometry (Brent et al. 1988 cavity)
    W = 0.0889   # 8.89 cm
    H = 0.0635   # 6.35 cm
    Ar = H / W   # 0.714286

    # Physical properties of gallium
    rho = 6095.0
    Cp = 381.5
    k = 32.0
    mu = 1.81e-3
    nu = mu / rho
    alpha = k / (rho * Cp)
    beta = 1.2e-4
    g = 9.81

    # Temperatures & Latent Heat
    T_hot = 311.15
    T_m = 302.93
    dT = T_hot - T_m
    L = 80160.0

    # Dimensionless numbers
    Ste = Cp * dT / L
    Fo = alpha * t_seconds / (W**2)
    tau = Fo * Ste
    Ra = g * beta * dT * (H**3) / (nu * alpha)

    # Correlation Eq. 6:
    # phi_n = 2.708 * tau^0.712 * Ra^0.0243 * Ar^-0.0515
    term_tau = tau**0.712
    term_Ra = Ra**0.0243
    term_Ar = Ar**(-0.0515)
    phi_n = 2.708 * term_tau * term_Ra * term_Ar

    return {
        "t_sec": t_seconds,
        "t_min": t_seconds / 60.0,
        "W": W,
        "H": H,
        "Ar": Ar,
        "alpha": alpha,
        "nu": nu,
        "Ste": Ste,
        "Fo": Fo,
        "tau": tau,
        "Ra": Ra,
        "term_tau": term_tau,
        "term_Ra": term_Ra,
        "term_Ar": term_Ar,
        "phi_n": phi_n
    }

if __name__ == "__main__":
    res = evaluate_ye_arici(1020.0)
    print("=" * 70)
    print(" YE & ARICI (2023) BENCHMARK CORRELATION EVALUATION (Eq. 6)")
    print("=" * 70)
    print(f"Time:                       {res['t_sec']:.1f} s ({res['t_min']:.1f} min)")
    print(f"Cavity Width (W):           {res['W']*100:.2f} cm")
    print(f"Cavity Height (H):          {res['H']*100:.2f} cm")
    print(f"Aspect Ratio (Ar = H/W):    {res['Ar']:.6f}")
    print(f"Thermal Diffusivity alpha:  {res['alpha']:.6e} m^2/s")
    print(f"Kinematic Viscosity nu:     {res['nu']:.6e} m^2/s")
    print(f"Stefan Number (Ste):        {res['Ste']:.5f}")
    print(f"Fourier Number (Fo):        {res['Fo']:.5f}")
    print(f"Dimensionless Time tau:     {res['tau']:.6f}")
    print(f"Rayleigh Number (Ra):       {res['Ra']:.3e}")
    print("-" * 70)
    print(f"Component: 2.708            = 2.7080")
    print(f"Component: tau^0.712        = {res['term_tau']:.6f}")
    print(f"Component: Ra^0.0243        = {res['term_Ra']:.6f}")
    print(f"Component: Ar^-0.0515       = {res['term_Ar']:.6f}")
    print("-" * 70)
    print(f"Evaluated Liquid Fraction:  phi_n = {res['phi_n']:.4f} ({res['phi_n']:.6f})")
    print("=" * 70)
    print("\nPhysical Note on Scaling:")
    print("Because the Rayleigh exponent is only 0.0243, varying Ra by an order of magnitude (10x)")
    print(f"alters the Ra factor by only 10^0.0243 = {10**0.0243:.4f} (+5.7%).")
    print("However, pure 1D Stefan conduction yields f_L ~ 0.37 at 17 min vs 0.5703 from the correlation (+54%),")
    print("confirming that tau^0.712 and prefactor 2.708 embody substantial convective enhancement over conduction.")

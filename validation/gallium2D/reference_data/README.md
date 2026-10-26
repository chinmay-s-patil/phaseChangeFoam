# Reference Benchmark Data: Brent, Voller, Reid (1988)

This directory contains digitized benchmark data extracted directly from:

> **A. D. Brent, V. R. Voller, and K. J. Reid (1988)**  
> *"Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal"*  
> **Numerical Heat Transfer**, 13(3), pp. 297–318.  
> DOI: [10.1080/10407788808913615](http://dx.doi.org/10.1080/10407788808913615)  
> Paper copy: `validation/papers/brent1988.pdf`

---

## Benchmark Parameters (Table 1 & Section 4)

| Parameter | Symbol | Value | Unit |
| :--- | :--- | :--- | :--- |
| Cavity Width | $W$ | $0.0889$ (8.89 cm) | $\text{m}$ |
| Cavity Height | $H$ | $0.0635$ (6.35 cm) | $\text{m}$ |
| Aspect Ratio | $A = W / H$ | $1.40$ | - |
| Hot Wall Temperature (Left, $x=0$) | $T_{hot}$ | $38.0$ (311.15 K) | $^\circ\text{C}$ |
| Cold Wall Temperature (Right, $x=W$) | $T_{cold}$ | $28.3$ (301.45 K) | $^\circ\text{C}$ |
| Initial Temperature | $T_{init}$ | $28.3$ (301.45 K) | $^\circ\text{C}$ |
| Melting Temperature | $T_m$ | $29.78$ (302.93 K) | $^\circ\text{C}$ |
| Latent Heat of Fusion | $L$ | $80160$ | $\text{J/kg}$ |
| Specific Heat Capacity | $c_p$ | $381.5$ | $\text{J/(kg}\cdot\text{K)}$ |
| Thermal Conductivity | $k$ | $32.0$ | $\text{W/(m}\cdot\text{K)}$ |
| Dynamic Viscosity | $\mu$ | $1.81 \times 10^{-3}$ | $\text{Pa}\cdot\text{s}$ |
| Liquid Density | $\rho$ | $6093$ | $\text{kg/m}^3$ |
| Reference Density | $\rho_0$ | $6095$ | $\text{kg/m}^3$ |
| Thermal Expansion Coefficient | $\beta$ | $1.2 \times 10^{-4}$ | $\text{K}^{-1}$ |
| Darcy Morphology Constant | $C$ | $1.6 \times 10^6$ | $\text{kg/(m}^3\cdot\text{s)}$ |
| Darcy Porosity Denominator Constant | $q$ | $1.0 \times 10^{-3}$ | - |
| Prandtl Number | $\text{Pr}$ | $0.0216$ | - |
| Stefan Number | $\text{Ste}$ | $0.039$ | - |
| Rayleigh Number | $\text{Ra}$ | $\approx 6 \times 10^5$ | - |

---

## Digitized Files

1. **`brent1988_front_2min.csv`**: Numerical melt front position at $t = 2\text{ min}$ (120 s) digitized from Figure 9 (solid line).
2. **`brent1988_front_6min.csv`**: Numerical melt front position at $t = 6\text{ min}$ (360 s) digitized from Figure 9 (solid line).
3. **`brent1988_front_10min.csv`**: Numerical melt front position at $t = 10\text{ min}$ (600 s) digitized from Figure 9 (solid line).
4. **`brent1988_front_17min.csv`**: Numerical melt front position at $t = 17\text{ min}$ (1020 s) digitized from Figure 9 (solid line).
5. **`brent1988_liquid_fraction.csv`**: Liquid phase fraction $f_L = \frac{1}{W H}\int x_{front}(y)\,dy$ vs time.

Columns in front CSVs:
- `x_m`: Horizontal position of front in meters
- `y_m`: Vertical position of front in meters
- `x_norm`: Normalized horizontal position $x / W \in [0, 1]$
- `y_norm`: Normalized vertical position $y / H \in [0, 1]$

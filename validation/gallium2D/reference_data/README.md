# Reference Benchmark Data: Brent, Voller, Reid (1988)

This directory contains digitized benchmark data extracted directly from:

> **A. D. Brent, V. R. Voller, and K. J. Reid (1988)**  
> *"Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal"*  
> **Numerical Heat Transfer**, 13(3), pp. 297–318.  
> DOI: [10.1080/10407788808913615](https://doi.org/10.1080/10407788808913615)  

> These files contain independently digitised approximate coordinates extracted from the published figures. They should not be interpreted as the original experimental or machine-readable numerical datasets.

---

## Digitized Files

All numerical melt front curves were traced directly from Figure 9 (solid lines) of Brent, Voller, and Reid (1988). The digitised horizontal front positions carry an estimated tracing uncertainty of $\pm 1\text{ mm}$ (approx. $\pm 1.1\%$ of cavity width $W = 88.9\text{ mm}$), attributable to line thickness and resolution in the published figure:

1. **`brent1988_front_2min.csv`**: Numerical melt front position at $t = 2\text{ min}$ (120 s) traced from Figure 9 ($\pm 1\text{ mm}$ horizontal position uncertainty).
2. **`brent1988_front_6min.csv`**: Numerical melt front position at $t = 6\text{ min}$ (360 s) traced from Figure 9 ($\pm 1\text{ mm}$ horizontal position uncertainty).
3. **`brent1988_front_10min.csv`**: Numerical melt front position at $t = 10\text{ min}$ (600 s) traced from Figure 9 ($\pm 1\text{ mm}$ horizontal position uncertainty).
4. **`brent1988_front_17min.csv`**: Numerical melt front position at $t = 17\text{ min}$ (1020 s) traced from Figure 9 ($\pm 1\text{ mm}$ horizontal position uncertainty).
5. **`brent1988_liquid_fraction.csv`**: Liquid phase fraction $f_L = \frac{1}{W H}\int x_{front}(y)\,dy$ vs time digitized from Figure 7 of Brent et al. (1988).

Columns in front CSVs:
- `x_m`: Horizontal position of front in meters
- `y_m`: Vertical position of front in meters
- `x_norm`: Normalized horizontal position $x / W \in [0, 1]$
- `y_norm`: Normalized vertical position $y / H \in [0, 1]$

---

## Attribution & Legal Notice

The digitised coordinate dataset is copyright © 2026 Chinmay S Patil and licensed under GPL-3.0-or-later for benchmark reproduction and regression testing. No license or copyright is claimed or granted on the original published figures, plots, or underlying article contents of Brent, Voller, and Reid (1988) / Taylor & Francis.

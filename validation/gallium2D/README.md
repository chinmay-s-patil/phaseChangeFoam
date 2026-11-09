# 2D Pure Gallium Melting Benchmark Verification

Numerical benchmark verification suite for 2D natural convection melting of pure gallium in a rectangular cavity ($W = 88.9\text{ mm}$, $H = 63.5\text{ mm}$, Aspect Ratio $A = 0.714$) heated from the vertical side.

Solvers: `phaseChangeMultiRegionFoam` (OpenFOAM-v2412).

---

## 1. Physical & Numerical Benchmark References

- **Brent, Voller, and Reid (1988)**: "Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal", *Numerical Heat Transfer*, 13(3), pp. 297–318. [DOI: 10.1080/10407788808913615](https://doi.org/10.1080/10407788808913615).
- **Hannoun, Alexiades, and Mai (2003)**: "Resolving the Controversy over Tin and Gallium Melting in a Rectangular Cavity Heated from the Side", *Numerical Heat Transfer, Part B: Fundamentals*, 44(3), pp. 253–276. [DOI: 10.1080/10407790390226974](https://doi.org/10.1080/10407790390226974).
- **Ye and Arıcı (2023)**: "Redefined interface error, 2D verification and validation for pure solid-gallium phase change modeling by enthalpy-porosity methodology", *Int. Commun. Heat Mass Transfer*, 147, 106952. [DOI: 10.1016/j.icheatmasstransfer.2023.106952](https://doi.org/10.1016/j.icheatmasstransfer.2023.106952).
- **Gau and Viskanta (1986)**: "Melting and Solidification of a Pure Metal on a Vertical Wall", *J. Heat Transfer*, 108(1), pp. 174–181.

---

## 2. Numerical Configuration

The baseline benchmark uses the `EHCDarcy` formulation with a $42 \times 32$ structured mesh and a time step of $\Delta t = 0.1\text{ s}$.

For this benchmark configuration, `EHCDarcy` accounts for sensible heat-capacity variation through the path-based effective heat capacity with linearised latent heat sources, achieving discrete energy closure ($|\epsilon_{\text{energy}}| \le 0.38\%$) within the specified regression tolerance using three outer correctors ($n_{\text{outer}} = 3$).

The alternative `enthalpyPorosity` formulation (Brent et al. 1988) updates liquid fraction via Newton projection. For this benchmark problem, under-relaxation with $\lambda = 0.2$ leads to incomplete latent heat addition per time step unless compensated by additional outer correctors, whereas the full Newton step ($\lambda = 1.0$) requires higher corrector counts ($n_{\text{outer}} \ge 12$) to achieve comparable energy conservation.

Detailed comparisons between `EHCDarcy` and `enthalpyPorosity`, including relaxation-factor sensitivity, backward-Euler energy quadrature, and computational cost, are documented in [ANALYSIS_REPORT.md](analysis/ANALYSIS_REPORT.md).

---

## 3. Benchmark Execution

The baseline comparative summary is committed in [`gallium_validation_results.txt`](gallium_validation_results.txt) and can be reproduced or extended using `runTest.py`:

### Quick Verification Run (2 Minutes, Baseline Grid $42 \times 32$)

```bash
# Run default EHCDarcy (passes all criteria, |energy closure| <= 0.38%, ~35 s runtime):
python3 runTest.py

# Run alternative enthalpyPorosity with full Newton step (lambda = 1.0):
python3 runTest.py --model enthalpyPorosity --lambda 1.0

# Run full 17-minute benchmark transient (~3.5--5 min runtime):
python3 runTest.py --model EHCDarcy --endTime 1020 --plot
```

`runTest.py` updates `gallium_validation_results.txt` upon completion.

### Authoritative Report Generation

All tables, metrics, and narrative numbers are dynamically compiled from committed simulation datasets:

```bash
python3 analysis/generate_report.py
```

Generated report:
- [ANALYSIS_REPORT.md](analysis/ANALYSIS_REPORT.md)

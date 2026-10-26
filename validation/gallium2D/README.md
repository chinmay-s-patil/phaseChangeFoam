# 2D Pure Gallium Melting Benchmark Verification

Numerical benchmark verification suite for 2D natural convection melting of pure gallium in a rectangular cavity ($W = 88.9\text{ mm}$, $H = 63.5\text{ mm}$, Aspect Ratio $A = 0.714$) heated from the vertical side.

Solvers: `phaseChangeMultiRegionFoam` (OpenFOAM-v2412).

---

## 1. Physical & Numerical Benchmark References

- **Brent, Voller, and Reid (1988)**: "Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal", *Numerical Heat Transfer*, 13(3), pp. 297–318. [DOI: 10.1080/10407788808913615](https://doi.org/10.1080/10407788808913615).
- **Hannoun, Alexiades, and Mai (2003)**: "Resolving the Controversy over Tin and Gallium Melting in a Rectangular Cavity Heated from the Side", *Numerical Heat Transfer, Part B: Fundamentals*, 44(3), pp. 253–276. [DOI: 10.1080/10407790390226974](https://doi.org/10.1080/10407790390226974).
- **Ye and Arıcı (2023)**: "Redefined-interface-error-based 2D verification and validation study on natural convection melting: Benchmark with pure Gallium", *Int. Commun. Heat Mass Transfer*, 147, 106952. [DOI: 10.1016/j.icheatmasstransfer.2023.106952](https://doi.org/10.1016/j.icheatmasstransfer.2023.106952).
- **Gau and Viskanta (1986)**: "Melting and Solidification of a Pure Metal on a Vertical Wall", *J. Heat Transfer*, 108(1), pp. 174–181.

---

## 2. Modeling Recommendations & Formulations

### Recommended Default: `EHCDarcy`

For fluid-region liquid metal melting with natural convection, **`EHCDarcy` is the documented default formulation**:

- **Formulation**: Path-integral effective heat capacity $C_p^{\text{eff}}(T)$ with Carman-Kozeny Darcy porosity damping.
- **Energy Conservation**: Directly incorporates latent heat into the left-hand matrix capacity, achieving discrete energy conservation ($|\epsilon_{\text{energy}}| \le 0.38\%$) with only $n_{\text{outer}} = 3$ ($15.6\text{ s}$ CPU time for 1,200 time steps), and $-0.13\%$ over the full 17-minute transient.
- **Monotonic Convergence**: Energy closure is strictly monotonic and bounded within $[-0.38\%, -0.34\%]$ for all $n_{\text{outer}} \ge 3$.
- **Efficiency**: $3.6\times$ cheaper in CPU execution time and requires $75\%$ fewer outer corrector passes than Enthalpy-Porosity to achieve comparable discrete conservation.

### Alternative: `enthalpyPorosity` (Brent et al. 1988)

Enthalpy-Porosity updates the liquid fraction explicitly after the sensible enthalpy solve via Newton projection:

$$
f_{\text{new}} = \min\left(\max\left(f + \lambda \frac{T - T_f}{L/C_p + \Delta T_{\text{window}}}, 0\right), 1\right)
$$

#### Practical Numerical Guidance:

1. **Newton Relaxation Factor ($\lambda$)**:
   - Always maintain code-default **`lambda 1.0`** (the exact Newton step for a linear melting window).
   - **Never use under-relaxation $\lambda = 0.2$** in standard production runs. Under $\lambda = 0.2$, each outer pass corrects only $20\%$ of the Newton increment, causing incomplete latent heat addition and massive non-conservation ($+7.35\%$ energy surplus at 2 min, $+7.49\%$ at 17 min with $n_{\text{outer}}=8$; and $+2.86\%$ even at $n_{\text{outer}}=12$).
2. **PIMPLE Outer Correctors ($n_{\text{outer}}$)**:
   - Under $\lambda = 1.0$, explicit latent heat updates exhibit damped, under-relaxed Newton oscillations: $+56.8\%$ surplus at $n=2$, $+11.3\%$ at $n=3$, overshooting to $-7.3\%$ deficit at $n=5$, crossing zero near $n=7\text{--}8$ ($-0.06\%$).
   - To reliably settle within $\le 0.3\text{--}0.5\%$ closure across all conditions, `enthalpyPorosity` requires **$n_{\text{outer}} \ge 10\text{--}12$** ($56.6\text{ s}$ CPU time for 14,400 outer passes).

---

## 3. Discrete Energy Conservation Audit & Quadrature

In OpenFOAM's finite-volume formulation with implicit backward-Euler time discretization (`Euler ddtScheme`), the discrete storage change satisfies:

$$
\Delta E_{\text{stored}}^n = \dot{Q}_{\text{net}}(t^n)\,\Delta t_n
$$

The exact discrete boundary-flux conjugate is the **right-endpoint Riemann sum**:

$$
E_{\text{net, in}}(t_N) = \sum_{n=1}^N \dot{Q}_{\text{net}}(t_n)\,\Delta t_n
$$

Using composite trapezoidal quadrature adds an artificial offset $+\frac{\dot{Q}_1 - \dot{Q}_N}{2}\Delta t \approx +3.6\text{ J}$ due to the intense early conduction spike ($\approx 81\text{ W}$ at $0.1\text{ s}$). Right-endpoint quadrature eliminates this artifact, matching OpenFOAM's discrete matrix equations exactly.

---

## 4. Benchmark Execution

### Quick Verification Run (2 Minutes, Baseline Grid $42 \times 32$)

```bash
# Run default EHCDarcy (passes all criteria, |energy closure| <= 0.38%):
python3 runTest.py

# Run alternative enthalpyPorosity with full Newton step (lambda = 1.0):
python3 runTest.py --model enthalpyPorosity --lambda 1.0

# Run with custom parameters:
python3 runTest.py --model EHCDarcy --endTime 120 --nx 42 --ny 32
```

### Authoritative Report Generation

All tables, metrics, and narrative numbers are dynamically compiled from committed simulation datasets:

```bash
python3 analysis/generate_report.py
```

Generated reports:

- [ANALYSIS_REPORT.md](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/validation/gallium2D/analysis/ANALYSIS_REPORT.md)
- Artifact: [gallium_validation_analysis_report.md](file:///home/lavender/.gemini/antigravity-ide/brain/6e10d676-3179-4ccd-8848-a30db5d9e81a/gallium_validation_analysis_report.md)

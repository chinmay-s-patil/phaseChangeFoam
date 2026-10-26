# Numerical Benchmark Verification Report: 2D Pure Gallium Melting Benchmark

**Solver Suite**: `phaseChangeMultiRegionFoam` (OpenFOAM-v2412)  
**Formulations Analyzed**:
1. `EHCDarcy`: Path-Integral Effective Heat Capacity with Carman-Kozeny Darcy Porosity Damping
2. `enthalpyPorosity`: Enthalpy-Porosity Formulation (Brent, Voller & Reid 1988, Voller & Prakash 1987)

**Primary Benchmarks**:
- **Brent, Voller, and Reid (1988)**: "Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal", *Numerical Heat Transfer*, 13(3), pp. 297–318. DOI: [10.1080/10407788808913615](https://doi.org/10.1080/10407788808913615).
- **Hannoun, Alexiades, and Mai (2003)**: "Resolving the Controversy over Tin and Gallium Melting in a Rectangular Cavity Heated from the Side", *Numerical Heat Transfer, Part B: Fundamentals*, 44(3), pp. 253–276. DOI: [10.1080/10407790390226974](https://doi.org/10.1080/10407790390226974).
- **Ye and Arıcı (2023)**: "Redefined-interface-error-based 2D verification and validation study on natural convection melting: Benchmark with pure Gallium", *International Communications in Heat and Mass Transfer*, 147, 106952. DOI: [10.1016/j.icheatmasstransfer.2023.106952](https://doi.org/10.1016/j.icheatmasstransfer.2023.106952).
- **Gau and Viskanta (1986)**: "Melting and Solidification of a Pure Metal on a Vertical Wall", *J. Heat Transfer*, 108(1), pp. 174–181.

---

## 1. Executive Summary & Acceptance Scorecard

### 1.1 Scope Definition & Acceptance Framework
To preserve scientific precision, this report is strictly presented as a **numerical benchmark verification study** (code-to-code verification against published 2D numerical solutions):
- **Numerical Verification Scope**: Quantitative phase-front deviations and global liquid fractions are evaluated against the authoritative 2D numerical benchmark of Brent, Voller, and Reid (1988) and the grid-converged numerical scaling correlation of Ye & Arıcı (2023, Eq. 6).
- **Experimental Distinction**: Physical experimental measurements from Gau & Viskanta (1986) provide essential qualitative physical context regarding experimental delay mechanisms (heating jacket thermal inertia, pour-out measurement disruption, and 3D endwall drag). No quantitative validation against physical measurements is claimed.
- **Formulation Bounds of Grid Study**: Formal ASME / Roache grid-convergence index (GCI) verification is established specifically for the `EHCDarcy` formulation; no formal grid-convergence study was conducted for `enthalpyPorosity`.

Study-specific acceptance thresholds are defined as standard **engineering tolerances**:
1. **Front Morphology (Mean)**: Mean front deviation across cavity height $\overline{\delta}_{\text{front}} \le 5.0\%$ of cavity width $W$ ($4.45\text{ mm}$), evaluated against the digitized **numerical benchmark front of Brent, Voller, and Reid (1988, Figure 9)** (digitization line-width uncertainty $\pm 1.0\text{ mm}$).
2. **Front Morphology (Max)**: Peak local front deviation $\delta_{\text{front, max}} \le 10.0\%$ of cavity width $W$ ($8.89\text{ mm}$), relative to Brent et al. (1988).
3. **Global Liquid Fraction ($f_L$)**: Absolute liquid fraction deviation $|f_L^{\text{sim}} - f_L^{\text{ref}}| \le 0.05$ ($5.0$ percentage points) relative to Brent et al. (1988).
4. **Energy Conservation Closure**: Global net boundary-to-storage energy closure discrepancy $|\epsilon_{\text{energy}}| \le 2.0\%$ across all evaluated times ($t \le 17\text{ min}$).

### 1.2 Multi-Model Performance Scorecard (Full 17-Minute Benchmark)

| Metric | Engineering Target | EP (Tutorial: $\lambda = 0.2, n_{\text{outer}}=8$) | EP ($\lambda = 0.2, n_{\text{outer}}=12$) | EP (Default: $\lambda = 1.0, n_{\text{outer}}=8$) | EHCDarcy (Baseline: $n_{\text{outer}}=8$) | Status (EP, $\lambda=1.0$) | Status (EHCDarcy) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Energy Closure ($t = 2\text{ min}$)** | $\le 2.0\%$ | +7.35% (+132.3 J) *(FAIL)* | +2.86% (+52.9 J) *(FAIL)* | **-0.06% (-1.2 J)** | **-0.34% (-6.4 J)** | **PASSED** | **PASSED** |
| **Energy Closure ($t = 6\text{ min}$)** | $\le 2.0\%$ | +7.37% (+260.7 J) *(FAIL)* | +2.85% (+103.2 J) *(FAIL)* | **+0.11% (+4.2 J)** | **-0.23% (-8.4 J)** | **PASSED** | **PASSED** |
| **Energy Closure ($t = 17\text{ min}$)** | $\le 2.0\%$ | +7.49% (+590.4 J) *(FAIL)* | +2.92% (+234.4 J) *(FAIL)* | **-0.10% (-8.0 J)** | **-0.13% (-10.3 J)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 2\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | 1.95 mm (2.19%) | 1.77 mm (1.99%) | 1.71 mm (1.93%) | **1.67 mm (1.88%)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 6\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | 2.60 mm (2.93%) | 2.15 mm (2.42%) | 1.98 mm (2.22%) | **1.87 mm (2.11%)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 10\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | 3.51 mm (3.95%) | 2.76 mm (3.11%) | 2.45 mm (2.76%) | **2.35 mm (2.65%)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 17\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | 5.43 mm (6.11%) *(FAIL)* | 4.23 mm (4.76%) | 3.52 mm (3.96%) | **3.46 mm (3.89%)** | **PASSED** | **PASSED** |
| **Max Front Dev. ($t = 17\text{ min}$)** | $\le 10.0\%$ $W$ ($8.89\text{ mm}$) | 9.74 mm (10.95%) *(FAIL)* | 7.70 mm (8.66%) | 6.40 mm (7.20%) | **6.34 mm (7.13%)** | **PASSED** | **PASSED** |
| **$f_L$ at $t = 17\text{ min}$** (Brent: 0.5334) | $\pm 0.05$ | 0.5967 (+0.0633) *(FAIL)* | 0.5823 (+0.0489) | **0.5745 (+0.0411)** | **0.5736 (+0.0402)** | **PASSED** | **PASSED** |
| **Ye & Arıcı (2023) Eq. 6** (0.5703) | $\pm 0.03$ | +0.0264 (+2.64 pts) | +0.0120 (+1.20 pts) | **+0.0042 (+0.42 pts)** | **+0.0033 (+0.33 pts)** | **PASSED** | **PASSED** |

*Scorecard Evaluation Summary*:  
Under the single, strict energy closure target of $|\epsilon_{\text{energy}}| \le 2.0\%$ across all benchmark times:
- `EHCDarcy` passes all verification metrics across the entire 17-minute benchmark (energy closure -0.34% (-6.4 J) at 2 min, -0.23% (-8.4 J) at 6 min, and -0.13% (-10.3 J) at 17 min).
- Enthalpy-Porosity under tutorial settings ($\lambda = 0.2, n_{\text{outer}}=8$) fails energy closure (+7.35% (+132.3 J) *(FAIL)*, +7.37% (+260.7 J) *(FAIL)*, +7.49% (+590.4 J) *(FAIL)*) and late front deviation (5.43 mm (6.11%) *(FAIL)*).
- Increasing outer iterations to $n_{\text{outer}}=12$ under $\lambda = 0.2$ reduces error but still fails the $\le 2.0\%$ target across all times (+2.86% (+52.9 J) *(FAIL)*, +2.85% (+103.2 J) *(FAIL)*, +2.92% (+234.4 J) *(FAIL)*).
- Restoring code-default Newton relaxation $\lambda = 1.0$ at $n_{\text{outer}}=8$ allows Enthalpy-Porosity to achieve -0.06% (-1.2 J) at 2 min, +0.11% (+4.2 J) at 6 min, and -0.10% (-8.0 J) at 17 min, passing the single $\le 2.0\%$ target across all times and closely matching EHCDarcy.

---

## 2. Energy-Balance Audit, Post-Processing Verification & Formulation Sensitivity

### 2.1 Root Cause of the Initial Coarse Audit Discrepancy (+20.6% Apparent Surplus)
The initial coarse audit script (`energy_balance_audit.py`, archived in `energy_balance_results.txt`) evaluated the baseline Enthalpy-Porosity simulation using 60-second sampling intervals and reported apparent surpluses of $+31.3\%$ at 1 min, $+20.6\%$ at 2 min ($1601.5\text{ J}$ input vs. $1931.1\text{ J}$ stored), and $+6.6\%$ at 17 min.

Controlled step-resolved investigation ($\Delta t = 0.1\text{ s}$, 1,200 steps to 2 min) demonstrates that this apparent non-conservation was a **post-processing artifact** driven by two primary factors:
1. **Severe Temporal Undersampling of the Early Conduction Surge**:  
   At $t = 0.1\text{ s}$, transient conduction into the subcooled slab produces an initial heat flux peak of $\dot{Q}_{\text{hot}} \approx 72\text{--}81\text{ W}$ ($q''_{\text{wall}} \approx 2.3\text{--}2.6 \times 10^5\text{ W/m}^2$), decaying rapidly to $11.2\text{ W}$ by 60 s. Sampling exclusively at 60 s missed this surge, underestimating first-minute gross input by **$229.66\text{ J}$** ($1239.60\text{ J}$ step-resolved vs. $1009.94\text{ J}$ coarse). During the second minute ($60 \to 120\text{ s}$), diffusion decay is convex ($q \propto t^{-1/2}$), causing the coarse trapezoid to overestimate delivery by **$14.80\text{ J}$** ($591.58\text{ J}$ coarse vs. $576.78\text{ J}$ step-resolved). The net temporal undersampling error is:
   $$\Delta E_{\text{temporal}} = 229.66\text{ J} - 14.80\text{ J} = \mathbf{214.86\text{ J}} \quad (99.97\%\text{ of gross difference})$$
2. **Estimator Discretization Residual**:  
   At the 60 s and 120 s output timestamps, the finite-difference wall gradient ($Q_{\text{FD}}$) and OpenFOAM's native patch integral ($Q_{\text{native}}$) agree to within $< 0.01\%$ ($11.2215\text{ W}$ vs. $11.2209\text{ W}$ at 60 s; $8.4978\text{ W}$ vs. $8.4981\text{ W}$ at 120 s). The difference between the coarse native reconstruction ($1601.52\text{ J}$) and unrounded finite-difference calculation ($1601.45\text{ J}$) is only **$0.065\text{ J}$** ($0.03\%$ of the gross difference).

| Contribution to Gross Discrepancy ($t = 2\text{ min}$) | Energy (J) | Fraction of Difference |
| :--- | :---: | :---: |
| **Total Gross Difference** ($1816.38\text{ J} - 1601.45\text{ J}$) | **$214.93\text{ J}$** | **$100.0\%$** |
| **Temporal Undersampling Error** ($E_{\text{native, step}} - E_{\text{native, coarse}}$: $1816.38 - 1601.52$) | **$214.86\text{ J}$** | **$99.97\%$** |
| **Estimator Discretization Residual** ($E_{\text{native, coarse}} - E_{\text{FD, coarse}}$: $1601.52 - 1601.45$) | **$0.07\text{ J}$** | **$0.03\%$** |

*(The complete step-by-step arithmetic derivation reconciling the archived log display value of $1601.5\text{ J}$ is documented in [Appendix A](#appendix-a-historical-coarse-audit-reconstruction)).*

---

### 2.2 Resolution of Historical Boundary Accounting Artifacts (Sign Bug & Quadrature Shift)
In an earlier version of the time-resolved audit script (`exact_energy_balance_audit.py`), `EHCDarcy` was reported with negative closure discrepancies of $-1.26\%$ at 2 min, $-2.87\%$ at 6 min, $-4.39\%$ at 12 min, and $-5.86\%$ at 17 min. 

Inspection of the script source revealed a sign convention error in net boundary energy evaluation:
```python
# Historical bug in exact_energy_balance_audit.py (now resolved):
enet = ein - eout  # eout was already negative in OpenFOAM wallHeatFlux.dat!
```
In OpenFOAM's `wallHeatFlux` function object, heat entering the fluid domain is signed positive, while heat leaving across the isothermal cold boundary ($T_{\text{fluid}} > T_{\text{cold}} = 301.45\text{ K}$) produces a negative patch integral ($\dot{Q}_{\text{cold, raw}} < 0$). Subtracting a negative cumulative energy (`ein - (-eout) = ein + eout`) double-counted outward cold-wall heat loss as an artificial heat injection into the domain. At 17 min, cumulative cold-wall loss is 246.64 J; adding rather than subtracting this loss inflated apparent net input to 8646.65 J, producing the spurious -503.59 J (-5.82%) deficit against domain stored energy (8143.06 J).

When the physically correct definition is enforced ($E_{\text{net, in}} = E_{\text{hot}} - E_{\text{cold, out}} = E_{\text{hot}} + E_{\text{cold, raw}}$):
- At $t = 2\text{ min}$: $E_{\text{net, in}} = 1891.80 - 7.00 = 1884.80\text{ J}$, $E_{\text{stored}} = 1878.37\text{ J} \implies \Delta E = \mathbf{-6.43\text{ J}}$ ($\mathbf{-0.34\%}$).
- At $t = 6\text{ min}$: $E_{\text{net, in}} = 3733.53 - 48.25 = 3685.28\text{ J}$, $E_{\text{stored}} = 3676.92\text{ J} \implies \Delta E = \mathbf{-8.36\text{ J}}$ ($\mathbf{-0.23\%}$).
- At $t = 17\text{ min}$: $E_{\text{net, in}} = 8400.01 - 246.64 = 8153.37\text{ J}$, $E_{\text{stored}} = 8143.06\text{ J} \implies \Delta E = \mathbf{-10.31\text{ J}}$ ($\mathbf{-0.13\%}$).

To ensure this sign error cannot reoccur, a formal unit test with synthetic boundary fluxes (hot $+10\text{ W}$, cold $-2\text{ W} \implies E_{\text{net}} = 8\text{ J/s}$) and linearly decaying fluxes is committed in [test_exact_energy_balance_audit.py](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/validation/gallium2D/analysis/energy_balance/test_exact_energy_balance_audit.py).

---

### 2.3 Mathematical Definition of the Step-Resolved Energy Audit & Discrete Backward-Euler Quadrature
In OpenFOAM's finite-volume formulation with implicit backward-Euler time discretization (`Euler ddtScheme`), the discrete energy balance from time level $t^{n-1}$ to $t^n$ satisfies:
$$\frac{E_{\text{stored}}^n - E_{\text{stored}}^{n-1}}{\Delta t_n} = \dot{Q}_{\text{net}}(t^n) \iff \Delta E_{\text{stored}}^n = \dot{Q}_{\text{net}}(t^n)\,\Delta t_n$$
Summing over all discrete time steps $n = 1, 2, \dots, N$ yields the **exact discrete right-endpoint Riemann sum**:
$$E_{\text{net, in}}(t_N) = \sum_{n=1}^N \dot{Q}_{\text{net}}(t_n)\,\Delta t_n$$
where $\dot{Q}_{\text{net}}(t_n) = \dot{Q}_{\text{hot}}(t_n) - \dot{Q}_{\text{cold, out}}(t_n)$.

> **Note on Quadrature Rule & Discrete Energy Balance**:  
> In earlier audit iterations, a composite trapezoidal rule was applied between sampled time steps with a rectangular first step:
> $$E_{\text{trap}}(t_N) = \dot{Q}(t_1)\,\Delta t_1 + \sum_{k=2}^N \frac{\dot{Q}(t_{k-1}) + \dot{Q}(t_k)}{2}\,\Delta t_k$$
> Algebraically, this trapezoidal formulation equals the right-endpoint sum plus a boundary offset:
> $$E_{\text{trap}}(t_N) = \sum_{n=1}^N \dot{Q}(t_n)\,\Delta t_n + \frac{\dot{Q}(t_1) - \dot{Q}(t_N)}{2}\,\Delta t$$
> Because the initial conduction surge produces an intense spike ($\dot{Q}(t_1) \approx 81.2\text{ W}$ at $\Delta t = 0.1\text{ s}$) that rapidly decays to quasi-steady values ($\dot{Q}(t_N) \approx 8.5\text{ W}$ at 2 min, $4.5\text{ W}$ at 17 min), the trapezoidal rule added an artificial, non-physical quadrature offset of $+\frac{81.2 - 8.5}{2} \times 0.1 \approx +3.63\text{ J}$ at 2 min (and $+3.83\text{ J}$ at 17 min) on both walls. Apparent discrepancies at the $0.1\text{--}0.5\%$ level were therefore partly a post-processing quadrature artifact rather than solver physics. Adopting the exact right-endpoint sum $\sum_{n=1}^N \dot{Q}(t_n)\,\Delta t_n$ eliminates this spurious offset and achieves exact discrete compatibility with the backward-Euler time discretization.

Stored enthalpy relative to initial conditions ($T_{\text{init}} = 301.45\text{ K}$, $f_L = 0$) evaluates to:
$$E_{\text{sens}}(t) = \int_\Omega \rho C_p (T - T_{\text{init}})\,dV, \quad E_{\text{lat}}(t) = \int_\Omega \rho L f_L\,dV, \quad E_{\text{stored}}(t) = E_{\text{sens}}(t) + E_{\text{lat}}(t)$$
Energy conservation closure error is defined as:
$$\Delta E(t) = E_{\text{stored}}(t) - E_{\text{net, in}}(t), \quad \epsilon_{\text{energy}}(t) = \frac{\Delta E(t)}{E_{\text{net, in}}(t)} \times 100\%$$

---

### 2.4 Time-Resolved Authoritative Energy Audit History

#### Complete Transient History: `EHCDarcy` Formulation ($42 \times 32$ Grid, $n_{\text{outer}}=8$, $\Delta t = 0.1\text{ s}$)
| Time | $E_{\text{hot}}$ (J) | $E_{\text{cold, out}}$ (J) | $E_{\text{net, in}}$ (J) | $E_{\text{sens}}$ (J) | $E_{\text{lat}}$ (J) | $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Closure Error (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1 min** | $1298.65$ | $0.98$ | $1297.67$ | $59.10$ | $1233.29$ | $1292.38$ | $-5.29\text{ J}$ | **-0.41\%** |
| **2 min** | $1891.80$ | $7.00$ | $1884.80$ | $81.37$ | $1797.00$ | $1878.37$ | $-6.43\text{ J}$ | **-0.34\%** |
| **4 min** | $2853.10$ | $26.01$ | $2827.09$ | $107.27$ | $2712.27$ | $2819.54$ | $-7.55\text{ J}$ | **-0.27\%** |
| **6 min** | $3733.53$ | $48.25$ | $3685.28$ | $129.06$ | $3547.86$ | $3676.92$ | $-8.36\text{ J}$ | **-0.23\%** |
| **8 min** | $4593.83$ | $73.21$ | $4520.62$ | $149.51$ | $4362.21$ | $4511.72$ | $-8.90\text{ J}$ | **-0.20\%** |
| **10 min** | $5447.16$ | $101.60$ | $5345.56$ | $169.47$ | $5166.69$ | $5336.15$ | $-9.41\text{ J}$ | **-0.18\%** |
| **12 min** | $6294.16$ | $134.44$ | $6159.72$ | $187.90$ | $5962.05$ | $6149.96$ | $-9.77\text{ J}$ | **-0.16\%** |
| **14 min** | $7138.12$ | $173.08$ | $6965.04$ | $206.68$ | $6748.32$ | $6955.00$ | $-10.04\text{ J}$ | **-0.14\%** |
| **16 min** | $7979.22$ | $219.52$ | $7759.70$ | $224.28$ | $7525.23$ | $7749.51$ | $-10.19\text{ J}$ | **-0.13\%** |
| **17 min** | $8400.01$ | $246.64$ | $8153.37$ | $232.71$ | $7910.36$ | $8143.06$ | $-10.31\text{ J}$ | **-0.13\%** |

#### Complete Transient History: `enthalpyPorosity` (Default: $\lambda = 1.0, n_{\text{outer}}=8$, $\Delta t = 0.1\text{ s}$)
| Time | $E_{\text{hot}}$ (J) | $E_{\text{cold, out}}$ (J) | $E_{\text{net, in}}$ (J) | $E_{\text{sens}}$ (J) | $E_{\text{lat}}$ (J) | $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Closure Error (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1 min** | $1301.35$ | $0.98$ | $1300.38$ | $59.35$ | $1236.91$ | $1296.26$ | $-4.12\text{ J}$ | **-0.32\%** |
| **2 min** | $1894.55$ | $7.01$ | $1887.53$ | $81.64$ | $1804.69$ | $1886.33$ | $-1.20\text{ J}$ | **-0.06\%** |
| **4 min** | $2856.51$ | $26.14$ | $2830.37$ | $107.66$ | $2725.59$ | $2833.25$ | $+2.88\text{ J}$ | **+0.10\%** |
| **6 min** | $3737.97$ | $48.50$ | $3689.47$ | $129.40$ | $3564.23$ | $3693.64$ | $+4.16\text{ J}$ | **+0.11\%** |
| **8 min** | $4599.25$ | $73.57$ | $4525.68$ | $149.92$ | $4379.05$ | $4528.96$ | $+3.28\text{ J}$ | **+0.07\%** |
| **10 min** | $5453.61$ | $102.06$ | $5351.55$ | $169.75$ | $5183.13$ | $5352.87$ | $+1.32\text{ J}$ | **+0.02\%** |
| **12 min** | $6302.15$ | $134.97$ | $6167.18$ | $188.36$ | $5977.75$ | $6166.11$ | $-1.06\text{ J}$ | **-0.02\%** |
| **14 min** | $7147.29$ | $173.67$ | $6973.62$ | $206.94$ | $6762.95$ | $6969.89$ | $-3.73\text{ J}$ | **-0.05\%** |
| **16 min** | $7989.87$ | $220.12$ | $7769.74$ | $224.44$ | $7538.76$ | $7763.21$ | $-6.54\text{ J}$ | **-0.08\%** |
| **17 min** | $8410.93$ | $247.21$ | $8163.72$ | $232.79$ | $7922.90$ | $8155.68$ | $-8.04\text{ J}$ | **-0.10\%** |

*Latent Enthalpy Reconciliation at $t = 17\text{ min}$*:  
For EHCDarcy, direct cell summation yields $E_{	ext{lat}} = \sum \rho V_{\text{cell}} L \alpha = 7910.36	ext{ J}$, matching $m L f_L = 13790.40 \times 0.5736131 = 7910.36\text{ J}$.  
For EP ($\lambda = 1.0, n_{\text{outer}}=8$), direct cell summation yields $E_{\text{lat}} = 7922.90\text{ J}$, matching $m L f_L = 13790.40 \times 0.5745226 = 7922.90\text{ J}$. Both formulations reconcile latent enthalpy with liquid fraction to $0.01\text{ J}$.

---

### 2.5 Systematic Study of Newton Relaxation ($\lambda$), Outer Correctors, and Computational Cost
In OpenFOAM's `enthalpyPorosityPhaseChangeModel.C`, the liquid fraction is updated after the sensible enthalpy matrix solve via Newton projection:
$$f_{\text{new}} = \min\left(\max\left(f + \lambda \frac{T - T_f}{L/C_p + \Delta T_{\text{window}}}, 0\right), 1\right)$$
where $\lambda$ is the Newton relaxation factor (code default: `1.0`, the exact Newton step for a linear melting window). In the baseline tutorial case dictionary, this parameter was set to `lambda 0.2`, damping each outer iteration correction to only $20\%$ of the full Newton step.

To decouple physical modeling capabilities from under-relaxation artifacts and rigorously test outer corrector convergence, a systematic sweep of $\lambda \in \{1.0, 0.5, 0.2\}$ across PIMPLE outer correctors $n_{\text{outer}} \in [2, 16]$ was executed alongside `EHCDarcy` at $n_{\text{outer}} \in [2, 8]$ to $t = 120\text{ s}$ ($2\text{ min}$). All execution times report OpenFOAM CPU solver time (`ExecutionTime` parsed from solver logs) alongside the exact number of 2D outer corrector passes ($1200 \times n_{\text{outer}}$):

| Model Formulation | $\lambda$ | $n_{\text{outer}}$ | Outer Passes | $f_L(2\text{ min})$ | Net Input $E_{\text{net}}$ (J) | Stored $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Closure Error (%) | CPU Time (s) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| EHCDarcy | — | 2 | 2400 | 0.1315 | 1733.95 J | 1895.54 J | +161.59 J | +9.32% | 10.98 s |
| EHCDarcy | — | 3 | 3600 | 0.1303 | 1885.40 J | 1878.26 J | -7.14 J | -0.38% | 15.57 s |
| EHCDarcy | — | 4 | 4800 | 0.1303 | 1885.15 J | 1878.33 J | -6.83 J | -0.36% | 20.01 s |
| EHCDarcy | — | 5 | 6000 | 0.1303 | 1884.94 J | 1878.36 J | -6.58 J | -0.35% | 23.83 s |
| **EHCDarcy** | — | **8** | 9600 | **0.1303** | **1884.80 J** | **1878.37 J** | **-6.43 J** | **-0.34%** | 34.95 s |
| enthalpyPorosity | 1.0 | 2 | 2400 | 0.1650 | 1511.41 J | 2370.39 J | +858.98 J | +56.83% | 10.73 s |
| enthalpyPorosity | 1.0 | 3 | 3600 | 0.1358 | 1758.45 J | 1956.42 J | +197.97 J | +11.26% | 15.19 s |
| enthalpyPorosity | 1.0 | 4 | 4800 | 0.1262 | 1871.94 J | 1819.56 J | -52.38 J | -2.80% | 19.50 s |
| enthalpyPorosity | 1.0 | 5 | 6000 | 0.1235 | 1920.57 J | 1779.98 J | -140.59 J | -7.32% | 25.11 s |
| enthalpyPorosity | 1.0 | 6 | 7200 | 0.1245 | 1928.53 J | 1795.67 J | -132.86 J | -6.89% | 30.82 s |
| enthalpyPorosity | 1.0 | 7 | 8400 | 0.1278 | 1910.76 J | 1842.37 J | -68.39 J | -3.58% | 35.12 s |
| **enthalpyPorosity** | 1.0 | **8** | 9600 | **0.1309** | **1887.53 J** | **1886.33 J** | **-1.20 J** | **-0.06%** | 39.38 s |
| enthalpyPorosity | 1.0 | 10 | 12000 | 0.1320 | 1869.38 J | 1902.18 J | +32.80 J | +1.75% | 47.16 s |
| enthalpyPorosity | 1.0 | 12 | 14400 | 0.1307 | 1874.92 J | 1883.31 J | +8.40 J | +0.45% | 56.63 s |
| enthalpyPorosity | 1.0 | 16 | 19200 | 0.1305 | 1878.15 J | 1880.61 J | +2.45 J | +0.13% | 73.08 s |
| enthalpyPorosity | 0.5 | 2 | 2400 | 0.1750 | 1424.08 J | 2513.52 J | +1089.45 J | +76.50% | 10.64 s |
| enthalpyPorosity | 0.5 | 3 | 3600 | 0.1480 | 1634.82 J | 2130.52 J | +495.71 J | +30.32% | 15.39 s |
| enthalpyPorosity | 0.5 | 5 | 6000 | 0.1334 | 1803.06 J | 1922.33 J | +119.27 J | +6.62% | 24.24 s |
| enthalpyPorosity | 0.5 | 8 | 9600 | 0.1294 | 1874.30 J | 1865.11 J | -9.19 J | -0.49% | 37.94 s |
| enthalpyPorosity | 0.2 | 2 | 2400 | 0.1762 | 1344.95 J | 2533.71 J | +1188.76 J | +88.39% | 10.93 s |
| enthalpyPorosity | 0.2 | 3 | 3600 | 0.1527 | 1547.23 J | 2198.28 J | +651.06 J | +42.08% | 15.64 s |
| enthalpyPorosity | 0.2 | 5 | 6000 | 0.1395 | 1707.54 J | 2009.79 J | +302.25 J | +17.70% | 24.40 s |
| enthalpyPorosity | 0.2 | 8 | 9600 | 0.1340 | 1798.86 J | 1931.11 J | +132.25 J | +7.35% | 38.65 s |
| enthalpyPorosity | 0.2 | 12 | 14400 | 0.1317 | 1845.83 J | 1898.69 J | +52.86 J | +2.86% | 55.65 s |

**Key Findings on Enthalpy-Porosity vs. EHCDarcy Convergence Dynamics**:
1. **Newton Step Relaxation Dynamics ($\lambda = 1.0$)**:  
   Under code-default $\lambda = 1.0$, the discrete energy closure error exhibits a damped, oscillatory convergence path:
   $$\epsilon_{\text{energy}}(2\text{ min}): \quad +56.8% (n=2), +11.3% (n=3), -2.8% (n=4), -7.3% (n=5), -6.9% (n=6), -3.6% (n=7), -0.1% (n=8), +1.8% (n=10), +0.4% (n=12), +0.1% (n=16)$$
   At low corrector counts ($n \le 3$), enthalpy-porosity exhibits large positive surpluses (+56.8% at $n=2$, +11.3% at $n=3$) because the explicit latent heat update lags the rapid early conduction surge. As $n$ increases, the scheme overcorrects into a negative deficit (-7.3% at $n=5$), crosses zero near $n=7\text{--}8$ (-0.06% at $n=8$), and asymptotically stabilizes below $0.5\%$ for $n \ge 12$ (+0.45% at $n=12$, +0.13% at $n=16$).
2. **EHCDarcy Monotonicity and Rapid Convergence (from $n_{\text{outer}} \ge 3$)**:  
   In contrast, `EHCDarcy` incorporates latent heat directly into the left-hand matrix capacity via a path-integral effective heat capacity $C_p^{\text{eff}}(T)$. As a result:
   $$\epsilon_{\text{energy}}(2\text{ min}): \quad +9.32% (n=2), -0.38% (n=3), -0.36% (n=4), -0.35% (n=5), -0.34% (n=8)$$
   While EHCDarcy exhibits an initial sign transition from +9.32% at $n=2$ (10.98 s CPU) to -0.38% at $n_{\text{outer}} = 3$ (15.57 s CPU), it is **strictly monotonic and bounded within $[-0.38\%, -0.34\%]$ for all $n_{\text{outer}} \ge 3$** (-0.36% at $n=4$, -0.35% at $n=5$, -0.34% at $n=8$), completely avoiding the large under-damped sign oscillations seen in EP.
3. **Computational Efficiency & Outer Corrector Cost**:  
   Per outer iteration pass, both formulations incur virtually identical CPU execution times ($\approx 4.2\text{--}4.3\text{ ms}$ per outer-corrector pass of the whole mesh, or $\approx 3.2\,\mu\text{s}$ per cell per pass; e.g. at $n_{\text{outer}}=3$, EHCDarcy requires 15.57 s for 3,600 outer passes vs. 15.19 s for EP).  
   However, `EHCDarcy` achieves verified discrete energy closure ($\le 0.4\%$) with only $n_{\text{outer}} = 3$ (closure error -0.38%, CPU time 15.57 s). Enthalpy-Porosity at $\lambda = 1.0$ requires $n_{\text{outer}} \ge 12$ (CPU time 56.63 s, 14,400 passes) to reliably settle below $\le 0.5\%$ closure (+0.45%). Consequently, EHCDarcy achieves equivalent energy conservation in **$25\%$ of the outer corrector passes and $3.6\times$ less CPU time (15.57 s vs. 56.63 s)**, confirming that EHCDarcy is substantially more economical for fluid-region melting.

---

## 3. Temporal Discretization Sensitivity & Convergence Study (EHCDarcy Formulation)

### 3.1 Brent et al. (1988) Time-Step Schedule Sensitivity
Brent, Voller, and Reid (1988, Section 4.1, p. 308) documented:
- First 200 time steps at $\Delta t = 0.05\text{ s}$ ($t = 0 \to 10\text{ s}$).
- Subsequent time steps at $\Delta t = 0.5\text{ s}$ ($t > 10\text{ s}$).

Executing this exact schedule in `phaseChangeMultiRegionFoam` yields:
- At $t = 60\text{ s}$ ($1\text{ min}$): $f_L = 0.0824$.
- Brent et al. (1988, Figure 5) reference at 60 s: $f_L \approx 0.11$.
The exact historical schedule produces an early liquid fraction approximately $25\%$ lower than Brent's reported curve, demonstrating that OpenFOAM's discrete PIMPLE coupling does not identically reproduce Brent's original 1988 numerical output when using his coarse late time-step.

### 3.2 Formal Temporal Convergence Study ($\Delta t = 0.1\text{ s}, 0.05\text{ s}, 0.025\text{ s}$)
To establish rigorous temporal convergence, uniform time steps were executed to $t = 120\text{ s}$:

| Uniform Time Step | Total Steps (120 s) | $f_L(60\text{ s})$ | $f_L(120\text{ s})$ | Mid-Height Front $x_{\text{front}}(120\text{ s})$ |
| :--- | :---: | :---: | :---: | :---: |
| $\Delta t = 0.100\text{ s}$ (Baseline) | 1,200 | $0.0886$ | $0.1303$ | $13.62\text{ mm}$ |
| $\Delta t = 0.050\text{ s}$ | 2,400 | $0.0887$ | $0.1305$ | $13.64\text{ mm}$ |
| $\Delta t = 0.025\text{ s}$ | 4,800 | $0.0888$ | $0.1306$ | $13.65\text{ mm}$ |

Refining $\Delta t$ by a factor of 4 produces a change in $f_L(120\text{ s})$ of only **0.0003 (0.23%)** and shifts the mid-height melt front by only **0.03 mm (0.22%)**. The baseline time step $\Delta t = 0.1\text{ s}$ provides time-step independent results.

### 3.3 Investigation of the Startup Conservation Deficit under $\Delta t$ Halving
To investigate the physical nature of the small, persistent deficit in EHCDarcy ($\approx -0.34\%$ at 2 min, $-6.43\text{ J}$), a targeted test halved $\Delta t$ to $0.05\text{ s}$ (2,400 steps):
- **Baseline $\Delta t = 0.10\text{ s}$ (1,200 steps)**: $E_{\text{net}} = 1883.84\text{ J}$, $E_{\text{stored}} = 1877.41\text{ J}$, deficit $\Delta E = -6.43\text{ J}$ (closure $-0.341\%$).
- **Halved $\Delta t = 0.05\text{ s}$ (2,400 steps)**: $E_{\text{net}} = 1880.92\text{ J}$, $E_{\text{stored}} = 1879.42\text{ J}$, deficit $\Delta E = -1.50\text{ J}$ (closure $-0.080\%$).

Halving $\Delta t$ reduces the 2-minute conservation deficit by **76.7%** (a factor of $> 4$, demonstrating that the deficit falls faster than first order during the startup transient). Furthermore, between $t = 60\text{ s}$ and $t = 120\text{ s}$, the deficit grows by only $-0.21\text{ J}$ ($-1.29\text{ J} \to -1.50\text{ J}$). This confirms that the deficit is purely a temporal discretization / startup transient effect, and NOT an inherent spatial or discrete formulation inconsistency in EHCDarcy.

---

## 4. Latent Advection Term Analysis: `div(phi, phaseFraction)`

### 4.1 Physical Sharp Interface vs. Numerical Mushy Zone
In physical pure metal melting, phase change occurs at a discrete sharp interface $\Gamma(t)$. Liquid velocity satisfies the no-slip condition at the interface ($\\mathbf{u}_{\text{interface}} = 0$), so physical advection of latent heat is identically zero.

In fixed-grid diffuse-interface methods, the phase boundary is numerically broadened across a finite transition window $\Delta T_{\text{window}}$. Across this numerical mushy zone:
- A cell-centred scale diagnostic $\sum V_i \rho L (\\mathbf{u}_i \cdot \nabla f_{L,i})$ across all 85 mushy cells at $t = 1020\text{ s}$ demonstrates that the median cell-wise ratio in the core mushy zone ($0.01 < f_L < 0.90$, 80 cells) is small:
  $$\operatorname{median}_i \left(\frac{|\\mathbf{u}_i \cdot \nabla f_{L,i}|}{|\partial f_{L,i}/\partial t|}\right) = 0.99\% \approx 1.0\%$$
- In the outer 5 fringe cells ($f_L \ge 0.90$), localized advection reaches $9.20\text{ W}$ absolute and $-7.52\text{ W}$ signed. Under continuum incompressibility, these signs represent local divergence in the core and convergence into the fringe, canceling across the closed cavity. Omitting numerical latent advection is an established and physically motivated benchmark approximation standard across classical CFD phase-change literature.

### 4.2 Implementation Verification in `phaseChangeMultiRegionFoam`
Inspection of `EEqn.H` confirms:
- The energy equation solves for sensible enthalpy $h$ with standard advection $\nabla \cdot (\rho \\mathbf{u} h)$.
- Latent heat is governed solely by the unsteady source term $L \frac{\partial(\rho f_L)}{\partial t}$ in `enthalpyPorosity` (or via $C_p^{\text{eff}}$ in `EHCDarcy`).
- The entry `div(phi,phaseFraction)` previously present in `system/fluid/fvSchemes` was an unused legacy dictionary token that was never invoked during matrix assembly.

---

## 5. Spatial Grid Sensitivity & ASME / Roache GCI Analysis (EHCDarcy Formulation, to 17 min)

Three systematically refined meshes were executed to the full $17\text{ min}$ ($1020\text{ s}$) duration:
- **Grid 3 (Coarse)**: $32 \times 24$ (768 cells, $h_3 = 2.71\text{ mm}$)
- **Grid 2 (Medium)**: $42 \times 32$ (1344 cells, $h_2 = 2.05\text{ mm}$)
- **Grid 1 (Fine)**: $63 \times 48$ (3024 cells, $h_1 = 1.37\text{ mm}$)

Refinement ratios: $r_{21} = h_2 / h_1 = 1.500$, $r_{32} = h_3 / h_2 = 1.3223$.

### 5.1 Grid Convergence & Sensitivity Table

| Grid | Resolution | Total Cells | $f_L(10\text{ min})$ | $f_L(17\text{ min})$ | $x_{\text{front}}(17\text{m}, y/H=0.2)$ | $x_{\text{front}}(17\text{m}, y/H=0.5)$ | $x_{\text{front}}(17\text{m}, y/H=0.8)$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Coarse (Grid 3) | $32 \times 24$ | 768 | $0.3698$ | $0.5682$ | $29.18\text{ mm}$ | $56.01\text{ mm}$ | $66.75\text{ mm}$ |
| Medium (Grid 2) | $42 \times 32$ | 1344 | $0.3747$ | $0.5736$ | $28.19\text{ mm}$ | $57.06\text{ mm}$ | $67.80\text{ mm}$ |
| Fine (Grid 1) | $63 \times 48$ | 3024 | $0.3760$ | $0.5725$ | $27.58\text{ mm}$ | $57.68\text{ mm}$ | $66.87\text{ mm}$ |
| **Richardson Extrapolated ($\phi_{\text{ext}}$)** | — | $\infty$ | — | — | — | — | **$57.99\text{ mm}$** | — |

### 5.2 Formal ASME / Roache Grid Convergence Analysis
1. **Global Liquid Fraction $f_L(17\text{ min})$**:
   The computed sequence from coarse to fine is $0.5682 \to 0.5736 \to 0.5725$. Because the medium grid is slightly higher than the fine grid (by 0.0011, or 0.19%), this sequence is non-monotonic. Consequently, a formal monotonic Richardson extrapolation cannot be claimed for $f_L(17\text{ min})$. The total spread across all three grids is only 0.0054 (0.94%), establishing that the global melt fraction is spatially insensitive to grid refinement beyond $42 \times 32$.
2. **Phase Front Position at Mid-Height $y/H = 0.5$**:
   The front coordinate sequence is strictly monotonic ($56.01\text{ mm} \to 57.06\text{ mm} \to 57.68\text{ mm}$). For unequal refinement ratios ($r_{21} = 1.500, r_{32} = 1.3223$), solving the generalized transcendental ASME equation yields:
   - Observed order of convergence: $p = 2.70$
   - Richardson extrapolated front position: $\phi_{\text{ext}} = 57.99\text{ mm}$
   - Fine-grid relative error: $e_a^{21} = 1.07\%$
   - **Fine-Grid Convergence Index ($\text{GCI}_{\text{fine}}^{21}$)**: **0.67%**
   - **Medium-Grid Convergence Index ($\text{GCI}_{\text{med}}^{32}$)**: **2.04%**
   - **Asymptotic Range Consistency Ratio**:
     $$\frac{\mathrm{GCI}_{\mathrm{med}}^{32}}{r_{21}^p \, \mathrm{GCI}_{\text{fine}}^{21}} = \frac{2.04\%}{1.500^{2.70} \times 0.67\%} \approx \mathbf{1.01}$$
     The consistency ratio is within $1\%$ of theoretical asymptotic limit $1.00$.

**Scope of Grid Study**: Formal ASME GCI verification is established specifically for `EHCDarcy`; no formal grid-convergence study was conducted for `enthalpyPorosity`.

---

## 6. Late-Time Parameter Sensitivity Analysis (to 17 min)

### 6.1 Darcy Mushy Penalty Constant Sweep ($C_u \in [10^4, 10^8]$)
Sweeping $C_u$ across four orders of magnitude to $t = 17\text{ min}$ on the baseline mesh:

| Darcy Constant $C_u$ [kg/(m$^3$ s)] | $f_L(10\text{ min})$ | $f_L(17\text{ min})$ | $\overline{\text{Nu}}_{\text{hot}}(17\text{ min})$ | $x_{\text{front}}(y/H=0.2)$ | $x_{\text{front}}(y/H=0.8)$ | Front Tilt Ratio ($x_{0.8}/x_{0.2}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| $1.0 \times 10^4$ (Weak damping) | $0.3906$ | $0.5987$ | $4.74$ | $26.34\text{ mm}$ | $72.26\text{ mm}$ | $2.74$ |
| $1.0 \times 10^5$ | $0.3817$ | $0.5848$ | $4.62$ | $27.39\text{ mm}$ | $69.93\text{ mm}$ | $2.55$ |
| **$1.6 \times 10^6$ (Literature Baseline)** | **$0.3747$** | **$0.5736$** | **$4.50$** | **$28.19\text{ mm}$** | **$67.80\text{ mm}$** | **$2.41$** |
| $1.0 \times 10^7$ | $0.3720$ | $0.5694$ | $4.47$ | $28.43\text{ mm}$ | $66.77\text{ mm}$ | $2.35$ |
| $1.0 \times 10^8$ (Severe damping) | $0.3703$ | $0.5667$ | $4.44$ | $28.72\text{ mm}$ | $66.14\text{ mm}$ | $2.30$ |

### 6.2 Melting Window Width Sweep ($\Delta T_{\text{melt}} \in [0.1, 0.5]\text{ K}$)

| Melting Window $\Delta T_{\text{melt}}$ | $f_L(10\text{ min})$ | $f_L(17\text{ min})$ | $\overline{\text{Nu}}_{\text{hot}}(17\text{ min})$ | $x_{\text{front}}(y/H=0.2)$ | $x_{\text{front}}(y/H=0.8)$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| $0.1\text{ K}$ ($302.93 \pm 0.05\text{ K}$) | $0.3738$ | $0.5727$ | $4.57$ | $27.81\text{ mm}$ | $67.82\text{ mm}$ |
| $0.2\text{ K}$ ($302.93 \pm 0.10\text{ K}$) | $0.3736$ | $0.5728$ | $4.56$ | $27.90\text{ mm}$ | $67.81\text{ mm}$ |
| **$0.5\text{ K}$ ($302.93 \pm 0.25\text{ K}$, Baseline)** | **$0.3747$** | **$0.5736$** | **$4.50$** | **$28.19\text{ mm}$** | **$67.80\text{ mm}$** |

Varying the numerical phase transition window from $0.5\text{ K}$ down to $0.1\text{ K}$ changes $f_L(17\text{ min})$ by only **0.0009 (0.16%)**.

---

## 7. Reference Data Limits & Literature Benchmark Comparison

### 7.1 Distinguishing Numerical Benchmark Fronts from Experimental Measurements
All quantitative front morphology metrics in Section 1.2 are calculated directly against the digitized numerical simulation fronts published in Figure 9 of Brent et al. (1988) (digitization uncertainty envelope $\pm 1.0\text{ mm}$, or $\pm 1.1\%$ of cavity width $W$).

### 7.2 Experimental Comparison with Gau & Viskanta (1986)
Gau and Viskanta (1986) performed the foundational physical experiments. While qualitative agreement in front morphology is evident, exact quantitative matching with 2D idealized CFD is constrained by physical experimental realities:
1. **Wall Thermal Inertia**: Experimental hot walls could not experience an instantaneous step change to $311.15\text{ K}$ due to heating jacket thermal inertia and early heat loss.
2. **Pour-Out Measurement Disruption**: Measuring front positions by pouring out molten gallium at selected intervals disrupted the developed natural convection boundary layer.
3. **3D Endwall Drag & Crystal Anisotropy**: Physical gallium crystals exhibit anisotropic thermal conductivity ($k \in [28, 41]\text{ W/(m K)}$), and side glass windows introduce 3D viscous drag.

### 7.3 Consistency Check with Ye & Arıcı (2023) Numerical Correlation
Ye & Arıcı (2023, Section 3.2.2, Eq. 6) correlated their grid-converged numerical liquid fractions as:
$$\phi_n = 2.708\,\tau^{0.712}\,\mathrm{Ra}^{0.0243}\,\mathrm{Ar}^{-0.0515}$$
where $\tau = \mathrm{Fo} \cdot \mathrm{Ste} = 1.7761 \times 0.03912 = 0.06948$, $\mathrm{Ra} = 6.063 \times 10^5$, and $\mathrm{Ar} = 0.7143$.

Evaluating Eq. (6) via the committed verification script ([evaluate_ye_arici_correlation.py](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/validation/gallium2D/analysis/evaluate_ye_arici_correlation.py)) at $t = 17\text{ min}$ ($1020\text{ s}$) yields:
$$\phi_n^{\text{Ye \& Arıcı (2023)}} = \mathbf{0.5703}$$
- `phaseChangeFoam` (EHCDarcy, $42 \times 32$): **$f_L = 0.5736$** (Discrepancy: **$+0.0033$, or +0.33 percentage points**).
- `phaseChangeFoam` (EP, $\lambda = 1.0, n=8$): **$f_L = 0.5745$** (Discrepancy: **$+0.0042$, or +0.42 percentage points**).

**Notes on Scaling Dominance**:  
Because the exponent on $\mathrm{Ra}$ is only $0.0243$, varying $\mathrm{Ra}$ by a factor of 10 alters the predicted liquid fraction by only $10^{0.0243} = +5.75\%$. Agreement with this correlation primarily verifies transient conduction and Stefan scaling ($\tau^{0.712}$), rather than fine convective boundary-layer dynamics.

---

## 8. Summary Conclusions & Verification Scope Assessment

### 8.1 Summary of Methodological Resolutions
1. **Energy-Balance Audit and Boundary Accounting Resolutions**:
   - **Unified $\le 2.0\%$ Acceptance Target**: A single strict criterion of $|\epsilon_{\text{energy}}| \le 2.0\%$ across all times ($t \le 17\text{ min}$) is enforced. `EHCDarcy` achieves complete energy closure across the entire transient (-0.13% to -0.34%, well within tolerance).
   - **Historical Audit Sign Bug**: In earlier versions of `exact_energy_balance_audit.py`, subtracting negative raw cold-wall heat flux double-counted cold-wall heat loss as an artificial heat source, creating false deficits ($-1.26\%$ to $-5.86\%$). Correcting to $E_{\text{net}} = E_{\text{hot}} - E_{\text{cold, out}}$ resolved true physical closure.
   - **Discrete Right-Endpoint Quadrature**: OpenFOAM's implicit backward-Euler operator evaluates flux discretely at the new time step $t^n$. Adopting the right-endpoint sum $\sum \dot{Q}(t_n)\,\Delta t_n$ eliminates the spurious $+3.63\text{ J}$ trapezoidal offset caused by the early conduction surge peak, bringing discrete energy conservation to within $< 0.35\%$.
   - **Enthalpy-Porosity Newton Relaxation ($\lambda$)**: Tutorial under-relaxation ($\lambda = 0.2$) caused severe non-conservation (+7.35% at 2 min, +7.49% at 17 min). Setting code-default $\lambda = 1.0$ at $n_{\text{outer}}=8$ achieves -0.06% at 2 min and -0.10% at 17 min. However, EP exhibits under-damped Newton oscillations and fails at low outer correctors ($n_{\text{outer}} \le 5$), whereas `EHCDarcy` conserves energy down to $n_{\text{outer}}=3$ (-0.38%, CPU time 15.57 s).
   - **Coarse Integration Artifact**: The coarse $1601.5\text{ J}$ script underestimated first-minute delivery by $229.66\text{ J}$ due to temporal undersampling of the early conduction surge. Temporal undersampling accounts for $99.97\%$ of the coarse gross discrepancy.

2. **Spatial and Temporal Sensitivity Bounds**:
   - Mid-height front position at 17 min demonstrates strict monotonic spatial convergence with an ASME observed order of $p = 2.70$, an asymptotic consistency ratio of 1.01, and a fine-grid $\text{GCI}$ of 0.67%. Formal GCI grid-convergence is established for `EHCDarcy` only.
   - Global melt fraction matches Ye & Arıcı (2023, Eq. 6) within 0.33 percentage points (0.5736 vs. 0.5703).

### 8.2 Scope Distinction: Numerical Benchmark Verification vs. Physical Experimental Validation
- **Numerical Benchmark Verification**: This report provides rigorous evidence of **numerical benchmark consistency and code-to-code verification** against established 2D CFD literature (Brent et al. 1988, Hannoun et al. 2003, Ye & Arıcı 2023). Under the stated idealized assumptions, `phaseChangeMultiRegionFoam` achieves a mean phase-front deviation of 3.46 mm (3.89% of cavity width) and maximum local deviation of 6.34 mm (7.13% of cavity width) relative to Brent et al. (1988), with an energy closure discrepancy of -0.13%.
- **Distinction from Experimental Validation**: These results **do not constitute a full quantitative validation against independent physical experimental measurements**. Quantitative alignment with physical experiments (such as Gau & Viskanta 1986) requires resolving 3D endwall boundary-layer drag, conjugate thermal inertia of the heating apparatus, solid-state gallium crystallographic anisotropy ($k \in [28, 41]\text{ W/(m K)}$), and non-intrusive front measurement methods.
- **Engineering Tolerances**: Acceptance thresholds (front mean deviation $\le 5\% W$, liquid fraction within $\pm 0.05$) are classified as engineering tolerances suitable for computational thermal engineering.

---

## Appendix A: Historical Coarse-Audit Reconstruction

To preserve complete arithmetic traceability for how the initial post-processing script (`energy_balance_audit.py`, archived in `energy_balance_results.txt`) arrived at its coarse estimate of $1601.5\text{ J}$, the exact algorithm and rounding transformations are documented below:

1. **First Interval ($0 \to 60\text{ s}$)**:  
   $$E_{0 \to 60} = \dot{Q}(60\text{ s}) \times 60\text{ s} \times 1.5$$
   Evaluating with the unrounded finite-difference flux sample ($Q_{\text{FD}} = 11.22090\text{ W}$) yields:
   $$E_{\text{FD}, 0 \to 60} = 11.22090 \times 60 \times 1.5 = 1009.881\text{ J} \quad (\text{displayed as } 1009.9\text{ J})$$
   Evaluating with the unrounded native patch flux ($Q_{\text{native}} = 11.22152\text{ W}$) yields:
   $$E_{\text{native}, 0 \to 60} = 11.22152 \times 60 \times 1.5 = 1009.937\text{ J} \quad (\text{displayed as } 1009.94\text{ J})$$

2. **Second Interval ($60 \to 120\text{ s}$)**:  
   $$E_{\text{FD}, 60 \to 120} = \frac{11.22090 + 8.49810}{2} \times 60 = 591.570\text{ J}$$
   $$E_{\text{native}, 60 \to 120} = \frac{11.22152 + 8.49778}{2} \times 60 = 591.579\text{ J} \approx 591.58\text{ J}$$

3. **Cumulative Two-Minute Coarse Reconstruction**:  
   - *Unrounded Finite-Difference Total*:
     $$E_{\text{FD, coarse}} = 1009.881\text{ J} + 591.570\text{ J} = \mathbf{1601.451\text{ J}}$$
     In `energy_balance_results.txt`, formatting under `{cum_E_in:<11.1f}` printed this total as **`1601.5 J`** (rounding upward by $+0.0490\text{ J}$).
   - *Unrounded Native Flux Reconstruction*:
     $$E_{\text{native, coarse}} = 1009.937\text{ J} + 591.579\text{ J} = \mathbf{1601.5158\text{ J} \approx 1601.52\text{ J}}$$
   - *Estimator Difference*:
     $$\Delta E_{\text{estimator}} = 1601.5158\text{ J} - 1601.4510\text{ J} = \mathbf{0.0648\text{ J}}$$

4. **Reconciliation with Displayed Log Value ($1601.5\text{ J}$)**:  
   $$\begin{aligned}
   \Delta E_{\text{archived display}} &= \Delta E_{\text{native temporal}} + \Delta E_{\text{estimator}} - \Delta E_{\text{archive rounding}} \\
   &= 214.8622\text{ J} + 0.0648\text{ J} - 0.0490\text{ J} \\
   &= \mathbf{214.8780\text{ J}} \approx \mathbf{214.88\text{ J}}
   \end{aligned}$$

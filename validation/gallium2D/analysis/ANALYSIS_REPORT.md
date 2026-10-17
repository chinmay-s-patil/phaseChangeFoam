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
| **Energy Closure ($t = 2\text{ min}$)** | $\le 2.0\%$ | $+6.72\%$ ($+121.6\text{ J}$) *(FAIL)* | $+2.30\%$ ($+42.6\text{ J}$) *(FAIL)* | **$-0.26\%$ ($-4.9\text{ J}$)** | **$-0.10\%$ ($-1.94\text{ J}$)** | **PASSED** | **PASSED** |
| **Energy Closure ($t = 6\text{ min}$)** | $\le 2.0\%$ | $+3.54\%$ ($+128.2\text{ J}$) *(FAIL)* | $+1.75\%$ ($+64.1\text{ J}$) | **$+0.01\%$ ($+0.4\text{ J}$)** | **$-0.11\%$ ($-4.0\text{ J}$)** | **PASSED** | **PASSED** |
| **Energy Closure ($t = 17\text{ min}$)** | $\le 2.0\%$ | $+5.18\%$ ($+416.7\text{ J}$) *(FAIL)* | $+2.97\%$ ($+238.5\text{ J}$) *(FAIL)* | **$-0.14\%$ ($-11.8\text{ J}$)** | **$-0.073\%$ ($-5.92\text{ J}$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 2\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $1.95\text{ mm}$ ($2.19\%$) | $1.72\text{ mm}$ ($1.93\%$) | $1.70\text{ mm}$ ($1.91\%$) | **$1.67\text{ mm}$ ($1.88\%$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 6\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $2.60\text{ mm}$ ($2.93\%$) | $2.12\text{ mm}$ ($2.38\%$) | $1.89\text{ mm}$ ($2.13\%$) | **$1.87\text{ mm}$ ($2.11\%$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 10\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $3.51\text{ mm}$ ($3.95\%$) | $2.84\text{ mm}$ ($3.19\%$) | $2.38\text{ mm}$ ($2.68\%$) | **$2.35\text{ mm}$ ($2.65\%$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 17\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $5.43\text{ mm}$ ($6.11\%$) *(FAIL)* | $4.18\text{ mm}$ ($4.70\%$) | $3.52\text{ mm}$ ($3.96\%$) | **$3.46\text{ mm}$ ($3.89\%$)** | **PASSED** | **PASSED** |
| **Max Front Dev. ($t = 17\text{ min}$)** | $\le 10.0\%$ $W$ ($8.89\text{ mm}$) | $9.74\text{ mm}$ ($10.95\%$) *(FAIL)* | $7.85\text{ mm}$ ($8.83\%$) | $6.48\text{ mm}$ ($7.29\%$) | **$6.34\text{ mm}$ ($7.13\%$)** | **PASSED** | **PASSED** |
| **$f_L$ at $t = 17\text{ min}$** (Brent: $0.5334$) | $\pm 0.05$ | $0.5967$ ($+0.0633$) *(FAIL)* | $0.5822$ ($+0.0488$) | **$0.5745$ ($+0.0411$)** | **$0.5736$ ($+0.0402$)** | **PASSED** | **PASSED** |
| **Ye & Arıcı (2023) Eq. 6** ($0.5703$) | $\pm 0.03$ | $+0.0264$ | $+0.0119$ | **$+0.0042$ ($0.42$ pts)** | **$+0.0033$ ($0.33$ pts)** | **PASSED** | **PASSED** |

*Scorecard Evaluation Summary*:  
Under the unified energy closure target of $|\epsilon_{\text{energy}}| \le 2.0\%$, `EHCDarcy` passes all verification metrics across the entire 17-minute benchmark (energy closure $-0.10\%$ at 2 min, $-0.11\%$ at 6 min, and $-0.073\%$ at 17 min). Enthalpy-Porosity with tutorial settings ($\lambda = 0.2, n_{\text{outer}}=8$) fails energy closure ($+6.72\%$, $+3.54\%$, $+5.18\%$) and late front deviation. While increasing outer iterations to $n_{\text{outer}}=12$ under $\lambda = 0.2$ reduces error, it still fails the $\le 2.0\%$ target at 2 min ($+2.30\%$) and 17 min ($+2.97\%$). However, when code-default Newton relaxation $\lambda = 1.0$ is restored, Enthalpy-Porosity achieves $-0.26\%$ at 2 min, $+0.01\%$ at 6 min, and $-0.14\%$ at 17 min, passing the single $\le 2.0\%$ target across all times and closely matching EHCDarcy.

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

### 2.2 Resolution of Historical Post-Processing Artifacts (EHCDarcy -1.26%, -2.87%, -4.39%)
In an earlier version of the time-resolved audit script (`exact_energy_balance_audit.py`), `EHCDarcy` was reported with negative closure discrepancies of $-1.26\%$ at 2 min, $-2.87\%$ at 6 min, $-4.39\%$ at 12 min, and $-5.86\%$ at 17 min. 

Inspection of the script source revealed a sign convention error in net boundary energy evaluation:
```python
# Historical bug in exact_energy_balance_audit.py (now resolved):
enet = ein - eout  # eout was already negative in OpenFOAM wallHeatFlux.dat!
```
In OpenFOAM's `wallHeatFlux` function object, heat entering the fluid domain is signed positive, while heat leaving across the isothermal cold boundary ($T_{\text{fluid}} > T_{\text{cold}} = 301.45\text{ K}$) produces a negative patch integral ($\dot{Q}_{\text{cold, raw}} < 0$). Subtracting a negative cumulative energy (`ein - (-eout) = ein + eout`) double-counted outward cold-wall heat loss as an artificial heat injection into the domain. At 17 min, cumulative cold-wall loss is $246.61\text{ J}$; adding rather than subtracting this loss inflated apparent net input to $8650.34\text{ J}$, producing the spurious $-507.28\text{ J}$ ($-5.86\%$) deficit against domain stored energy ($8143.07\text{ J}$).

When the physically correct definition is enforced ($E_{\text{net, in}} = E_{\text{hot}} - E_{\text{cold, out}} = E_{\text{hot}} + E_{\text{cold, raw}}$):
- At $t = 2\text{ min}$: $E_{\text{net, in}} = 1887.30 - 6.99 = 1880.31\text{ J}$, $E_{\text{stored}} = 1878.37\text{ J} \implies \Delta E = \mathbf{-1.94\text{ J}}$ ($\mathbf{-0.10\%}$).
- At $t = 6\text{ min}$: $E_{\text{net, in}} = 3729.13 - 48.24 = 3680.89\text{ J}$, $E_{\text{stored}} = 3676.92\text{ J} \implies \Delta E = \mathbf{-3.97\text{ J}}$ ($\mathbf{-0.11\%}$).
- At $t = 17\text{ min}$: $E_{\text{net, in}} = 8395.60 - 246.61 = 8148.99\text{ J}$, $E_{\text{stored}} = 8143.07\text{ J} \implies \Delta E = \mathbf{-5.92\text{ J}}$ ($\mathbf{-0.073\%}$).

The corrected script is committed in `validation/gallium2D/analysis/energy_balance/exact_energy_balance_audit.py`.

---

### 2.3 Mathematical Definition of the Step-Resolved Energy Audit
The authoritative energy audit evaluates boundary fluxes at every time step ($\Delta t = 0.1\text{ s}$, 10,200 steps to 17 min) against domain-integrated sensible and latent enthalpy:
$$E_{\text{net, in}}(t_n) = E_{\text{net}, 0 \to t_1} + \sum_{k=2}^n \frac{1}{2}\left[\dot{Q}_{\text{net}}(t_{k-1}) + \dot{Q}_{\text{net}}(t_k)\right]\Delta t_k$$
where $\dot{Q}_{\text{net}}(t) = \dot{Q}_{\text{hot}}(t) - \dot{Q}_{\text{cold, out}}(t)$. For the unrecorded initial startup step ($0 \le t \le t_1 = 0.1\text{ s}$), gross hot-wall input is evaluated as a rectangular step:
$$E_{\text{hot}, 0 \to 0.1} = \dot{Q}_{\text{hot}}(t_1)\Delta t_1$$
yielding $7.173\text{ J}$ for baseline EP ($\dot{Q}_{\text{hot}} = 71.73\text{ W}$) and $8.123\text{ J}$ for EHCDarcy ($\dot{Q}_{\text{hot}} = 81.23\text{ W}$). At $t_1 = 0.1\text{ s}$, the hot-wall thermal disturbance has not produced a meaningful cold-wall response on this timescale, and the computed cold-wall heat flux is negligible at $t = 0.1\text{ s}$ ($\dot{Q}_{\text{cold, out}}(0.1\text{ s}) \approx 4.5 \times 10^{-11}\text{ W}$, cumulative loss $< 10^{-11}\text{ J}$), so $E_{\text{net}, 0 \to 0.1}$ matches gross hot-wall startup delivery to within $< 10^{-11}\text{ J}$.

Stored enthalpy relative to initial conditions ($T_{\text{init}} = 301.45\text{ K}$, $f_L = 0$) evaluates to:
$$E_{\text{sens}}(t) = \int_\Omega \rho C_p (T - T_{\text{init}})\,dV, \quad E_{\text{lat}}(t) = \int_\Omega \rho L f_L\,dV, \quad E_{\text{stored}}(t) = E_{\text{sens}}(t) + E_{\text{lat}}(t)$$
Energy conservation closure error is defined as:
$$\Delta E(t) = E_{\text{stored}}(t) - E_{\text{net, in}}(t), \quad \epsilon_{\text{energy}}(t) = \frac{\Delta E(t)}{E_{\text{net, in}}(t)} \times 100\%$$

---

### 2.4 Time-Resolved Authoritative Energy Audit History

#### Complete Transient History: `EHCDarcy` Formulation ($42 \times 32$ Grid, $n_{\text{outer}}=8$, $\Delta t = 0.1\text{ s}$)
| Time | $\dot{Q}_{\text{hot}}$ (W) | $\dot{Q}_{\text{cold}}$ (W) | $E_{\text{hot}}$ (J) | $E_{\text{cold, out}}$ (J) | $E_{\text{net, in}}$ (J) | $E_{\text{sens}}$ (J) | $E_{\text{lat}}$ (J) | $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Closure Error (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1 min** | $11.630$ | $0.0578$ | $1294.0$ | $0.98$ | $1293.0$ | $59.1$ | $1233.3$ | $1292.4$ | $-0.65\text{ J}$ | **$-0.05\%$** |
| **2 min** | $8.802$ | $0.1318$ | $1887.3$ | $6.99$ | $1880.3$ | $81.4$ | $1797.0$ | $1878.4$ | $-1.94\text{ J}$ | **$-0.10\%$** |
| **4 min** | $7.503$ | $0.1748$ | $2848.7$ | $26.00$ | $2822.7$ | $107.3$ | $2712.3$ | $2819.5$ | $-3.12\text{ J}$ | **$-0.11\%$** |
| **6 min** | $7.215$ | $0.1959$ | $3729.1$ | $48.24$ | $3680.9$ | $129.1$ | $3547.9$ | $3676.9$ | $-3.95\text{ J}$ | **$-0.11\%$** |
| **8 min** | $7.130$ | $0.2210$ | $4589.4$ | $73.20$ | $4516.2$ | $149.5$ | $4362.2$ | $4511.7$ | $-4.49\text{ J}$ | **$-0.10\%$** |
| **10 min** | $7.057$ | $0.2536$ | $5442.7$ | $101.59$ | $5341.2$ | $169.5$ | $5166.7$ | $5336.2$ | $-5.01\text{ J}$ | **$-0.09\%$** |
| **12 min** | $7.095$ | $0.2957$ | $6289.7$ | $134.42$ | $6155.3$ | $187.9$ | $5962.1$ | $6150.0$ | $-5.36\text{ J}$ | **$-0.09\%$** |
| **14 min** | $7.033$ | $0.3514$ | $7133.7$ | $173.06$ | $6960.6$ | $206.7$ | $6748.3$ | $6955.0$ | $-5.65\text{ J}$ | **$-0.08\%$** |
| **16 min** | $7.000$ | $0.4281$ | $7974.8$ | $219.50$ | $7755.3$ | $224.3$ | $7525.2$ | $7749.5$ | $-5.80\text{ J}$ | **$-0.07\%$** |
| **17 min** | $6.992$ | $0.4782$ | $8395.60$ | $246.61$ | $8148.99$ | $232.71$ | $7910.36$ | $8143.07$ | $-5.92\text{ J}$ | **$-0.073\%$** |

#### Complete Transient History: `enthalpyPorosity` (Default: $\lambda = 1.0, n_{\text{outer}}=8$, $\Delta t = 0.1\text{ s}$)
| Time | $\dot{Q}_{\text{hot}}$ (W) | $\dot{Q}_{\text{cold}}$ (W) | $E_{\text{hot}}$ (J) | $E_{\text{cold, out}}$ (J) | $E_{\text{net, in}}$ (J) | $E_{\text{sens}}$ (J) | $E_{\text{lat}}$ (J) | $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Closure Error (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1 min** | $11.668$ | $0.0573$ | $1304.88$ | $0.97$ | $1303.91$ | $59.35$ | $1236.91$ | $1296.26$ | $-7.65\text{ J}$ | **$-0.59\%$** |
| **2 min** | $8.824$ | $0.1317$ | $1898.22$ | $7.01$ | $1891.21$ | $81.64$ | $1804.69$ | $1886.33$ | $-4.88\text{ J}$ | **$-0.26\%$** |
| **4 min** | $7.509$ | $0.1755$ | $2860.25$ | $26.13$ | $2834.11$ | $107.66$ | $2725.59$ | $2833.25$ | $-0.86\text{ J}$ | **$-0.03\%$** |
| **6 min** | $7.218$ | $0.1969$ | $3741.72$ | $48.49$ | $3693.23$ | $129.40$ | $3564.23$ | $3693.64$ | $+0.41\text{ J}$ | **$+0.01\%$** |
| **8 min** | $7.132$ | $0.2222$ | $4603.00$ | $73.56$ | $4529.44$ | $149.92$ | $4379.05$ | $4528.96$ | $-0.48\text{ J}$ | **$-0.01\%$** |
| **10 min** | $7.059$ | $0.2549$ | $5457.37$ | $102.05$ | $5355.32$ | $169.75$ | $5183.13$ | $5352.87$ | $-2.45\text{ J}$ | **$-0.05\%$** |
| **12 min** | $7.094$ | $0.2971$ | $6305.90$ | $134.96$ | $6170.94$ | $188.36$ | $5977.75$ | $6166.11$ | $-4.83\text{ J}$ | **$-0.08\%$** |
| **14 min** | $7.032$ | $0.3528$ | $7151.05$ | $173.65$ | $6977.40$ | $206.94$ | $6762.95$ | $6969.89$ | $-7.50\text{ J}$ | **$-0.11\%$** |
| **16 min** | $6.999$ | $0.4294$ | $7993.62$ | $220.10$ | $7773.52$ | $224.44$ | $7538.76$ | $7763.21$ | $-10.31\text{ J}$ | **$-0.13\%$** |
| **17 min** | $6.991$ | $0.4795$ | $8414.69$ | $247.19$ | $8167.51$ | $232.79$ | $7922.90$ | $8155.68$ | $-11.82\text{ J}$ | **$-0.14\%$** |

*Latent Enthalpy Reconciliation at $t = 17\text{ min}$*:  
For EHCDarcy, $f_L = 0.5736135 \implies E_{\text{lat}} = m L f_L = 13790.40 \times 0.5736135 = \mathbf{7910.36\text{ J}}$. Direct cell summation yields $\sum \rho V_{\text{cell}} L \alpha = 7910.36\text{ J}$.  
For EP ($\lambda = 1.0, n_{\text{outer}}=8$), $f_L = 0.5745226 \implies E_{\text{lat}} = 13790.40 \times 0.5745226 = \mathbf{7922.90\text{ J}}$. Direct cell summation yields $7922.90\text{ J}$. Both formulations reconcile latent enthalpy with liquid fraction to $0.01\text{ J}$.

---

### 2.5 Systematic Study of Newton Relaxation ($\lambda$) and Outer Iterations in Enthalpy-Porosity
In OpenFOAM's `enthalpyPorosityPhaseChangeModel.C`, the liquid fraction is updated after the sensible enthalpy matrix solve via Newton projection:
$$f_{\text{new}} = \min\left(\max\left(f + \lambda \frac{T - T_f}{L/C_p + \Delta T_{\text{window}}}, 0\right), 1\right)$$
where $\lambda$ is the Newton relaxation factor (code default: `1.0`, the exact Newton step for a linear melting window). In the baseline tutorial case dictionary, this parameter was set to `lambda 0.2`, damping each outer iteration correction to only $20\%$ of the full Newton step.

To decouple physical modeling capabilities from under-relaxation artifacts, a systematic sweep of $\lambda \in \{1.0, 0.5, 0.2\}$ across PIMPLE outer correctors $n_{\text{outer}} \in \{2, 3, 5, 8\}$ was executed to $t = 120\text{ s}$ ($2\text{ min}$) using the committed runner `validation/gallium2D/analysis/energy_balance/run_ep_relaxation_study.py`:

| Model Formulation | $\lambda$ | $n_{\text{outer}}$ | $f_L(2\text{ min})$ | Net Input $E_{\text{net}}$ (J) | Stored $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Closure Error (%) | Clock Time (s) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **EHCDarcy (Baseline)** | — | **8** | **0.1303** | **1880.31 J** | **1878.37 J** | **-1.94 J** | **-0.10%** | **35.0 s** |
| Enthalpy-Porosity | 1.0 | 2 | 0.1650 | 1514.33 J | 2370.39 J | +856.06 J | +56.53% | 9.7 s |
| Enthalpy-Porosity | 1.0 | 3 | 0.1358 | 1761.65 J | 1956.42 J | +194.77 J | +11.06% | 14.0 s |
| Enthalpy-Porosity | 1.0 | 5 | 0.1235 | 1924.17 J | 1779.98 J | -144.19 J | -7.49% | 22.9 s |
| **Enthalpy-Porosity** | **1.0** | **8** | **0.1309** | **1891.21 J** | **1886.33 J** | **-4.88 J** | **-0.26%** | **34.5 s** |
| Enthalpy-Porosity | 0.5 | 2 | 0.1750 | 1426.91 J | 2513.52 J | +1086.61 J | +76.15% | 9.5 s |
| Enthalpy-Porosity | 0.5 | 3 | 0.1480 | 1637.79 J | 2130.52 J | +492.73 J | +30.08% | 13.9 s |
| Enthalpy-Porosity | 0.5 | 5 | 0.1334 | 1806.34 J | 1922.33 J | +115.99 J | +6.42% | 22.4 s |
| Enthalpy-Porosity | 0.2 | 2 | 0.1762 | 1347.75 J | 2533.71 J | +1185.96 J | +88.00% | 9.8 s |
| Enthalpy-Porosity | 0.2 | 3 | 0.1527 | 1550.06 J | 2198.28 J | +648.22 J | +41.82% | 13.8 s |
| Enthalpy-Porosity | 0.2 | 5 | 0.1395 | 1710.51 J | 2009.79 J | +299.28 J | +17.50% | 22.2 s |
| Enthalpy-Porosity (Tutorial) | 0.2 | 8 | 0.1340 | 1802.03 J | 1931.11 J | +129.09 J | +7.16% | 35.0 s |
| Enthalpy-Porosity | 0.2 | 12 | 0.1317 | 1849.18 J | 1898.69 J | +49.51 J | +2.68% | 52.8 s |

**Key Findings on Enthalpy-Porosity Formulation**:
1. **The Apparent Non-Conservation Was an Under-Relaxation Artifact**:  
   Under tutorial settings ($\lambda = 0.2$), error decayed geometrically ($88\% \to 42\% \to 17.5\% \to 7.2\% \to 2.7\%$) because each pass applied only $20\%$ of the required Newton correction. Setting $\lambda = 1.0$ (exact Newton step) allows Enthalpy-Porosity to reach **$-0.26\%$** closure at $n_{\text{outer}} = 8$, easily passing the single $\le 2.0\%$ target and matching EHCDarcy's $-0.10\%$ closure.
2. **Enthalpy-Porosity Still Requires $\ge 6\text{--}8$ Correctors**:  
   At small corrector counts ($n_{\text{outer}} \le 5$), Enthalpy-Porosity fails the $\le 2.0\%$ target across all tested relaxation factors ($+56.5\%$ at $n=2$, $+11.1\%$ at $n=3$, $-7.49\%$ at $n=5$ for $\lambda = 1.0$). Because latent heat is introduced via an explicit RHS source term ($S_h = L \frac{\partial(\rho f)}{\partial t}$) with $f$ projected after the matrix solve, updating $f$ changes the thermal field for subsequent iterations; at low corrector counts, iteration lag during rapid conduction surges prevents discrete energy closure.
3. **Formulation Trade-off**:  
   Enthalpy-Porosity is **not** a cheaper alternative at low corrector counts because it requires $n_{\text{outer}} \ge 8$ to conserve energy. However, when run with code-default $\lambda = 1.0$ and $n_{\text{outer}} = 8$, Enthalpy-Porosity achieves full code-to-code agreement with `EHCDarcy`:
   - Liquid fraction at 17 min: $f_L = 0.5745$ (EP) vs. $0.5736$ (EHCDarcy), a difference of only $0.09$ percentage points.
   - Front position at mid-height ($y/H=0.5$): $56.91\text{ mm}$ (EP) vs. $57.06\text{ mm}$ (EHCDarcy), agreeing within $0.15\text{ mm}$.
   - Energy closure at 17 min: $-0.14\%$ (EP) vs. $-0.073\%$ (EHCDarcy).

---

## 3. Temporal Discretization Sensitivity & Convergence Study (EHCDarcy Formulation)

### 3.1 Brent et al. (1988) Time-Step Schedule Sensitivity
Brent, Voller, and Reid (1988, Section 4.1, p. 308) documented:
> *"A time step of 5 s was used for the first 20 s of the simulation and 10 s thereafter."*

To test the sensitivity of early-time melting to coarse time steps in the `EHCDarcy` formulation ($42 \times 32$ grid, $t \le 120\text{ s}$), runs were conducted using backward Euler time stepping:

| Time-Stepping Strategy | Time Step $\Delta t$ | Total Steps to 120 s | Computed $f_L(2\text{ min})$ | Difference vs. Brent ($0.1102$) |
| :--- | :---: | :---: | :---: | :---: |
| Fine Time Step | $0.1\text{ s}$ | 1200 | $0.1303$ | $+0.0201$ ($+18.2\%$) |
| Moderate Time Step | $1.0\text{ s}$ | 120 | $0.1288$ | $+0.0186$ ($+16.9\%$) |
| Coarse Constant Step | $5.0\text{ s}$ | 24 | $0.1071$ | $-0.0031$ ($-2.8\%$) |
| **Brent Exact Schedule** ($5\text{ s} \to 10\text{ s}$) | $5\text{ s} \ (0\text{--}20\text{s}) + 10\text{ s} \ (20\text{--}120\text{s})$ | **14** | **$0.0824$** | **$-0.0278$ ($-25.2\%$)** |

**Interpretation**:  
While a constant $5\text{ s}$ step fortuitously produces $f_L = 0.1071$ (close to Brent's 0.1102), executing Brent's exact 14-step schedule yields $f_L = 0.0824$ (25.2% below 0.1102). This demonstrates that early-time melting is highly sensitive to temporal discretization. However, the schedule alone does not reproduce Brent's original 1988 result, because Brent et al. used a 2D SIMPLER formulation with distinct under-relaxation factors, cell-face interpolations, and iterative convergence thresholds.

### 3.2 Formal Temporal Convergence Study ($\Delta t = 0.1\text{ s}, 0.05\text{ s}, 0.025\text{ s}$)
To establish whether the $\Delta t = 0.1\text{ s}$ solution is temporally resolved for the `EHCDarcy` formulation, systematic step halving was performed on the $42 \times 32$ mesh to $t = 120\text{ s}$:

| Time Step $\Delta t$ | Total Steps | Computed $f_L(2\text{ min})$ | Mid-Height Front $x_{\text{front}}(y/H=0.5)$ | Stored Enthalpy $E_{\text{stored}}$ |
| :---: | :---: | :---: | :---: | :---: |
| **$0.100\text{ s}$** | 1200 | **$0.1303$** | **$11.12\text{ mm}$** | $1878.37\text{ J}$ |
| **$0.050\text{ s}$** | 2400 | **$0.1304$** | **$11.13\text{ mm}$** | $1879.42\text{ J}$ |
| **$0.025\text{ s}$** | 4800 | **$0.1304$** | **$11.13\text{ mm}$** | $1879.90\text{ J}$ |

**Conclusion & Scope Limitation**:  
Halving the time step from $0.1\text{ s} \to 0.05\text{ s} \to 0.025\text{ s}$ in `EHCDarcy` changes $f_L(2\text{ min})$ by only **$0.0001$ ($0.07\%$)** and the phase front position by **$0.01\text{ mm}$**. This confirms temporal convergence of these reported quantities at early times ($t = 120\text{ s}$, $2\text{ min}$) under consistent numerical tolerances.

*Scope Limitation*: Because the full benchmark extends to $1020\text{ s}$ ($17\text{ min}$) where buoyancy-driven natural convection is fully developed and the discrete energy-balance discrepancy is largest, this early-time study does not establish temporal convergence for the entire 17-minute solution. Late-time temporal convergence in the convective regime remains unassessed and represents an avenue for future dedicated refinement studies.

---

## 4. Latent Advection Term Analysis: `div(phi, phaseFraction)`

Brent et al. (1988, p. 301, Eq. 7) stated:
> *"In the case of an isothermal phase change, the term $\mathrm{div}(\rho \mathbf{u} \Delta H)$ vanishes..."*

### 4.1 Physical Sharp Interface vs. Numerical Mushy Zone
A critical distinction must be drawn between the physical continuum problem and its fixed-grid numerical realization:

1. **Physical Isothermal Phase Change (Zero-Thickness Interface)**:  
   In a pure metal melting isothermally at temperature $T_m$ ($\Delta T_{\text{melt}} \to 0$), the phase boundary is a sharp surface:
   - In the liquid domain, $\Delta H = L = \text{constant}$. For an incompressible fluid ($\nabla \cdot \mathbf{u} = 0$), $\nabla \cdot (\rho \mathbf{u} L) = \rho L (\nabla \cdot \mathbf{u}) \equiv 0$.
   - Across the sharp interface, jump conditions enforce no-slip ($\mathbf{u} = 0$ on the solidus). Hence, physically, latent heat is released/absorbed solely through interface propagation (Stefan condition) rather than bulk convective transport of partially melted metal.

2. **Numerical Mushy Zone ($0 < f_L < 1$) & Latent Advection Diagnostic**:  
   In fixed-grid enthalpy formulations with a numerical smoothing window $\Delta T_{\text{melt}}$, a narrow artificial "mushy zone" spans $2\text{--}3$ grid cells where $\nabla f_L \neq 0$. In the continuum formulation, under constant density $\rho$, constant latent heat $L$, and incompressibility ($\nabla \cdot \mathbf{u} = 0$), the latent advection divergence expands mathematically to:
   $$\nabla \cdot (\rho \mathbf{u} L f_L) = \rho L (\mathbf{u} \cdot \nabla f_L) + \rho L f_L (\nabla \cdot \mathbf{u}) = \rho L (\mathbf{u} \cdot \nabla f_L)$$
   Because the magnitude of this term depends both on fluid velocity and on the phase-fraction gradient $\nabla f_L$, low velocity alone does not guarantee that the term is negligible everywhere.

   **Methodological Distinction: Continuum Diagnostic vs. Discrete Finite-Volume Operator**:  
   It is vital to establish the exact numerical operator evaluated:
   - *Continuum-Expanded Cell-Centred Diagnostic*: The results reported below are computed from cell-centred fields $\mathbf{u}_i$ and reconstructed cell-centred phase-fraction gradients $\nabla f_{L,i}$ via:
     $$I_{\text{signed}} = \sum_{i \in \mathcal{S}} V_i \rho L (\mathbf{u}_i \cdot \nabla f_{L,i}), \quad I_{\text{abs}} = \sum_{i \in \mathcal{S}} V_i \rho L |\mathbf{u}_i \cdot \nabla f_{L,i}|$$
     where the evaluation mask $\mathcal{S}$ is explicitly defined by the instantaneous mushy-zone cell set at $t = 1020\text{ s}$ ($0.01 < f_L(1020\text{ s}) < 0.99$, containing 85 cells; partitioned into the core mushy region $0.01 < f_L(1020\text{ s}) < 0.90$ with 80 cells and the outer interfacial fringe $0.90 \le f_L(1020\text{ s}) < 0.99$ with 5 cells).
   - *Discrete Finite-Volume Face-Flux Operator*: In OpenFOAM, the finite-volume divergence of a conserved face-flux quantity is obtained from the signed sum of face volume fluxes divided by cell volume:
     $$\left[\nabla_h \cdot (\boldsymbol{\phi} f_L)\right]_i = \frac{1}{V_i} \sum_{f \in \partial i} \phi_{i,f} f_{L,f} \quad [\text{s}^{-1}]$$
     and, for constant $\rho$ and $L$, the divergence of latent enthalpy flux per unit volume is:
     $$\left[\nabla_h \cdot (\rho L \boldsymbol{\phi} f_L)\right]_i = \frac{\rho L}{V_i} \sum_{f \in \partial i} \phi_{i,f} f_{L,f} \quad [\text{W/m}^3]$$
     where $\phi_{i,f} = \mathbf{u}_f \cdot \mathbf{S}_{i,f}$ denotes the face volume flux (in $\text{m}^3/\text{s}$, with outward-facing normal relative to cell $i$), and $f_{L,f}$ is interpolated to cell faces using the chosen convective discretization scheme (e.g., `upwind` or `linear`). The cell-integrated contribution to latent-advection power (in Watts) is then:
     $$V_i \left[\nabla_h \cdot (\rho L \boldsymbol{\phi} f_L)\right]_i = \rho L \sum_{f \in \partial i} \phi_{i,f} f_{L,f} \quad [\text{W}]$$
   Due to face-flux interpolation, cell gradient reconstruction, and discrete continuity residuals ($\sum_{f \in \partial i} \phi_{i,f} \ne 0$ at linear-solver tolerance), the cell-centred sum $\sum_i V_i \rho L (\mathbf{u}_i \cdot \nabla f_{L,i})$ is not strictly identical to the discrete finite-volume divergence.  
   For the whole closed domain $\Omega$ with impermeable walls ($\phi_{i,f} = 0$ on $\partial \Omega$), the discrete volume integral of a conservative finite-volume divergence algebraically reduces to the net boundary flux, which vanishes identically by interior face cancellation ($\sum_{i \in \Omega} V_i \left[\nabla_h \cdot (\rho L \boldsymbol{\phi} f_L)\right]_i = \rho L \sum_{i \in \Omega} \sum_{f \in \partial i} \phi_{i,f} f_{L,f} \equiv 0$). In contrast, sums of cell-centred products may retain small discrete reconstruction errors. Neither operator integrates to zero over an arbitrary mushy subregion $\mathcal{S}$.  
   Therefore, we describe the reported figures ($+3.01\text{ W}$, $-7.52\text{ W}$, $-4.52\text{ W}$) as a **continuum-expanded, cell-centred scale diagnostic** characterizing the order-of-magnitude intensity of advective transport, rather than definitive measurements of the discrete finite-volume face-flux operator.

   - **Signed Volume Integral** ($I_{\text{signed}} = \sum_{i \in \mathcal{S}} V_i \rho L (\mathbf{u}_i \cdot \nabla f_{L,i})$):  
     Under the continuum incompressibility identity, by Gauss's divergence theorem $\int_{\mathcal{S}} \nabla \cdot (\rho L \mathbf{u} f_L)\,dV = \int_{\partial \mathcal{S}} \rho L f_L (\mathbf{u} \cdot \mathbf{n}_{\text{out}})\,dA$. Thus, under continuum assumptions, $I_{\text{signed}} > 0$ is consistent with net outward advective flux (net divergence) across the boundary of subregion $\mathcal{S}$, whereas $I_{\text{signed}} < 0$ is consistent with net advective **flux convergence (net influx)** of latent enthalpy into subregion $\mathcal{S}$ from the adjacent fluid domain (in Watts).
   - **Absolute-Volume Integral** ($I_{\text{abs}} = \sum_{i \in \mathcal{S}} V_i \rho L |\mathbf{u}_i \cdot \nabla f_{L,i}|$):  
     Measures the total volume-integrated intensity of local advective redistribution (in Watts), independent of spatial sign cancellations.
   - **Transient Latent Uptake/Storage Rates**:  
     To make the comparison mathematically consistent, we evaluate both the signed net transient rate and the absolute transient rate:
     $$P_{\text{trans, signed}} = \sum_{i \in \mathcal{S}} V_i \rho L \left(\frac{\partial f_L}{\partial t}\right)_i, \quad P_{\text{trans, abs}} = \sum_{i \in \mathcal{S}} V_i \rho L \left|\frac{\partial f_L}{\partial t}\right|_i$$
     evaluated over the output interval $960 \to 1020\text{ s}$ ($\Delta t = 60\text{ s}$) across all cells belonging to the mushy-zone mask $\mathcal{S}$ defined at $t = 1020\text{ s}$. Positive rates represent melting (latent heat uptake/storage). The computed per-cell phase-fraction differences between 960 s and 1020 s are nonnegative across all cells in this evaluated mask ($\Delta f_{L,i} = f_{L,i}(1020\text{ s}) - f_{L,i}(960\text{ s}) \ge 0$ for all $i \in \mathcal{S}$). Consequently, the signed and absolute interval-averaged transient latent uptake rates coincide for this diagnostic interval:
     $$P_{\text{trans, signed}} \equiv P_{\text{trans, abs}} = 2.61\text{ W} \ (\text{core}), \quad 0.47\text{ W} \ (\text{fringe}), \quad 3.08\text{ W} \ (\text{entire mushy zone})$$

   | Subregion | Cell Count | Signed Advection $I_{\text{signed}}$ | Absolute Advection $I_{\text{abs}}$ | Signed Transient $P_{\text{trans, signed}}$ | Absolute Transient $P_{\text{trans, abs}}$ | Ratio ($I_{\text{abs}} / P_{\text{trans, abs}}$) |
   | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
   | **Core Mushy Region** ($0.01 < f_L < 0.90$) | 80 cells ($94.1\%$) | **$+3.01\text{ W}$** | **$3.69\text{ W}$** | **$2.61\text{ W}$** | **$2.61\text{ W}$** | **$1.41$** |
   | **Outer Interfacial Fringe** ($0.90 \le f_L < 0.99$) | 5 cells ($5.9\%$) | **$-7.52\text{ W}$** | **$9.20\text{ W}$** | **$0.47\text{ W}$** | **$0.47\text{ W}$** | **$19.6$** |
   | **Entire Mushy Zone** ($0.01 < f_L < 0.99$) | 85 cells ($100.0\%$) | **$-4.52\text{ W}$** | **$12.90\text{ W}$** | **$3.08\text{ W}$** | **$3.08\text{ W}$** | **$4.19$** |

   - **Integrated Magnitude Comparison vs. Cell-Wise Median in the Core**:  
     In the core mushy zone ($0.01 < f_L < 0.90$, 80 cells), Darcy damping strongly quenches flow, yielding a median velocity of **$1.62 \times 10^{-7}\text{ m/s}$**. The median of cell-wise ratios evaluates to:
     $$\operatorname{median}_{i \in \text{core}} \left( \frac{|\mathbf{u}_i \cdot \nabla f_{L,i}|}{|\partial f_{L,i}/\partial t|} \right) = \mathbf{0.99\%} \approx \mathbf{1.0\%} \quad (\text{ratio of median magnitudes: } 0.50\%)$$
     While this low median confirms that the *typical* core cell experiences negligible advection relative to local melting, the integrated core absolute advection ($3.69\text{ W}$) is **$1.41\times$** the integrated transient melting rate ($2.61\text{ W}$), and signed advection evaluates to $+3.01\text{ W}$ (consistent with net advective outflux under continuum assumptions). This disparity reveals that advection is highly non-uniform: a small subset of cells within the core possessing moderate residual flow and steep phase-fraction gradients contributes disproportionately to the total advective transport, precluding a claim that integrated core advection is negligible.
   - **Time-Level Disparity Note**:  
     The advection diagnostic is evaluated instantaneously at $t = 1020\text{ s}$, whereas the transient latent uptake rate $\frac{\partial f_L}{\partial t}$ is backward-differenced (averaged) over $960 \le t \le 1020\text{ s}$ ($\Delta t = 60\text{ s}$). This temporal disparity between an instantaneous velocity state and a 60-second interval-averaged phase change limits a strict term-by-term instantaneous balance, but provides an authoritative order-of-magnitude scale comparison.
   - **Outer Fringe Behavior ($f_L \ge 0.90$) & Cavity Power Scale Comparison**:  
     In the 5 outermost cells adjacent to the convective liquid boundary layer, fluid velocities are only partially damped ($\sim 2.3 \times 10^{-4}$ to $4.7 \times 10^{-3}\text{ m/s}$) while steep phase-fraction gradients persist ($|\nabla f_L| \sim 400\text{--}500\text{ m}^{-1}$). Here, the localized absolute advection integral is **$9.20\text{ W}$**, and the signed integral is **$-7.52\text{ W}$** (consistent with net advective influx / convergence from the liquid boundary layer into the fringe under continuum assumptions).  
     To place these numbers in context, at $t = 17\text{ min}$, cavity-level power scales are:
     - Instantaneous hot-wall heat transfer rate: $\dot{Q}_{\text{hot}} \approx \mathbf{6.99\text{ W}}$.
     - Domain-wide latent energy storage rate over $16 \to 17\text{ min}$: $\frac{\Delta E_{\text{lat}}}{\Delta t} = \frac{7910.4 - 7525.2}{60} \approx \mathbf{6.42\text{ W}}$.  
     Across the entire mushy zone (85 cells), absolute advection totals **$12.90\text{ W}$**, which is **$4.19\times$** the total mushy-zone transient uptake rate ($3.08\text{ W}$) and nearly double the cavity power input ($6.99\text{ W}$).
   - **Physical Role of Local Redistribution vs. Global Conservation**:  
     The reported ratios ($1.41$ in the core, $4.19$ domain-wide) are comparisons of **integrated local magnitudes**—not evidence that the omitted term would alter the global energy balance by those amounts. For an incompressible fluid in a closed domain with no-slip walls ($\mathbf{u} = 0$ on all boundaries), the global volume integral of advective divergence vanishes identically by Gauss's theorem:
     $$\int_\Omega \nabla \cdot (\rho \mathbf{u} L f_L)\,dV = \int_{\partial \Omega} \rho L f_L (\mathbf{u} \cdot \mathbf{n})\,dA \equiv 0$$
     Therefore, omitting the term does not introduce a net external energy sink or source for the domain as a whole. Rather, the term acts as an **internal spatial redistribution of latent enthalpy**. Under the continuum incompressibility identity, the negative fringe diagnostic ($I_{\text{signed}} = -7.52\text{ W}$) is consistent with net latent-enthalpy flux convergence, while the positive core diagnostic ($I_{\text{signed}} = +3.01\text{ W}$) is consistent with net flux divergence. These signs are indicative of the reconstructed cell-centred field; definitive confirmation of the regional flux directions requires evaluating the corresponding flux across the boundaries of each subregion using a consistent face-flux representation.
   - **Methodological Status of the Omission**:  
     Brent et al. (1988, Section 2, Eq. 7) justified dropping $\nabla \cdot (\rho \mathbf{u} \Delta H)$ under the assumption of an isothermal phase change, where a sharp physical interface with no-slip on the solidus strictly precludes fluid motion within the phase-change region. While this argument is exact for the continuum sharp-interface Stefan problem, **it does not establish that the term vanishes in a diffuse numerical interface of finite width**.  
     Consequently, omitting $\nabla \cdot (\rho \mathbf{u} L f_L)$ is an **established benchmark modeling approximation** standard in classical CFD phase-change literature (e.g., Brent et al. 1988, Hannoun et al. 2003, Ye & Arıcı 2023) that preserves consistency with published benchmark formulations. However, our diagnostic demonstrates that this omission leaves a potentially important localized redistribution across the outer fringe ($9.20\text{ W}$ absolute) and parts of the core ($3.69\text{ W}$ absolute) whose detailed quantitative effect on the local phase-front trajectory has not yet been isolated. A controlled comparison against a consistently implemented solver incorporating the full discrete finite-volume divergence term represents an important avenue for future research.

### 4.2 Implementation Verification in `phaseChangeMultiRegionFoam`
Inspection of `EEqn.H` confirms:
- The energy equation solves for sensible enthalpy $h$ with standard advection $\nabla \cdot (\rho \mathbf{u} h)$.
- Latent heat is governed solely by the unsteady source term $L \frac{\partial(\rho f_L)}{\partial t}$ in `enthalpyPorosity` (or via $C_p^{\text{eff}}$ in `EHCDarcy`).
- The entry `div(phi,phaseFraction)` previously present in `system/fluid/fvSchemes` was an unused legacy dictionary token that was never invoked during matrix assembly. Omitting numerical latent advection is an established and physically motivated benchmark approximation, while recognizing that in diffuse-interface methods it represents a modeling choice rather than an exact mathematical identity.

---

## 5. Spatial Grid Sensitivity & ASME / Roache GCI Analysis (EHCDarcy Formulation, to 17 min)

Three systematically refined meshes were executed to the full $17\text{ min}$ ($1020\text{ s}$) duration:
- **Grid 3 (Coarse)**: $32 \times 24$ ($768$ cells, $h_3 = 2.71\text{ mm}$)
- **Grid 2 (Medium)**: $42 \times 32$ ($1344$ cells, $h_2 = 2.05\text{ mm}$)
- **Grid 1 (Fine)**: $63 \times 48$ ($3024$ cells, $h_1 = 1.37\text{ mm}$)

Refinement ratios: $r_{21} = h_2 / h_1 = 1.500$, $r_{32} = h_3 / h_2 = 1.3223$.

### 5.1 Grid Convergence & Sensitivity Table

| Grid | Resolution | Total Cells | $f_L(10\text{ min})$ | $f_L(17\text{ min})$ | $x_{\text{front}}(17\text{m}, y/H=0.2)$ | $x_{\text{front}}(17\text{m}, y/H=0.5)$ | $x_{\text{front}}(17\text{m}, y/H=0.8)$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Coarse (Grid 3) | $32 \times 24$ | 768 | $0.3698$ | $0.5682$ | $29.18\text{ mm}$ | $56.01\text{ mm}$ | $66.75\text{ mm}$ |
| Medium (Grid 2) | $42 \times 32$ | 1344 | $0.3747$ | $0.5736$ | $28.19\text{ mm}$ | $57.06\text{ mm}$ | $67.80\text{ mm}$ |
| Fine (Grid 1) | $63 \times 48$ | 3024 | $0.3760$ | $0.5725$ | $27.58\text{ mm}$ | $57.68\text{ mm}$ | $66.87\text{ mm}$ |
| **Richardson Extrapolated ($\phi_{\text{ext}}$)** | — | $\infty$ | — | — | — | **$57.99\text{ mm}$** | — |

### 5.2 Formal ASME / Roache Grid Convergence Analysis
1. **Global Liquid Fraction $f_L(17\text{ min})$**:
   The computed sequence from coarse to fine is $0.5682 \to 0.5736 \to 0.5725$. Because the medium grid is slightly higher than the fine grid (by $0.0011$, or $0.19\%$), this sequence is non-monotonic. Consequently, a formal monotonic Richardson extrapolation cannot be claimed for $f_L(17\text{ min})$. The total spread across all three grids is only $0.0054$ ($0.95\%$), establishing that the global melt fraction is spatially insensitive to grid refinement beyond $42 \times 32$.
2. **Phase Front Position at Mid-Height $y/H = 0.5$**:
   The front coordinate sequence is strictly monotonic ($56.01\text{ mm} \to 57.06\text{ mm} \to 57.68\text{ mm}$). For unequal refinement ratios ($r_{21} = 1.500, r_{32} = 1.3223$), solving the generalized transcendental ASME equation:
   $$\frac{\epsilon_{32}}{\epsilon_{21}} = \frac{57.06 - 56.01}{57.68 - 57.06} = 1.6935 \implies p = 2.70$$
   - Observed order of convergence: $p = 2.70$
   - Richardson extrapolated front position: $\phi_{\text{ext}} = 57.99\text{ mm}$
   - Fine-grid relative error: $e_a^{21} = 1.07\%$
   - **Fine-Grid Convergence Index ($\text{GCI}_{\text{fine}}^{21}$)**: **$0.67\%$**
   - **Medium-Grid Convergence Index ($\text{GCI}_{\text{med}}^{32}$)**: **$2.04\%$**
     $$\mathrm{GCI}_{\mathrm{med}}^{32} = F_s \frac{|\phi_2 - \phi_3|/|\phi_2|}{r_{32}^p - 1} \times 100\% = 1.25 \frac{|57.06 - 56.01|/57.06}{1.3223^{2.70} - 1} \times 100\% = \mathbf{2.04\%}$$
   - **Asymptotic Range Consistency Ratio**:
     $$\frac{\mathrm{GCI}_{\mathrm{med}}^{32}}{r_{21}^p \, \mathrm{GCI}_{\text{fine}}^{21}} = \frac{2.038\%}{1.500^{2.70} \times 0.672\%} = \frac{2.038\%}{2.007\%} \approx \mathbf{1.01}$$
     The consistency ratio is within $1\%$ of the theoretical asymptotic limit of $1.00$, confirming that the mid-height front calculation is well within the asymptotic range of grid convergence.

**Scope and Interpretation of Grid Convergence**:  
A spread of $0.0054$ ($0.95\%$) across all three grids for $f_L(17\text{ min})$ demonstrates low spatial sensitivity, but does not by itself establish global grid independence. The formal ASME GCI of $0.67\%$ applies specifically to the mid-height phase front coordinate at $17\text{ min}$ ($x_{\text{front}}(17\text{m}, y/H=0.5)$), confirming that spatial discretization error for this key morphological metric is bounded and well resolved.

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

**Physical Insight**:  
$C_u = 1.6 \times 10^6\text{ kg/(m}^3\text{ s)}$ serves as a literature-motivated baseline originating from Brent et al. (1988). Lower $C_u$ permits residual fluid velocities in the mushy zone, increasing upper front tilt ($x_{0.8}/x_{0.2} = 2.74$). Increasing $C_u$ toward $10^8$ extinguishes interfacial velocity, leading to an asymptotic liquid fraction with only a $1.2\%$ change between $1.6 \times 10^6$ and $1.0 \times 10^8$.

### 6.2 Melting Window Width Sweep ($\Delta T_{\text{melt}} \in [0.1, 0.5]\text{ K}$)

| Melting Window $\Delta T_{\text{melt}}$ | $f_L(10\text{ min})$ | $f_L(17\text{ min})$ | $\overline{\text{Nu}}_{\text{hot}}(17\text{ min})$ | $x_{\text{front}}(y/H=0.2)$ | $x_{\text{front}}(y/H=0.8)$ |
| :--- | :---: | :---: | :---: | :---: | :---: |
| $0.1\text{ K}$ ($302.93 \pm 0.05\text{ K}$) | $0.3738$ | $0.5727$ | $4.57$ | $27.81\text{ mm}$ | $67.82\text{ mm}$ |
| $0.2\text{ K}$ ($302.93 \pm 0.10\text{ K}$) | $0.3736$ | $0.5728$ | $4.56$ | $27.90\text{ mm}$ | $67.81\text{ mm}$ |
| **$0.5\text{ K}$ ($302.93 \pm 0.25\text{ K}$, Baseline)** | **$0.3747$** | **$0.5736$** | **$4.50$** | **$28.19\text{ mm}$** | **$67.80\text{ mm}$** |

**Physical Insight**:  
Varying the numerical phase transition window from $0.5\text{ K}$ down to $0.1\text{ K}$ changes $f_L(17\text{ min})$ by only **$0.0009$ ($0.15\%$)**. The upper front position ($y/H=0.8$) changes by only $0.02\text{ mm}$ ($67.80 \to 67.82\text{ mm}$), while the lower front near the floor ($y/H=0.2$) shifts by $0.38\text{ mm}$ ($28.19 \to 27.81\text{ mm}$), confirming that $\Delta T = 0.5\text{ K}$ provides a stable representation of the isothermal phase boundary.

---

## 7. Reference Data Limits & Literature Benchmark Comparison

### 7.1 Distinguishing Numerical Benchmark Fronts from Experimental Measurements
To ensure scientific rigor, a clear distinction is maintained between code-to-code verification and physical experimental validation:
1. **Numerical Benchmark Contours (Brent, Voller, and Reid 1988)**:  
   All quantitative front morphology metrics in Section 1.2 (mean and max front deviations) are calculated directly against the digitized numerical simulation fronts published in Figure 9 of Brent et al. (1988). Agreement with these contours demonstrates consistent numerical implementation of convection-diffusion phase change.
2. **Digitization Uncertainty**:  
   The reference melt fronts digitized from Figure 9 of Brent et al. (1988) carry an intrinsic scanning/resolution uncertainty of approximately **$\pm 1.0\text{ mm}$** ($\pm 1.1\%$ of cavity width $W$) due to contour line width. A shaded uncertainty envelope of $\pm 1\text{ mm}$ is included in all validation plots.

### 7.2 Experimental Comparison with Gau & Viskanta (1986)
Gau and Viskanta (1986) performed the foundational physical experiments for gallium melting in this geometry. While qualitative agreement in front morphology is evident, exact quantitative matching with 2D idealized CFD is constrained by physical experimental realities:
1. **Wall Thermal Inertia**: In the experiment, the hot wall could not experience an instantaneous step change to $311.15\text{ K}$ due to heating jacket thermal inertia and early heat loss. Consequently, experimental melt fronts at $t \le 2\text{ min}$ lag behind idealized isothermal numerical simulations.
2. **Pour-Out Measurement Disruption**: Gau & Viskanta measured front positions by pouring out the molten gallium at selected intervals, solidifying a shell, and restarting the test, disrupting the developed natural convection boundary layer.
3. **3D Endwall Drag & Crystal Anisotropy**: Physical gallium crystals exhibit anisotropic thermal conductivity in the solid phase ($k$ varying from $28$ to $41\text{ W/(m K)}$ across crystallographic axes), and non-empty z-direction glass windows introduce side-wall viscous drag, explaining why 2D simulations predict slightly faster upper melting than raw experimental data.

### 7.3 Consistency Check with Ye & Arıcı (2023) Numerical Correlation
Ye & Arıcı (2023, Section 3.2.2, Eq. 6) conducted an independent 2D verification study on natural convection melting of gallium, correlating their grid-converged numerical liquid fractions as:
$$\phi_n = 2.708\,\tau^{0.712}\,\mathrm{Ra}^{0.0243}\,\mathrm{Ar}^{-0.0515}$$
where:
- $\phi_n$ represents the numerical global liquid volume fraction ($f_L$).
- $\tau = \mathrm{Fo} \cdot \mathrm{Ste}$, with $\mathrm{Fo} = \frac{\alpha t}{W^2} = 1.7761$ and $\mathrm{Ste} = \frac{C_p (T_h - T_m)}{L} = 0.03912$ (at $t = 1020\text{ s}$).
- $\mathrm{Ra} = \frac{g \beta (T_h - T_m) H^3}{\nu \alpha} = 6.063 \times 10^5$.
- $\mathrm{Ar} = \frac{H}{W} = 0.7143$.

Evaluating Eq. (6) via the committed verification script (`validation/gallium2D/analysis/evaluate_ye_arici_correlation.py`) at $t = 17\text{ min}$ ($1020\text{ s}$) yields:
$$\phi_n^{\text{Ye \& Arıcı (2023)}} = \mathbf{0.5703}$$
- `phaseChangeFoam` (EHCDarcy, $42 \times 32$): **$f_L = 0.5736$** (Discrepancy: **$+0.0033$, or $0.33$ percentage points**).
- `phaseChangeFoam` (EP, $\lambda = 1.0, n=8$): **$f_L = 0.5745$** (Discrepancy: **$+0.0042$, or $0.42$ percentage points**).
- `phaseChangeFoam` (EHCDarcy, $63 \times 48$): **$f_L = 0.5725$** (Discrepancy: **$+0.0022$, or $0.22$ percentage points**).
- `phaseChangeFoam` (EHCDarcy, $32 \times 24$): **$f_L = 0.5682$** (Discrepancy: **$-0.0021$, or $0.21$ percentage points**).

**Notes on Equation Provenance & Physical Scaling**:
1. *Formula Clarification*: An earlier draft cited a different correlation from another Ye & Arıcı paper ($0.03816\,\mathrm{Fo}^{0.4907}\mathrm{Ra}^{0.3164}\mathrm{Ste}^{0.8407}$) fitted to mean melt-layer thickness, which evaluates to $\sim 0.23\text{--}0.38$ for these dimensionless numbers. The authoritative benchmark liquid-fraction correlation is Eq. (6) above, which is fully documented and evaluated by the committed script.
2. *Dominance of Conduction/Stefan Scaling*: The exponent on the Rayleigh number is only $\mathrm{Ra}^{0.0243}$. Consequently, varying $\mathrm{Ra}$ by an entire order of magnitude ($10\times$) modifies the predicted liquid fraction by only $10^{0.0243} = 1.0575$ ($+5.75\%$). Agreement with this correlation primarily verifies the transient conduction and Stefan scaling ($\tau^{0.712}$), rather than fine convective boundary layer dynamics.
3. *Interpretation*: Agreement with a numerical correlation fitted to independent CFD simulations provides valuable code-to-code verification, but is strictly distinguished from experimental validation against physical measurements.

---

## 8. Summary Conclusions & Verification Scope Assessment

### 8.1 Summary of Methodological Resolutions
1. **Energy-Balance Audit and Boundary Accounting Resolutions**:
   - **Unified $\le 2.0\%$ Acceptance Target**: A single, strict energy closure criterion of $|\epsilon_{\text{energy}}| \le 2.0\%$ across all times ($t \le 17\text{ min}$) is enforced. `EHCDarcy` achieves complete energy closure across the entire transient ($-0.10\%$ at 2 min, $-0.11\%$ at 6 min, and $-0.073\%$ at 17 min).
   - **Root Cause of Historical Audit Negative Deficits**: In earlier versions of `exact_energy_balance_audit.py`, negative raw cold-wall heat flux was subtracted, double-counting cold-wall heat loss as an artificial heat injection and producing false deficits ($-1.26\%$ to $-5.86\%$). Correcting to $E_{\text{net}} = E_{\text{hot}} - E_{\text{cold, out}}$ revealed the true physical closure of $-0.10\%$ to $-0.073\%$.
   - **Enthalpy-Porosity Newton Relaxation ($\lambda$) Resolution**: The apparent failure of Enthalpy-Porosity to conserve energy ($+6.72\%$ at 2 min, $+5.18\%$ at 17 min) was an artifact of the tutorial's `lambda 0.2` under-relaxation. Setting $\lambda = 1.0$ (exact Newton step) allows EP with $n_{\text{outer}}=8$ to achieve $-0.26\%$ closure at 2 min, $+0.01\%$ at 6 min, and $-0.14\%$ at 17 min, matching EHCDarcy's global conservation and front morphology ($f_L = 0.5745$ vs $0.5736$). However, EP still fails at low outer correctors ($n_{\text{outer}} \le 5$), confirming that EHCDarcy remains the superior, more robust formulation when corrector count is constrained.
   - **Coarse Integration Artifact**: The coarse $1601.5\text{ J}$ script underestimated first-minute delivery by $229.66\text{ J}$ due to temporal undersampling of the early conduction surge, partially offset by a $14.80\text{ J}$ trapezoidal overestimation during minute 2. Temporal undersampling accounts for $99.97\%$ ($214.86\text{ J}$) of the coarse gross discrepancy, while flux estimator discretization accounts for only $0.03\%$ ($0.065\text{ J}$).

2. **Methodological Qualification of Latent Heat Advection**:
   - Evaluated as a cell-centred continuum scale diagnostic $\sum V_i \rho L (\mathbf{u}_i \cdot \nabla f_{L,i})$ across 85 mushy cells at $t = 1020\text{ s}$, the median cell-wise ratio in the core mushy zone is small ($0.99\% \approx 1.0\%$).
   - In the outer 5 fringe cells ($f_L \ge 0.90$), localized advection reaches $9.20\text{ W}$ absolute and $-7.52\text{ W}$ signed. Under continuum incompressibility, these signs represent local divergence in the core and convergence into the fringe, canceling across the closed domain. Omitting $\nabla \cdot (\rho \mathbf{u} L f_L)$ is correctly classified as an established benchmark modeling approximation standard across classical CFD phase-change literature.

3. **Spatial and Temporal Sensitivity Bounds**:
   - Mid-height front position at 17 min demonstrates strict monotonic spatial convergence with an ASME observed order of $p = 2.70$, an asymptotic consistency ratio of $1.01$, and a fine-grid $\text{GCI}$ of $0.67\%$. This grid-convergence verification is established for `EHCDarcy` only.
   - Global melt fraction at 17 min matches the modern grid-converged numerical correlation of Ye & Arıcı (2023, Eq. 6) within $0.33$ percentage points ($0.5736$ vs. $0.5703$).

### 8.2 Scope Distinction: Numerical Benchmark Verification vs. Physical Experimental Validation
To preserve scientific precision, the scope and status of this report are explicitly defined:
- **Numerical Benchmark Verification**: This report provides rigorous evidence of **numerical benchmark consistency and code-to-code verification** against the established 2D CFD literature (Brent et al. 1988, Hannoun et al. 2003, Ye & Arıcı 2023). Under the stated idealized assumptions, `phaseChangeMultiRegionFoam` with the `EHCDarcy` formulation achieves a mean phase-front deviation of $3.46\text{ mm}$ ($3.89\%$ of cavity width) and maximum local deviation of $6.34\text{ mm}$ ($7.13\%$ of cavity width) relative to Brent et al. (1988), with an energy closure discrepancy of $-0.073\%$.
- **Distinction from Experimental Validation**: These results **do not constitute a full quantitative validation against independent physical experimental measurements**. Quantitative alignment with physical experiments (such as Gau & Viskanta 1986) requires resolving 3D endwall boundary-layer drag, conjugate thermal inertia of the heating jacket and wall apparatus, solid-state gallium crystallographic anisotropy ($k \in [28, 41]\text{ W/(m K)}$), and non-intrusive front measurement methods that avoid pour-out boundary layer disruption.
- **Engineering Tolerances**: Acceptance thresholds (front mean deviation $\le 5\% W$, liquid fraction within $\pm 0.05$) are classified as engineering tolerances suitable for computational thermal engineering.
- **Document Status**: Accordingly, this study is formally presented as a **high-fidelity numerical benchmark verification report**, providing a solid foundation for future 3D conjugated experimental validation studies.

---

## Appendix A: Historical Coarse-Audit Reconstruction

To preserve complete arithmetic traceability for how the initial post-processing script (`energy_balance_audit.py`, archived in `energy_balance_results.txt`) arrived at its coarse estimate of $1601.5\text{ J}$, the exact algorithm and rounding transformations are documented below:

1. **First Interval ($0 \to 60\text{ s}$)**:  
   Sampling exclusively at 60 s intervals without $t = 0$ data, the script applied an ad-hoc factor of $1.5$ to account for early conduction:
   $$E_{0 \to 60} = \dot{Q}(60\text{ s}) \times 60\text{ s} \times 1.5$$
   Evaluating with the unrounded finite-difference flux sample ($Q_{\text{FD}} = 11.22090\text{ W}$) yields:
   $$E_{\text{FD}, 0 \to 60} = 11.22090 \times 60 \times 1.5 = 1009.881\text{ J} \quad (\text{displayed as } 1009.9\text{ J})$$
   Evaluating with the unrounded native patch flux ($Q_{\text{native}} = 11.22152\text{ W}$) yields:
   $$E_{\text{native}, 0 \to 60} = 11.22152 \times 60 \times 1.5 = 1009.937\text{ J} \quad (\text{displayed as } 1009.94\text{ J})$$

2. **Second Interval ($60 \to 120\text{ s}$)**:  
   Standard trapezoidal quadrature between the two samples yields:
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
   Comparing the step-resolved baseline energy ($1816.378\text{ J}$) against the displayed log total ($1601.5\text{ J}$) yields an apparent difference of $214.8780\text{ J} \approx 214.88\text{ J}$. This arithmetic identity decomposes exactly as:
   $$\begin{aligned}
   \Delta E_{\text{archived display}} &= \Delta E_{\text{native temporal}} + \Delta E_{\text{estimator}} - \Delta E_{\text{archive rounding}} \\
   &= 214.8622\text{ J} + 0.0648\text{ J} - 0.0490\text{ J} \\
   &= \mathbf{214.8780\text{ J}} \approx \mathbf{214.88\text{ J}}
   \end{aligned}$$
   where $-0.0490\text{ J} = 1601.451\text{ J} - 1601.500\text{ J}$ is the one-decimal log display offset.

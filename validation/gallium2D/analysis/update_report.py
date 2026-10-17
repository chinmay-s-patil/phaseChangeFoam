#!/usr/bin/env python3
"""
Updates ANALYSIS_REPORT.md with:
1. Unified single <= 2.0% energy closure criterion across all sections.
2. Complete Enthalpy-Porosity Newton relaxation (lambda) and outer-corrector study.
3. Cold-wall raw flux sign artifact resolution (EHCDarcy -1.26%, -2.87%, -4.39%).
4. Streamlined Section 2 with historical 4-decimal coarse arithmetic moved to Appendix A.
5. Ye & Arıcı (2023, Eq. 6) correlation scripted check, scaling note (Ra^0.0243), and formula clarification.
6. Scope wording: Numerical benchmark verification, engineering tolerances, and EHCDarcy-only grid study scope.
"""

import os
import re

report_path = "validation/gallium2D/analysis/ANALYSIS_REPORT.md"
with open(report_path, "r") as f:
    text = f.read()

# -----------------------------------------------------------------------------
# 1. Section 2: Replace lines between "## 2. Energy-Balance Audit" and "## 3. Temporal Discretization"
# -----------------------------------------------------------------------------
sec2_start = text.find("## 2. Energy-Balance Audit")
sec3_start = text.find("## 3. Temporal Discretization Sensitivity")

if sec2_start == -1 or sec3_start == -1:
    print("Could not find section markers for Section 2 and 3!")
    exit(1)

new_sec2 = """## 2. Energy-Balance Audit, Post-Processing Verification & Formulation Sensitivity

### 2.1 Root Cause of the Initial Coarse Audit Discrepancy (+20.6% Apparent Surplus)
The initial coarse audit script (`energy_balance_audit.py`, archived in `energy_balance_results.txt`) evaluated the baseline Enthalpy-Porosity simulation using 60-second sampling intervals and reported apparent surpluses of $+31.3\%$ at 1 min, $+20.6\%$ at 2 min ($1601.5\\text{ J}$ input vs. $1931.1\\text{ J}$ stored), and $+6.6\%$ at 17 min.

Controlled step-resolved investigation ($\\Delta t = 0.1\\text{ s}$, 1,200 steps to 2 min) demonstrates that this apparent non-conservation was a **post-processing artifact** driven by two primary factors:
1. **Severe Temporal Undersampling of the Early Conduction Surge**:  
   At $t = 0.1\\text{ s}$, transient conduction into the subcooled slab produces an initial heat flux peak of $\\dot{Q}_{\\text{hot}} \\approx 72\\text{--}81\\text{ W}$ ($q''_{\\text{wall}} \\approx 2.3\\text{--}2.6 \\times 10^5\\text{ W/m}^2$), decaying rapidly to $11.2\\text{ W}$ by 60 s. Sampling exclusively at 60 s missed this surge, underestimating first-minute gross input by **$229.66\\text{ J}$** ($1239.60\\text{ J}$ step-resolved vs. $1009.94\\text{ J}$ coarse). During the second minute ($60 \\to 120\\text{ s}$), diffusion decay is convex ($q \\propto t^{-1/2}$), causing the coarse trapezoid to overestimate delivery by **$14.80\\text{ J}$** ($591.58\\text{ J}$ coarse vs. $576.78\\text{ J}$ step-resolved). The net temporal undersampling error is:
   $$\\Delta E_{\\text{temporal}} = 229.66\\text{ J} - 14.80\\text{ J} = \\mathbf{214.86\\text{ J}} \\quad (99.97\\%\\text{ of gross difference})$$
2. **Estimator Discretization Residual**:  
   At the 60 s and 120 s output timestamps, the finite-difference wall gradient ($Q_{\\text{FD}}$) and OpenFOAM's native patch integral ($Q_{\\text{native}}$) agree to within $< 0.01\\%$ ($11.2215\\text{ W}$ vs. $11.2209\\text{ W}$ at 60 s; $8.4978\\text{ W}$ vs. $8.4981\\text{ W}$ at 120 s). The difference between the coarse native reconstruction ($1601.52\\text{ J}$) and unrounded finite-difference calculation ($1601.45\\text{ J}$) is only **$0.065\\text{ J}$** ($0.03\\%$ of the gross difference).

| Contribution to Gross Discrepancy ($t = 2\\text{ min}$) | Energy (J) | Fraction of Difference |
| :--- | :---: | :---: |
| **Total Gross Difference** ($1816.38\\text{ J} - 1601.45\\text{ J}$) | **$214.93\\text{ J}$** | **$100.0\\%$** |
| **Temporal Undersampling Error** ($E_{\\text{native, step}} - E_{\\text{native, coarse}}$: $1816.38 - 1601.52$) | **$214.86\\text{ J}$** | **$99.97\\%$** |
| **Estimator Discretization Residual** ($E_{\\text{native, coarse}} - E_{\\text{FD, coarse}}$: $1601.52 - 1601.45$) | **$0.07\\text{ J}$** | **$0.03\\%$** |

*(The complete step-by-step arithmetic derivation reconciling the archived log display value of $1601.5\\text{ J}$ is documented in [Appendix A](#appendix-a-historical-coarse-audit-reconstruction)).*

---

### 2.2 Resolution of Historical Post-Processing Artifacts (EHCDarcy -1.26%, -2.87%, -4.39%)
In an earlier version of the time-resolved audit script (`exact_energy_balance_audit.py`), `EHCDarcy` was reported with negative closure discrepancies of $-1.26\\%$ at 2 min, $-2.87\\%$ at 6 min, $-4.39\\%$ at 12 min, and $-5.86\\%$ at 17 min. 

Inspection of the script source revealed a sign convention error in net boundary energy evaluation:
```python
# Historical bug in exact_energy_balance_audit.py (now resolved):
enet = ein - eout  # eout was already negative in OpenFOAM wallHeatFlux.dat!
```
In OpenFOAM's `wallHeatFlux` function object, heat entering the fluid domain is signed positive, while heat leaving across the isothermal cold boundary ($T_{\\text{fluid}} > T_{\\text{cold}} = 301.45\\text{ K}$) produces a negative patch integral ($\\dot{Q}_{\\text{cold, raw}} < 0$). Subtracting a negative cumulative energy (`ein - (-eout) = ein + eout`) double-counted outward cold-wall heat loss as an artificial heat injection into the domain. At 17 min, cumulative cold-wall loss is $246.61\\text{ J}$; adding rather than subtracting this loss inflated apparent net input to $8650.34\\text{ J}$, producing the spurious $-507.28\\text{ J}$ ($-5.86\\%$) deficit against domain stored energy ($8143.07\\text{ J}$).

When the physically correct definition is enforced ($E_{\\text{net, in}} = E_{\\text{hot}} - E_{\\text{cold, out}} = E_{\\text{hot}} + E_{\\text{cold, raw}}$):
- At $t = 2\\text{ min}$: $E_{\\text{net, in}} = 1887.30 - 6.99 = 1880.31\\text{ J}$, $E_{\\text{stored}} = 1878.37\\text{ J} \\implies \\Delta E = \\mathbf{-1.94\\text{ J}}$ ($\\mathbf{-0.10\\%}$).
- At $t = 6\\text{ min}$: $E_{\\text{net, in}} = 3729.13 - 48.24 = 3680.89\\text{ J}$, $E_{\\text{stored}} = 3676.92\\text{ J} \\implies \\Delta E = \\mathbf{-3.97\\text{ J}}$ ($\\mathbf{-0.11\\%}$).
- At $t = 17\\text{ min}$: $E_{\\text{net, in}} = 8395.60 - 246.61 = 8148.99\\text{ J}$, $E_{\\text{stored}} = 8143.07\\text{ J} \\implies \\Delta E = \\mathbf{-5.92\\text{ J}}$ ($\\mathbf{-0.073\\%}$).

The corrected script is committed in `validation/gallium2D/analysis/energy_balance/exact_energy_balance_audit.py`.

---

### 2.3 Mathematical Definition of the Step-Resolved Energy Audit
The authoritative energy audit evaluates boundary fluxes at every time step ($\\Delta t = 0.1\\text{ s}$, 10,200 steps to 17 min) against domain-integrated sensible and latent enthalpy:
$$E_{\\text{net, in}}(t_n) = E_{\\text{net}, 0 \\to t_1} + \\sum_{k=2}^n \\frac{1}{2}\\left[\\dot{Q}_{\\text{net}}(t_{k-1}) + \\dot{Q}_{\\text{net}}(t_k)\\right]\\Delta t_k$$
where $\\dot{Q}_{\\text{net}}(t) = \\dot{Q}_{\\text{hot}}(t) - \\dot{Q}_{\\text{cold, out}}(t)$. For the unrecorded initial startup step ($0 \\le t \\le t_1 = 0.1\\text{ s}$), gross hot-wall input is evaluated as a rectangular step:
$$E_{\\text{hot}, 0 \\to 0.1} = \\dot{Q}_{\\text{hot}}(t_1)\\Delta t_1$$
yielding $7.173\\text{ J}$ for baseline EP ($\\dot{Q}_{\\text{hot}} = 71.73\\text{ W}$) and $8.123\\text{ J}$ for EHCDarcy ($\\dot{Q}_{\\text{hot}} = 81.23\\text{ W}$). At $t_1 = 0.1\\text{ s}$, the hot-wall thermal disturbance has not produced a meaningful cold-wall response on this timescale, and the computed cold-wall heat flux is negligible at $t = 0.1\\text{ s}$ ($\\dot{Q}_{\\text{cold, out}}(0.1\\text{ s}) \\approx 4.5 \\times 10^{-11}\\text{ W}$, cumulative loss $< 10^{-11}\\text{ J}$), so $E_{\\text{net}, 0 \\to 0.1}$ matches gross hot-wall startup delivery to within $< 10^{-11}\\text{ J}$.

Stored enthalpy relative to initial conditions ($T_{\\text{init}} = 301.45\\text{ K}$, $f_L = 0$) evaluates to:
$$E_{\\text{sens}}(t) = \\int_\\Omega \\rho C_p (T - T_{\\text{init}})\\,dV, \\quad E_{\\text{lat}}(t) = \\int_\\Omega \\rho L f_L\\,dV, \\quad E_{\\text{stored}}(t) = E_{\\text{sens}}(t) + E_{\\text{lat}}(t)$$
Energy conservation closure error is defined as:
$$\\Delta E(t) = E_{\\text{stored}}(t) - E_{\\text{net, in}}(t), \\quad \\epsilon_{\\text{energy}}(t) = \\frac{\\Delta E(t)}{E_{\\text{net, in}}(t)} \\times 100\\%$$

---

### 2.4 Time-Resolved Authoritative Energy Audit History

#### Complete Transient History: `EHCDarcy` Formulation ($42 \\times 32$ Grid, $n_{\\text{outer}}=8$, $\\Delta t = 0.1\\text{ s}$)
| Time | $\\dot{Q}_{\\text{hot}}$ (W) | $\\dot{Q}_{\\text{cold}}$ (W) | $E_{\\text{hot}}$ (J) | $E_{\\text{cold, out}}$ (J) | $E_{\\text{net, in}}$ (J) | $E_{\\text{sens}}$ (J) | $E_{\\text{lat}}$ (J) | $E_{\\text{stored}}$ (J) | Discrepancy $\\Delta E$ (J) | Closure Error (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1 min** | $11.630$ | $0.0578$ | $1294.0$ | $0.98$ | $1293.0$ | $59.1$ | $1233.3$ | $1292.4$ | $-0.65\\text{ J}$ | **$-0.05\\%$** |
| **2 min** | $8.802$ | $0.1318$ | $1887.3$ | $6.99$ | $1880.3$ | $81.4$ | $1797.0$ | $1878.4$ | $-1.94\\text{ J}$ | **$-0.10\\%$** |
| **4 min** | $7.503$ | $0.1748$ | $2848.7$ | $26.00$ | $2822.7$ | $107.3$ | $2712.3$ | $2819.5$ | $-3.12\\text{ J}$ | **$-0.11\\%$** |
| **6 min** | $7.215$ | $0.1959$ | $3729.1$ | $48.24$ | $3680.9$ | $129.1$ | $3547.9$ | $3676.9$ | $-3.95\\text{ J}$ | **$-0.11\\%$** |
| **8 min** | $7.130$ | $0.2210$ | $4589.4$ | $73.20$ | $4516.2$ | $149.5$ | $4362.2$ | $4511.7$ | $-4.49\\text{ J}$ | **$-0.10\\%$** |
| **10 min** | $7.057$ | $0.2536$ | $5442.7$ | $101.59$ | $5341.2$ | $169.5$ | $5166.7$ | $5336.2$ | $-5.01\\text{ J}$ | **$-0.09\\%$** |
| **12 min** | $7.095$ | $0.2957$ | $6289.7$ | $134.42$ | $6155.3$ | $187.9$ | $5962.1$ | $6150.0$ | $-5.36\\text{ J}$ | **$-0.09\\%$** |
| **14 min** | $7.033$ | $0.3514$ | $7133.7$ | $173.06$ | $6960.6$ | $206.7$ | $6748.3$ | $6955.0$ | $-5.65\\text{ J}$ | **$-0.08\\%$** |
| **16 min** | $7.000$ | $0.4281$ | $7974.8$ | $219.50$ | $7755.3$ | $224.3$ | $7525.2$ | $7749.5$ | $-5.80\\text{ J}$ | **$-0.07\\%$** |
| **17 min** | $6.992$ | $0.4782$ | $8395.60$ | $246.61$ | $8148.99$ | $232.71$ | $7910.36$ | $8143.07$ | $-5.92\\text{ J}$ | **$-0.073\\%$** |

#### Complete Transient History: `enthalpyPorosity` (Default: $\lambda = 1.0, n_{\\text{outer}}=8$, $\\Delta t = 0.1\\text{ s}$)
| Time | $\\dot{Q}_{\\text{hot}}$ (W) | $\\dot{Q}_{\\text{cold}}$ (W) | $E_{\\text{hot}}$ (J) | $E_{\\text{cold, out}}$ (J) | $E_{\\text{net, in}}$ (J) | $E_{\\text{sens}}$ (J) | $E_{\\text{lat}}$ (J) | $E_{\\text{stored}}$ (J) | Discrepancy $\\Delta E$ (J) | Closure Error (%) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1 min** | $11.668$ | $0.0573$ | $1304.88$ | $0.97$ | $1303.91$ | $59.35$ | $1236.91$ | $1296.26$ | $-7.65\\text{ J}$ | **$-0.59\\%$** |
| **2 min** | $8.824$ | $0.1317$ | $1898.22$ | $7.01$ | $1891.21$ | $81.64$ | $1804.69$ | $1886.33$ | $-4.88\\text{ J}$ | **$-0.26\\%$** |
| **4 min** | $7.509$ | $0.1755$ | $2860.25$ | $26.13$ | $2834.11$ | $107.66$ | $2725.59$ | $2833.25$ | $-0.86\\text{ J}$ | **$-0.03\\%$** |
| **6 min** | $7.218$ | $0.1969$ | $3741.72$ | $48.49$ | $3693.23$ | $129.40$ | $3564.23$ | $3693.64$ | $+0.41\\text{ J}$ | **$+0.01\\%$** |
| **8 min** | $7.132$ | $0.2222$ | $4603.00$ | $73.56$ | $4529.44$ | $149.92$ | $4379.05$ | $4528.96$ | $-0.48\\text{ J}$ | **$-0.01\\%$** |
| **10 min** | $7.059$ | $0.2549$ | $5457.37$ | $102.05$ | $5355.32$ | $169.75$ | $5183.13$ | $5352.87$ | $-2.45\\text{ J}$ | **$-0.05\\%$** |
| **12 min** | $7.094$ | $0.2971$ | $6305.90$ | $134.96$ | $6170.94$ | $188.36$ | $5977.75$ | $6166.11$ | $-4.83\\text{ J}$ | **$-0.08\\%$** |
| **14 min** | $7.032$ | $0.3528$ | $7151.05$ | $173.65$ | $6977.40$ | $206.94$ | $6762.95$ | $6969.89$ | $-7.50\\text{ J}$ | **$-0.11\\%$** |
| **16 min** | $6.999$ | $0.4294$ | $7993.62$ | $220.10$ | $7773.52$ | $224.44$ | $7538.76$ | $7763.21$ | $-10.31\\text{ J}$ | **$-0.13\\%$** |
| **17 min** | $6.991$ | $0.4795$ | $8414.69$ | $247.19$ | $8167.51$ | $232.79$ | $7922.90$ | $8155.68$ | $-11.82\\text{ J}$ | **$-0.14\\%$** |

*Latent Enthalpy Reconciliation at $t = 17\\text{ min}$*:  
For EHCDarcy, $f_L = 0.5736135 \\implies E_{\\text{lat}} = m L f_L = 13790.40 \\times 0.5736135 = \\mathbf{7910.36\\text{ J}}$. Direct cell summation yields $\\sum \\rho V_{\\text{cell}} L \\alpha = 7910.36\\text{ J}$.  
For EP ($\\lambda = 1.0, n_{\\text{outer}}=8$), $f_L = 0.5745226 \\implies E_{\\text{lat}} = 13790.40 \\times 0.5745226 = \\mathbf{7922.90\\text{ J}}$. Direct cell summation yields $7922.90\\text{ J}$. Both formulations reconcile latent enthalpy with liquid fraction to $0.01\\text{ J}$.

---

### 2.5 Systematic Study of Newton Relaxation ($\\lambda$) and Outer Iterations in Enthalpy-Porosity
In OpenFOAM's `enthalpyPorosityPhaseChangeModel.C`, the liquid fraction is updated after the sensible enthalpy matrix solve via Newton projection:
$$f_{\\text{new}} = \\min\\left(\\max\\left(f + \\lambda \\frac{T - T_f}{L/C_p + \\Delta T_{\\text{window}}}, 0\\right), 1\\right)$$
where $\\lambda$ is the Newton relaxation factor (code default: `1.0`, the exact Newton step for a linear melting window). In the baseline tutorial case dictionary, this parameter was set to `lambda 0.2`, damping each outer iteration correction to only $20\\%$ of the full Newton step.

To decouple physical modeling capabilities from under-relaxation artifacts, a systematic sweep of $\\lambda \\in \\{1.0, 0.5, 0.2\\}$ across PIMPLE outer correctors $n_{\\text{outer}} \\in \\{2, 3, 5, 8\\}$ was executed to $t = 120\\text{ s}$ ($2\\text{ min}$) using the committed runner `validation/gallium2D/analysis/energy_balance/run_ep_relaxation_study.py`:

| Model Formulation | $\\lambda$ | $n_{\\text{outer}}$ | $f_L(2\\text{ min})$ | Net Input $E_{\\text{net}}$ (J) | Stored $E_{\\text{stored}}$ (J) | Discrepancy $\\Delta E$ (J) | Closure Error (%) | Clock Time (s) |
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
   Under tutorial settings ($\\lambda = 0.2$), error decayed geometrically ($88\\% \\to 42\\% \\to 17.5\\% \\to 7.2\\% \\to 2.7\\%$) because each pass applied only $20\\%$ of the required Newton correction. Setting $\\lambda = 1.0$ (exact Newton step) allows Enthalpy-Porosity to reach **$-0.26\\%$** closure at $n_{\\text{outer}} = 8$, easily passing the single $\\le 2.0\\%$ target and matching EHCDarcy's $-0.10\\%$ closure.
2. **Enthalpy-Porosity Still Requires $\\ge 6\\text{--}8$ Correctors**:  
   At small corrector counts ($n_{\\text{outer}} \\le 5$), Enthalpy-Porosity fails the $\\le 2.0\\%$ target across all tested relaxation factors ($+56.5\\%$ at $n=2$, $+11.1\\%$ at $n=3$, $-7.49\\%$ at $n=5$ for $\\lambda = 1.0$). Because latent heat is introduced via an explicit RHS source term ($S_h = L \\frac{\\partial(\\rho f)}{\\partial t}$) with $f$ projected after the matrix solve, updating $f$ changes the thermal field for subsequent iterations; at low corrector counts, iteration lag during rapid conduction surges prevents discrete energy closure.
3. **Formulation Trade-off**:  
   Enthalpy-Porosity is **not** a cheaper alternative at low corrector counts because it requires $n_{\\text{outer}} \\ge 8$ to conserve energy. However, when run with code-default $\\lambda = 1.0$ and $n_{\\text{outer}} = 8$, Enthalpy-Porosity achieves full code-to-code agreement with `EHCDarcy`:
   - Liquid fraction at 17 min: $f_L = 0.5745$ (EP) vs. $0.5736$ (EHCDarcy), a difference of only $0.09$ percentage points.
   - Front position at mid-height ($y/H=0.5$): $56.91\\text{ mm}$ (EP) vs. $57.06\\text{ mm}$ (EHCDarcy), agreeing within $0.15\\text{ mm}$.
   - Energy closure at 17 min: $-0.14\\%$ (EP) vs. $-0.073\\%$ (EHCDarcy).

---

"""

text = text[:sec2_start] + new_sec2 + text[sec3_start:]

# -----------------------------------------------------------------------------
# 2. Section 5: Clarify that ASME GCI Grid Study Covers EHCDarcy Only
# -----------------------------------------------------------------------------
text = text.replace(
    "## 5. Spatial Grid Sensitivity & ASME / Roache GCI Analysis (to 17 min)",
    "## 5. Spatial Grid Sensitivity & ASME / Roache GCI Analysis (EHCDarcy Formulation, to 17 min)"
)
text = text.replace(
    "To establish spatial discretization uncertainty for the gallium benchmark,",
    "To establish spatial discretization uncertainty for the `EHCDarcy` benchmark formulation (note: no formal grid study was performed for `enthalpyPorosity`),"
)

# -----------------------------------------------------------------------------
# 3. Section 7.3: Ye & Arıcı Correlation Scripting and Scaling Note
# -----------------------------------------------------------------------------
old_sec73_pattern = re.compile(r"### 7\.3 Consistency Check with Ye & Arıcı \(2023\) Numerical Correlation.*?(?=---|\Z)", re.DOTALL)
new_sec73 = """### 7.3 Consistency Check with Ye & Arıcı (2023) Numerical Correlation
Ye & Arıcı (2023, Section 3.2.2, Eq. 6) conducted an independent 2D verification study on natural convection melting of gallium, correlating their grid-converged numerical liquid fractions as:
$$\\phi_n = 2.708\\,\\tau^{0.712}\\,\\mathrm{Ra}^{0.0243}\\,\\mathrm{Ar}^{-0.0515}$$
where:
- $\\phi_n$ represents the numerical global liquid volume fraction ($f_L$).
- $\\tau = \\mathrm{Fo} \\cdot \\mathrm{Ste}$, with $\\mathrm{Fo} = \\frac{\\alpha t}{W^2} = 1.7761$ and $\\mathrm{Ste} = \\frac{C_p (T_h - T_m)}{L} = 0.03912$ (at $t = 1020\\text{ s}$).
- $\\mathrm{Ra} = \\frac{g \\beta (T_h - T_m) H^3}{\\nu \\alpha} = 6.063 \\times 10^5$.
- $\\mathrm{Ar} = \\frac{H}{W} = 0.7143$.

Evaluating Eq. (6) via the committed verification script (`validation/gallium2D/analysis/evaluate_ye_arici_correlation.py`) at $t = 17\\text{ min}$ ($1020\\text{ s}$) yields:
$$\\phi_n^{\\text{Ye \\& Arıcı (2023)}} = \\mathbf{0.5703}$$
- `phaseChangeFoam` (EHCDarcy, $42 \\times 32$): **$f_L = 0.5736$** (Discrepancy: **$+0.0033$, or $0.33$ percentage points**).
- `phaseChangeFoam` (EP, $\\lambda = 1.0, n=8$): **$f_L = 0.5745$** (Discrepancy: **$+0.0042$, or $0.42$ percentage points**).
- `phaseChangeFoam` (EHCDarcy, $63 \\times 48$): **$f_L = 0.5725$** (Discrepancy: **$+0.0022$, or $0.22$ percentage points**).
- `phaseChangeFoam` (EHCDarcy, $32 \\times 24$): **$f_L = 0.5682$** (Discrepancy: **$-0.0021$, or $0.21$ percentage points**).

**Notes on Equation Provenance & Physical Scaling**:
1. *Formula Clarification*: An earlier draft cited a different correlation from another Ye & Arıcı paper ($0.03816\\,\\mathrm{Fo}^{0.4907}\\mathrm{Ra}^{0.3164}\\mathrm{Ste}^{0.8407}$) fitted to mean melt-layer thickness, which evaluates to $\\sim 0.23\\text{--}0.38$ for these dimensionless numbers. The authoritative benchmark liquid-fraction correlation is Eq. (6) above, which is fully documented and evaluated by the committed script.
2. *Dominance of Conduction/Stefan Scaling*: The exponent on the Rayleigh number is only $\\mathrm{Ra}^{0.0243}$. Consequently, varying $\\mathrm{Ra}$ by an entire order of magnitude ($10\\times$) modifies the predicted liquid fraction by only $10^{0.0243} = 1.0575$ ($+5.75\\%$). Agreement with this correlation primarily verifies the transient conduction and Stefan scaling ($\\tau^{0.712}$), rather than fine convective boundary layer dynamics.
3. *Interpretation*: Agreement with a numerical correlation fitted to independent CFD simulations provides valuable code-to-code verification, but is strictly distinguished from experimental validation against physical measurements.

"""

text = old_sec73_pattern.sub(lambda m: new_sec73, text)

# -----------------------------------------------------------------------------
# 4. Section 8: Summary Conclusions & Verification Scope Assessment
# -----------------------------------------------------------------------------
old_sec8_pattern = re.compile(r"## 8\. Summary Conclusions & Verification Scope Assessment.*", re.DOTALL)
new_sec8 = r"""## 8. Summary Conclusions & Verification Scope Assessment

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
"""

text = old_sec8_pattern.sub(lambda m: new_sec8, text)

with open(report_path, "w") as f:
    f.write(text)

print("ANALYSIS_REPORT.md updated successfully!")

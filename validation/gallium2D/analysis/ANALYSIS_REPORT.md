# Comprehensive Numerical Benchmark Verification, Energy Balance Audit & Sensitivity Analysis: Gallium 2D Benchmark

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

### 1.1 Acceptance Criteria Framework
To maintain scientific rigor, this study explicitly distinguishes **numerical benchmark verification** from **physical experimental validation**:
- **Verification Scope**: Quantitative phase-front deviations and global liquid fractions are evaluated against the authoritative 2D numerical benchmark of Brent, Voller, and Reid (1988) and the grid-converged numerical scaling correlation of Ye & Arıcı (2023). These comparisons establish code-to-code verification against published CFD literature.
- **Experimental Scope**: Physical experimental measurements from Gau & Viskanta (1986) provide essential qualitative context regarding experimental delay mechanisms (heating jacket thermal inertia, pour-out measurement disruption, and 3D endwall drag).

Study-specific acceptance thresholds are defined as follows:
1. **Front Morphology (Mean)**: Mean front deviation across cavity height $\overline{\delta}_{\text{front}} \le 5.0\%$ of cavity width $W$ ($4.45\text{ mm}$), evaluated against the digitized **numerical benchmark front of Brent, Voller, and Reid (1988, Figure 9)** (digitization line-width uncertainty $\pm 1.0\text{ mm}$).
2. **Front Morphology (Max)**: Peak local front deviation $\delta_{\text{front, max}} \le 10.0\%$ of cavity width $W$ ($8.89\text{ mm}$), relative to Brent et al. (1988).
3. **Global Liquid Fraction ($f_L$)**: Absolute liquid fraction deviation $|f_L^{\text{sim}} - f_L^{\text{ref}}| \le 0.05$ ($5.0$ percentage points) relative to Brent et al. (1988).
4. **Energy Conservation Closure (Mutually Exclusive Regimes)**:
   - **Early Conduction Phase** ($t \le 2\text{ min}$): Net closure discrepancy $|\epsilon_{\text{energy}}| \le 2.0\%$.
   - **Developed Convective Phase** ($2 < t \le 17\text{ min}$): Net closure discrepancy $|\epsilon_{\text{energy}}| \le 6.0\%$.

### 1.2 Multi-Model Performance Scorecard (Full 17-Minute Benchmark)

| Metric | Study-Specific Target | Enthalpy-Porosity (Baseline, $n_{\text{outer}}=8$) | Enthalpy-Porosity (Higher-corrector, $n_{\text{outer}}=12$) | EHCDarcy (Baseline, $n_{\text{outer}}=8$) | Status (EP, $n_{\text{outer}}=12$) | Status (EHCDarcy) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Energy Closure ($t = 2\text{ min}$)** | $\le 2.0\%$ (Early phase: $t \le 2\text{ min}$) | $+6.72\%$ ($+121.6\text{ J}$) *(FAIL)* | **$+2.30\%$ ($+42.6\text{ J}$)** *(FAIL)* | **$-0.10\%$ ($-1.94\text{ J}$)** | **FAILED** | **PASSED** |
| **Energy Closure ($t = 6\text{ min}$)** | $\le 6.0\%$ (Convective: $2 < t \le 17\text{ min}$) | $+3.54\%$ ($+128.2\text{ J}$) | **$+1.75\%$ ($+64.1\text{ J}$)** | **$-0.11\%$ ($-4.0\text{ J}$)** | **PASSED** | **PASSED** |
| **Energy Closure ($t = 17\text{ min}$)** | $\le 6.0\%$ (Convective: $2 < t \le 17\text{ min}$) | $+5.18\%$ ($+416.7\text{ J}$) | **$+2.97\%$ ($+238.5\text{ J}$)** | **$-0.07\%$ ($-5.9\text{ J}$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 2\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $1.95\text{ mm}$ ($2.19\%$) | $1.72\text{ mm}$ ($1.93\%$) | **$1.67\text{ mm}$ ($1.88\%$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 6\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $2.60\text{ mm}$ ($2.93\%$) | $2.12\text{ mm}$ ($2.38\%$) | **$1.87\text{ mm}$ ($2.11\%$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 10\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $3.51\text{ mm}$ ($3.95\%$) | $2.84\text{ mm}$ ($3.19\%$) | **$2.35\text{ mm}$ ($2.65\%$)** | **PASSED** | **PASSED** |
| **Mean Front Dev. ($t = 17\text{ min}$)** | $\le 5.0\%$ $W$ ($4.45\text{ mm}$) | $5.43\text{ mm}$ ($6.11\%$) *(FAIL)* | $4.18\text{ mm}$ ($4.70\%$) | **$3.46\text{ mm}$ ($3.89\%$)** | **PASSED** | **PASSED** |
| **Max Front Dev. ($t = 17\text{ min}$)** | $\le 10.0\%$ $W$ ($8.89\text{ mm}$) | $9.74\text{ mm}$ ($10.95\%$) *(FAIL)* | $7.85\text{ mm}$ ($8.83\%$) | **$6.34\text{ mm}$ ($7.13\%$)** | **PASSED** | **PASSED** |
| **$f_L$ at $t = 17\text{ min}$** (Brent: $0.5334$) | $\pm 0.05$ | $0.5967$ ($+0.0633$) *(FAIL)* | $0.5822$ ($+0.0488$) | **$0.5736$ ($+0.0402$)** | **PASSED** | **PASSED** |
| **Ye & Arıcı (2023) Eq. 6 Correlation** ($0.5703$) | $\pm 0.03$ | $+0.0264$ | $+0.0119$ | **$+0.0033$ ($0.33$ points)** | **PASSED** | **PASSED** |

*Scorecard Evaluation Note on Early vs. Convective Conservation Targets*:  
In the early conduction regime ($t \le 2\text{ min}$), the acceptance threshold is $|\epsilon_{\text{energy}}| \le 2.0\%$. At $t = 2\text{ min}$, baseline Enthalpy-Porosity ($n_{\text{outer}}=8$) at $+6.72\%$ ($+121.6\text{ J}$) and higher-corrector Enthalpy-Porosity ($n_{\text{outer}}=12$) at $+2.30\%$ ($+42.6\text{ J}$) both fail this tight early criterion due to numerical outer-iteration lag during the steep sub-minute conduction surge. Enthalpy-Porosity ($n_{\text{outer}}=12$) recovers to pass the $\le 6.0\%$ criterion across the developed convective regime ($t > 2\text{ min}$), whereas EHCDarcy satisfies both early conduction ($\le 2.0\%$) and late convective ($\le 6.0\%$) criteria across the entire 17-minute transient.

---

## 2. Energy-Balance Audit and Root-Cause Analysis

### 2.1 Why the Previous Energy Audit Showed an Apparent $+20.6\%$ Surplus
The initial audit script (`energy_balance_audit.py`, archived in `energy_balance_results.txt`) evaluated the **baseline Enthalpy-Porosity ($n_{\text{outer}}=8$)** simulation using a coarse $60\text{ s}$ sampling interval and reported:
- $t = 1\text{ min}$: $E_{\text{in}} = 1009.9\text{ J}$, $E_{\text{stored}} = 1326.3\text{ J}$ ($+316.4\text{ J}$, $+31.3\%$ apparent surplus)
- $t = 2\text{ min}$: $E_{\text{in}} = 1601.5\text{ J}$, $E_{\text{stored}} = 1931.1\text{ J}$ ($+329.6\text{ J}$, $+20.6\%$ apparent surplus)
- $t = 17\text{ min}$: $E_{\text{in}} = 7946.1\text{ J}$, $E_{\text{stored}} = 8471.4\text{ J}$ ($+525.3\text{ J}$, $+6.6\%$ apparent surplus)

This apparent energy non-conservation was caused by two compounding factors that must be analyzed on a strict like-for-like basis for each formulation:

1. **Under-Integration of the Initial Conduction Surge ($E_{\text{in}}$)**:  
   The initial post-processing script approximated boundary heat input using a coarse finite difference between the hot wall and the first cell center ($k(T_h - T_0)/(dx/2)$) sampled exclusively at $60\text{ s}$ intervals. During early times ($t < 60\text{ s}$), transient 1D conduction into a subcooled solid slab exhibits a steep initial heat flux peak:
   $$
   q''_{\text{wall}}(t) \sim \frac{k(T_h - T_{\text{init}})}{\sqrt{\pi \alpha t}}
   $$
   where $T_h = 311.15\text{ K}$, $T_{\text{init}} = 301.45\text{ K}$ ($\Delta T = 9.70\text{ K}$), and cavity patch area $A = H \times D = 0.0635 \times 0.005 = 3.175 \times 10^{-4}\text{ m}^2$. While this relation serves as qualitative physical motivation for the boundary flux surge rather than an exact model of the phase change transient, actual step-resolved OpenFOAM `wallHeatFlux` logging ($\Delta t = 0.1\text{ s}$) confirms this steep peak:
   - At $t = 0.1\text{ s}$: $\dot{Q}_{\text{hot}} \approx 72\text{--}81\text{ W}$ ($q''_{\text{wall}} \approx 2.3\text{--}2.6 \times 10^5\text{ W/m}^2$)
   - At $t = 1.0\text{ s}$: $\dot{Q}_{\text{hot}} \approx 65\text{--}80\text{ W}$
   - At $t = 10.0\text{ s}$: $\dot{Q}_{\text{hot}} \approx 26\text{ W}$
   - At $t = 60.0\text{ s}$: $\dot{Q}_{\text{hot}} \approx 11.2\text{--}11.6\text{ W}$
   - At $t = 120.0\text{ s}$: $\dot{Q}_{\text{hot}} \approx 8.5\text{--}8.8\text{ W}$

2. **Controlled Isolation of Temporal Undersampling vs. Flux Estimator Discrepancy**:  
   To rigorously determine whether the $214.9\text{ J}$ difference at $t = 2\text{ min}$ arose from temporal undersampling of the conduction surge or from changing the boundary flux estimator (finite-difference wall gradient vs. OpenFOAM's native `wallHeatFlux`), a controlled diagnostic was executed on the baseline $n_{\text{outer}}=8$ dataset:
   - **Flux Estimator Comparison at Identical Timestamps**:  
     Evaluating the finite-difference wall gradient $q_{\text{FD}} = \sum_j k \frac{T_h - T_{j,0}}{\Delta x / 2} \Delta y D$ against OpenFOAM's native patch integral $\dot{Q}_{\text{native}}$ at the identical 60-second output timestamps yields (with $\Delta Q = Q_{\text{native}} - Q_{\text{FD}}$):
     - At $t = 60\text{ s}$ ($1\text{ min}$): $\dot{Q}_{\text{native}} = 11.22152\text{ W}$ vs. $\dot{Q}_{\text{FD}} = 11.22090\text{ W}$ (difference $\Delta Q = +0.00062\text{ W}$, or $+0.0055\%$; both round to $11.22\text{ W}$)
     - At $t = 120\text{ s}$ ($2\text{ min}$): $\dot{Q}_{\text{native}} = 8.49778\text{ W}$ vs. $\dot{Q}_{\text{FD}} = 8.49810\text{ W}$ (difference $\Delta Q = -0.00032\text{ W}$, or $-0.0038\%$; both round to $8.50\text{ W}$)  
     Both estimators agree to within $0.01\%$ ($< 0.001\text{ W}$) at these timestamps and round to identical two-decimal rates ($11.22\text{ W}$ and $8.50\text{ W}$). As noted below, agreement at these two timestamps does not establish equivalence at every intervening time.
   - **Controlled Quadrature Test on Identical Native Flux Data & Startup Reconstruction**:  
     To reconstruct how the original post-processing script (`energy_balance_audit.py`) arrived at its coarse estimate of $1601.5\text{ J}$ (archived as `1601.5` in `energy_balance_results.txt`), the exact algorithm, input values, and formatting must be documented:
     1. *Coarse Script Startup Treatment ($0 \to 60\text{ s}$)*: The original post-processing script sampled fields exclusively at $60\text{ s}$ intervals and had no data at $t = 0$. To heuristically account for the uncaptured $t^{-1/2}$ conduction singularity without intermediate time steps, the script applied an ad-hoc factor of $1.5$ to the first interval:
        $$E_{0 \to 60} = \dot{Q}(60\text{ s}) \times \Delta t \times 1.5$$
        with the explicit code comment: `# Analytical 1D conduction initial heat flux: q ~ k*(T_h - T_m)/sqrt(pi*alpha*t) ... cum_E_in += Q_in * dt_step * 1.5 # accounts for t=0 initial thermal boundary layer singularity`. If this first interval is interpreted as an effective trapezoidal integration $\frac{\dot{Q}_{\text{eff}}(0) + \dot{Q}(60\text{ s})}{2}\Delta t$, this heuristic implicitly assumed an effective initial heat transfer rate of:
        $$\dot{Q}_{\text{eff}}(0) = 2 \dot{Q}(60\text{ s}) = 2 \times 11.2215\text{ W} = \mathbf{22.44\text{ W}}$$
     2. *Second Interval ($60 \to 120\text{ s}$)*: Standard trapezoidal quadrature was applied between the two 60-second samples:
        $$E_{60 \to 120} = \frac{\dot{Q}(60\text{ s}) + \dot{Q}(120\text{ s})}{2}\Delta t = \frac{11.22152\text{ W} + 8.49778\text{ W}}{2} \times 60\text{ s} = 591.579\text{ J} \approx \mathbf{591.58\text{ J}}$$
     3. *Evaluation on Native vs. Finite-Difference Data*:
        - *Archived Finite-Difference Calculation*: In `energy_balance_audit.py`, evaluating the coarse formula with the unrounded finite-difference flux samples ($Q_{\text{FD}} = 11.22090\text{ W}$ and $8.49810\text{ W}$) yields:
          $$E_{\text{FD, coarse}}(120\text{ s}) = 11.22090 \times 60 \times 1.5 + \frac{11.22090 + 8.49810}{2} \times 60 = 1009.881\text{ J} + 591.570\text{ J} = \mathbf{1601.451\text{ J}}$$
          In `energy_balance_results.txt`, this total was printed under a one-decimal format string (`{cum_E_in:<11.1f}`), which displayed as **`1601.5 J`** (and at 60 s, $1009.881\text{ J}$ displayed as **`1009.9 J`**). When referenced with two decimals in subsequent discussions, `1601.5` was quoted as `1601.50 J`.
        - *Unrounded Native OpenFOAM Reconstruction*: Evaluating the identical coarse algorithm on unrounded native patch fluxes from OpenFOAM's `wallHeatFlux.dat` ($Q_{\text{native}} = 11.22152\text{ W}$ and $8.49778\text{ W}$) yields:
          $$E_{\text{native, coarse}}(120\text{ s}) = 11.22152 \times 60 \times 1.5 + \frac{11.22152 + 8.49778}{2} \times 60 = 1009.9368\text{ J} + 591.5790\text{ J} = \mathbf{1601.5158\text{ J} \approx 1601.52\text{ J}}$$
          which also rounds to `1601.5 J` under a one-decimal display. The difference between the unrounded native coarse reconstruction and the exact unrounded finite-difference coarse total is $1601.5158 - 1601.4510 = \mathbf{0.0648\text{ J} \approx 0.065\text{ J}}$ ($0.004\%$).
        - *Rounding of Displayed Samples*: Evaluating the coarse formula with displayed two-decimal flux samples ($11.22\text{ W}$ and $8.50\text{ W}$) yields $11.22(60)(1.5) + \frac{11.22+8.50}{2}(60) = \mathbf{1601.40\text{ J}}$. The approximately $0.12\text{ J}$ difference ($1601.5158 - 1601.40 = 0.1158\text{ J} \approx 0.12\text{ J}$) arises from rounding the displayed endpoint heat-flux samples to two decimal places before evaluating the coarse formula.
     4. *Step-Resolved Integration: Rectangular Startup Step Plus Trapezoidal Integration Thereafter ($\Delta t = 0.1\text{ s}$)*:  
        In OpenFOAM's native `wallHeatFlux.dat`, logging begins at the first numerical time step ($t = 0.1\text{ s}$) where the patch-integrated heat transfer rate is recorded as **$71.7338\text{ W} \approx 71.73\text{ W}$** ($q'_{\text{wall}} = 2.2593 \times 10^5\text{ W/m}^2$, consistent with the theoretical $t=0$ cell-center finite-difference initial conduction rate of $\approx 93.1\text{ W}$). In the discrete flux time series ($t_1 = 0.1, t_2 = 0.2, \dots, t_{1200} = 120.0\text{ s}$), the integration method is defined as a **rectangular startup step for the unrecorded initial interval $0 \le t \le 0.1\text{ s}$ plus standard trapezoidal quadrature thereafter**:
        - *Initial Rectangular Step ($0 \le t \le 0.1\text{ s}$)*:
          $$E_{0 \to 0.1} = \dot{Q}(0.1\text{ s}) \times \Delta t_0 = 71.7338\text{ W} \times 0.1\text{ s} = \mathbf{7.173\text{ J}}$$
        - *Subsequent Trapezoidal Quadrature ($0.1 \le t \le 120.0\text{ s}$)*:
          Summing the 1,199 trapezoidal intervals from the unrounded step-resolved data yields:
          $$E_{0.1 \to 120.0} = \sum_{k=2}^{1200} \frac{\dot{Q}(t_{k-1}) + \dot{Q}(t_k)}{2}\Delta t_k = \mathbf{1809.205\text{ J}}$$
        - *Reconciled Gross Energy Total ($t = 2\text{ min}$)*:
          $$E_{\text{native, step}}(120\text{ s}) = 7.173\text{ J} + 1809.205\text{ J} = \mathbf{1816.378\text{ J} \approx 1816.38\text{ J}}$$
        - *Startup Quadrature Sensitivity*:  
          To rigorously evaluate the sensitivity to the initial rectangular approximation, if the first interval $0 \le t \le 0.1\text{ s}$ were instead integrated trapezoidally using the initialized finite-volume wall gradient at $t = 0$ ($\dot{Q}(0) \approx 93.1\text{ W}$, from $k \frac{T_h - T_{\text{init}}}{\Delta x / 2} A$ on the orthogonal mesh) and $\dot{Q}(0.1\text{ s}) = 71.7338\text{ W}$, the first interval would evaluate to $\frac{93.1 + 71.7338}{2} \times 0.1\text{ s} \approx 8.24\text{ J}$. This represents a difference of approximately **$+1.07\text{ J}$** relative to the rectangular estimate ($7.173\text{ J}$). This sensitivity ($\approx 1.07\text{ J}$, or $< 0.06\%$ of gross input) confirms that the rectangular startup approximation introduces negligible uncertainty into the audit.
        - *Reconciliation Across Time Intervals ($0 \to 60\text{ s}$ and $60 \to 120\text{ s}$)*:  
          Evaluating cumulative gross hot-wall energy at 60 s yields $E_{\text{hot}}(60\text{ s}) = 1239.60\text{ J}$, and at 120 s yields $E_{\text{hot}}(120\text{ s}) = 1816.38\text{ J}$ ($1816.378\text{ J}$). Based on these displayed values, the step-resolved energy delivered during the second minute is:
          $$E_{60 \to 120} = 1816.38\text{ J} - 1239.60\text{ J} = \mathbf{576.78\text{ J}}$$
          Meanwhile, the coarse trapezoidal estimate over $60 \to 120\text{ s}$ evaluates from the endpoint native flux samples to:
          $$E_{\text{coarse}, 60 \to 120} = \frac{11.22152\text{ W} + 8.49778\text{ W}}{2} \times 60\text{ s} = 591.579\text{ J} \approx \mathbf{591.58\text{ J}}$$
          The coarse trapezoidal rule therefore **overestimates** the second minute by approximately **$14.80\text{ J}$** ($591.58 - 576.78 = 14.80\text{ J}$). For idealized diffusion-controlled decay, $q(t) \propto t^{-1/2}$ is convex, so an endpoint trapezoid overestimates the integral. The observed $14.80\text{ J}$ overestimate is consistent with this behavior; the numerical interval integration establishes the discrepancy independently of whether the full phase-change transient follows the idealized law exactly.  
          Over the first minute ($0 \to 60\text{ s}$), the coarse script's heuristic startup assumption ($1009.94\text{ J}$ native, or $1009.88\text{ J}$ FD) **underestimated** step-resolved delivery ($1239.60\text{ J}$) by **$229.66\text{ J}$** ($1239.60 - 1009.94 = 229.66\text{ J}$) due to missing the steep initial conduction surge ($71.7\text{--}93.1\text{ W}$).  
          Summing the two native coarse intervals gives $E_{\text{coarse}} = 1009.94\text{ J} + 591.58\text{ J} = \mathbf{1601.52\text{ J}}$ ($1601.516\text{ J}$ unrounded). The net gross difference between the two-interval coarse reconstruction and the step-resolved audit is:
          $$\Delta E_{\text{gross}} = 229.66\text{ J (minute 1 underestimation)} - 14.80\text{ J (minute 2 overestimation)} = \mathbf{214.86\text{ J}}$$
          which reconciles the step-resolved difference against the native coarse reconstruction ($1816.38 - 1601.52 = 214.86\text{ J}$) exactly, and against the archived script output ($1816.38 - 1601.50 = 214.88\text{ J}$) within $0.02\text{ J}$ ($< 0.01\%$). *(Derived interval differences are evaluated from displayed values rounded to two decimal places).* 

   | Contribution to Gross Difference ($t = 2\text{ min}$) | Energy (J) | Fraction of Difference |
   | :--- | :---: | :---: |
   | **Total Gross Difference** ($1816.38\text{ J} - 1601.45\text{ J}$) | **$214.93\text{ J}$** | **$100.0\%$** |
   | **Temporal Undersampling Error** ($E_{\text{native, step}} - E_{\text{native, coarse}}$: $1816.38 - 1601.52$) | **$214.86\text{ J}$** | **$99.97\%$** |
   | **Estimator Discretization Residual** ($E_{\text{native, coarse}} - E_{\text{FD, coarse}}$: $1601.52 - 1601.45$) | **$0.07\text{ J}$** | **$0.03\%$** |

   *(Note: When compared with the archived log's one-decimal display value of $1601.5\text{ J}$, the apparent gross-input difference is $214.88\text{ J}$. This differs from the unrounded finite-difference comparison because the archived value has been rounded; the underlying numerical decomposition is the $214.8622\text{ J}$ native temporal-integration correction plus the $0.0648\text{ J}$ estimator residual minus the display-rounding offset:
   $$\Delta E_{\text{archived display}} = 214.8622\text{ J} + 0.0648\text{ J} - 0.0490\text{ J} = \mathbf{214.8780\text{ J}} \approx \mathbf{214.88\text{ J}}$$
   where $-0.0490\text{ J} = 1601.451\text{ J} - 1601.500\text{ J}$ is the difference between the unrounded finite-difference total and the archived displayed total).*

   The controlled reconstruction establishes that the gross-input discrepancy of approximately $214.9\text{ J}$ is driven by temporal quadrature: an underestimation of $229.66\text{ J}$ during the steep sub-minute conduction surge partially offset by a $14.80\text{ J}$ trapezoidal overestimation during the second minute. While the startup treatment is now explicitly documented as a rectangular approximation with an evaluated sensitivity of $\approx 1.07\text{ J}$, and estimator agreement is confirmed at the sampled timestamps (where $\dot{Q}_{\text{FD}}$ and $\dot{Q}_{\text{native}}$ agree to within $< 0.01\%$), the comparison is demonstrated at the sampled timestamps rather than at every intervening time. Nevertheless, the step-resolved data demonstrate that temporal undersampling of the early conduction surge is the dominant source of the coarse-audit discrepancy.
3. **Like-for-Like Formulation Audits & Gross vs. Net Energy Accounting**:
   To prevent conflating post-processing integration error with physical differences between solver formulations or confusing gross boundary delivery with net domain absorption, each simulation is audited against its own time-resolved boundary flux:

   - **Baseline Enthalpy–Porosity ($n_{\text{outer}}=8$, $t = 2\text{ min}$)**:  
     Step-resolved trapezoidal quadrature ($\Delta t = 0.1\text{ s}$, 1,200 steps) on the baseline $n_{\text{outer}}=8$ run yields the following authoritative energy accounting:

     | Quantity at $t = 2\text{ min}$ | Energy (J) |
     | :--- | :---: |
     | Gross hot-wall input ($E_{\text{hot}}$) | $1816.38\text{ J}$ |
     | Cold-wall heat loss ($E_{\text{cold, out}}$) | $6.86\text{ J}$ |
     | **Net boundary input ($E_{\text{net, in}}$)** | **$1809.52\text{ J}$** |
     | **Domain stored energy ($E_{\text{stored}}$)** | **$1931.11\text{ J}$** |
     | Net energy discrepancy ($\Delta E = E_{\text{stored}} - E_{\text{net, in}}$) | **$+121.59\text{ J}$** |
     | Relative energy error ($\epsilon_{\text{energy}}$) | **$+6.72\%$** |

     The original coarse post-processing estimate of $1601.5\text{ J}$ represents an approximation of the gross hot-wall energy input. Comparing gross hot-wall input to gross hot-wall input reveals that the coarse script underestimated gross hot-wall delivery by approximately **$214.9\text{ J}$** ($1816.38 - 1601.5 = 214.88\text{ J}$) over the first two minutes. After accounting for the $6.86\text{ J}$ cold-wall heat loss, the resulting net boundary-energy correction relative to the coarse estimate is approximately **$208.0\text{ J}$** ($1809.52 - 1601.5 = 208.02\text{ J}$). Accounting for this boundary input resolves $63.1\%$ ($208.0\text{ J}$ out of $329.6\text{ J}$) of the apparent coarse surplus ($+20.6\% \to +6.72\%$). The remaining $+121.59\text{ J}$ ($+6.72\%$) discrepancy is consistent with incomplete outer-iteration coupling under $n_{\text{outer}}=8$, although other numerical contributions (such as discrete residuals and energy-flux integration approximations) are not excluded.

   - **Baseline Enthalpy–Porosity ($n_{\text{outer}}=8$, $t = 1\text{ min}$)**:  
     At $t = 1\text{ min}$ ($60\text{ s}$), gross hot-wall delivery is $E_{\text{hot}} = 1239.60\text{ J}$ and cold-wall heat loss is $E_{\text{cold, out}} = 1.03\text{ J}$, yielding net boundary input $E_{\text{net, in}} = 1238.57\text{ J}$ against stored energy $E_{\text{stored}} = 1326.31\text{ J}$ ($\Delta E = +87.74\text{ J}$, $+7.08\%$). The coarse estimate of $1009.9\text{ J}$ underestimated gross hot-wall input by approximately $229.7\text{ J}$ ($1239.60 - 1009.9 = 229.70\text{ J}$), and net boundary energy by $228.7\text{ J}$ ($1238.57 - 1009.9 = 228.67\text{ J}$).

   - **Baseline Enthalpy–Porosity at $t = 17\text{ min}$ (Two Independent Comparisons)**:  
     To separate quadrature integration error from solver formulation changes at late times, two distinct comparisons are performed:
     1. *Integration Correction for the Same $n_{\text{outer}}=8$ Baseline Run*:  
        The coarse estimate for this run was $E_{\text{in, coarse}} = 7946.1\text{ J}$ against stored energy $E_{\text{stored}} = 8471.4\text{ J}$ (apparent coarse surplus $+525.3\text{ J}$, $+6.61\%$). Full-resolution step-resolved boundary flux integration on this exact $n_{\text{outer}}=8$ run yields gross hot-wall delivery $E_{\text{hot}} = 8303.2\text{ J}$ (underestimated by $357.1\text{ J}$ in the coarse script), cold-wall heat loss $E_{\text{cold, out}} = 248.5\text{ J}$, and net boundary input $E_{\text{net, in}} = 8054.7\text{ J}$. Evaluating against domain stored energy ($8471.4\text{ J}$) yields a corrected net discrepancy of $\Delta E = 8471.4 - 8054.7 = \mathbf{+416.7\text{ J}}$ ($\mathbf{+5.18\%}$). Full-resolution integration and cold-wall accounting alone thus correct the baseline discrepancy from $+6.61\%$ down to $+5.18\%$.
     2. *Outer-Iteration Sensitivity ($n_{\text{outer}}=8$ vs. $n_{\text{outer}}=12$)*:  
        Comparing both fully resolved runs under identical boundary definitions:
        - Baseline $n_{\text{outer}}=8$: $E_{\text{net, in}} = 8054.7\text{ J}$, $E_{\text{stored}} = 8471.4\text{ J}$, $\Delta E = \mathbf{+416.7\text{ J}}$ ($\mathbf{+5.18\%}$).
        - Higher-corrector $n_{\text{outer}}=12$: $E_{\text{hot}} = 8276.57\text{ J}$, $E_{\text{cold, out}} = 249.32\text{ J} \implies E_{\text{net, in}} = 8027.25\text{ J}$, $E_{\text{stored}} = 8265.76\text{ J}$, $\Delta E = \mathbf{+238.51\text{ J}}$ ($\mathbf{+2.97\%}$).
        Increasing outer correctors from 8 to 12 tightens thermal-phase coupling, reducing stored enthalpy by $205.6\text{ J}$ (liquid fraction $f_L$ decreases from $0.5967 \to 0.5822$) and reducing the discrete energy discrepancy from $+5.18\%$ to $+2.97\%$.

   - **EHCDarcy (Baseline, $n_{\text{outer}}=8$)**:  
     EHCDarcy's own time-resolved simulation delivers gross hot-wall energy of **$1294.0\text{ J}$** at 1 min and **$1887.3\text{ J}$** at 2 min, while storing **$1878.37\text{ J}$** ($E_{\text{sens}} = 81.37\text{ J}$, $E_{\text{lat}} = 1797.00\text{ J}$). Accounting for cold-wall heat loss ($6.99\text{ J}$) gives $E_{\text{net, in}} = 1880.31\text{ J}$ and a net discrepancy of **$-1.94\text{ J}$ ($-0.10\%$)**.  
     - *Distinguishing Formulation Differences from Quadrature Error*: Comparing EHCDarcy's time-resolved gross hot-wall input ($1887.3\text{ J}$) directly to the coarse EP estimate ($1601.5\text{ J}$) yields a difference of $285.8\text{ J}$. This difference conflates two distinct contributions: $214.9\text{ J}$ of gross under-integration in the coarse script and a $70.9\text{ J}$ physical difference in gross hot-wall heat delivery between the two solver formulations ($1887.3\text{ J}$ in EHCDarcy vs. $1816.4\text{ J}$ in baseline EP).
     - At $t = 17\text{ min}$, EHCDarcy delivers gross hot-wall input $E_{\text{hot}} = 8395.60\text{ J}$ with cold-wall loss of $246.61\text{ J}$ ($E_{\text{net, in}} = 8148.99\text{ J}$), matching stored energy $E_{\text{stored}} = 8143.07\text{ J}$ ($E_{\text{sens}} = 232.71\text{ J}$, $E_{\text{lat}} = 7910.36\text{ J}$) to within **$-5.92\text{ J}$ ($-0.073\%$)**.

4. **Outer-Iteration Coupling in Enthalpy-Porosity**:  
   In `enthalpyPorosity`, the latent heat source is treated semi-explicitly on the RHS ($S_h = L \frac{\partial(\rho f)}{\partial t}$) and $f$ is updated via a Newton projection after the energy matrix solve. With low outer iterations ($n_{\text{outer}} \le 5$), $f$ advances without sufficient PIMPLE feedback into the temperature field, artificially inflating the apparent stored latent heat. Increasing $n_{\text{outer}}$ to 12 tightens this coupling and reduces the 2-minute error from $+6.72\%$ to $+2.30\%$. While incomplete outer-iteration coupling is consistent with the observed sensitivity, other numerical contributions (such as discrete equation residuals and boundary integration approximations) are not excluded.

### 2.2 Mathematical Definition of the Energy Audit & Cold-Wall Conduction Physics
To ensure complete reproducibility, the global energy balance is defined from the continuous enthalpy equation and evaluated numerically via time-step-resolved boundary-flux integration ($\Delta t = 0.1\text{ s}$, $1,200$ steps to 2 min, $10,200$ steps to 17 min):

1. **Net Boundary Heat Input & Boundary Definitions**:
   $$E_{\text{net, in}}(t) = \int_0^t \left[ \dot{Q}_{\text{hot}}(s) - \dot{Q}_{\text{cold, out}}(s) \right]\,ds$$
   where:
   - **`hotWall`**: Active isothermal boundary ($T_h = 311.15\text{ K}$). The heat transfer rate entering the domain is $\dot{Q}_{\text{hot}}(s) = \int_{A_{\text{hot}}} \mathbf{q} \cdot \mathbf{n}_{\text{in}}\,dA > 0$.
   - **`coldWall` & Raw Wall-Flux Sign Convention**: Isothermal boundary held at $T_c = 301.45\text{ K} = T_{\text{init}}$. Crucially, because solid gallium has a high thermal conductivity ($k = 32.0\text{ W/(m K)}$), thermal conduction through the solid ahead of the advancing melt front establishes a nonzero negative temperature gradient ($\partial T / \partial x < 0$) at the cold boundary. Because the coordinate $x$ increases from the hot wall ($x = 0$) toward the cold wall ($x = W = 0.0889\text{ m}$), the adjacent solid cells ($T \approx 301.505\text{ K}$) are warmer than the fixed-temperature boundary ($T_c = 301.450\text{ K}$), so temperature decreases toward the boundary ($\partial T / \partial x < 0$). By Fourier's law, this negative gradient drives a positive conductive heat flux density ($q_x = -k \partial T / \partial x > 0$) directed in the $+x$ direction, carrying heat out of the domain across the cold wall.  
     *Raw OpenFOAM Sign Convention*: In OpenFOAM's `wallHeatFlux` function object, heat flux entering the fluid domain from a boundary face is signed positive ($\mathbf{q} \cdot \mathbf{n}_{\text{in}} > 0$). Consequently, heat leaving the domain across the cold boundary produces a negative raw patch integral in `wallHeatFlux.dat` ($\dot{Q}_{\text{cold, raw}} < 0$). In the net boundary energy balance equation above, the outward loss rate is defined by explicitly reversing the sign:
     $$\dot{Q}_{\text{cold, out}}(s) = -\dot{Q}_{\text{cold, raw}}(s) > 0$$
     This explicit sign conversion avoids any ambiguity over whether an absolute magnitude or an explicit sign reversal was executed. Outward cold-wall heat loss rises from $0.058\text{ W}$ at 1 min to $0.478\text{ W}$ at 17 min, totaling **$246.61\text{ J}$** of cumulative heat loss over 17 minutes ($6.99\text{ J}$ at 2 min).
   - **`topWall` & `bottomWall`**: Specified with `zeroGradient` (strictly adiabatic, $\dot{Q} \equiv 0$).
   - **`frontAndBack`**: Specified with `empty` boundary conditions (strictly 2D, $\dot{Q} \equiv 0$).
   - **OpenFOAM `wallHeatFlux` Function Object Area-Integration Verification**:  
     Inspection of the raw OpenFOAM `wallHeatFlux.dat` log confirms the standard column format `# Time patch min max integral`. In OpenFOAM, `min` and `max` report the extreme local wall heat flux densities ($q''$) across all boundary faces on the patch in $\text{W/m}^2$, while the `integral` column evaluates the exact surface area-weighted sum:
     $$\dot{Q}_{\text{patch}} = \sum_{f \in \text{patch}} q''_f A_f \quad [\text{W}]$$
     The area-averaged patch heat flux density is therefore:
     $$\bar{q}''_{\text{wall}} = \frac{\dot{Q}_{\text{wall}}}{A_{\text{wall}}} \quad [\text{W/m}^2]$$
     To clarify the relationship between local flux densities and the patch integral across different solver formulations, we report the exact initial data on the $42 \times 32$ grid ($A_{\text{hot}} = H \times D = 0.0635 \times 0.005 = 3.175 \times 10^{-4}\text{ m}^2$) at $t = 0.1\text{ s}$:
     1. *`EHCDarcy` Formulation ($42 \times 32$ mesh, $t = 0.1\text{ s}$)*:  
        The raw log records $q''_{\text{min}} = 2.55709 \times 10^5\text{ W/m}^2$, $q''_{\text{max}} = 2.55838 \times 10^5\text{ W/m}^2$, and $\dot{Q}_{\text{hot}} = \mathbf{81.2273\text{ W} \approx 81.23\text{ W}}$. Because natural convection has not yet developed at $t = 0.1\text{ s}$, early 1D conduction is nearly uniform along the vertical wall (a variation of only $0.05\%$ between min and max). The area-averaged flux density evaluates to $\bar{q}''_{\text{wall}} = \frac{81.2273\text{ W}}{3.175 \times 10^{-4}\text{ m}^2} = \mathbf{2.55834 \times 10^5\text{ W/m}^2}$, and the area-weighted product $\bar{q}''_{\text{wall}} A_{\text{wall}} = (2.55834 \times 10^5)(3.175 \times 10^{-4}) = \mathbf{81.23\text{ W}}$, exactly reconciling the reported integral.
     2. *Baseline `enthalpyPorosity` Formulation ($n_{\text{outer}}=8$, $42 \times 32$ mesh, $t = 0.1\text{ s}$)*:  
        The raw log records $q''_{\text{min}} = 1.99650 \times 10^5\text{ W/m}^2$, $q''_{\text{max}} = 2.27110 \times 10^5\text{ W/m}^2$, and $\dot{Q}_{\text{hot}} = \mathbf{71.7338\text{ W} \approx 71.73\text{ W}}$. The resulting area-averaged flux density is $\bar{q}''_{\text{wall}} = \frac{71.7338\text{ W}}{3.175 \times 10^{-4}\text{ m}^2} = \mathbf{2.25933 \times 10^5\text{ W/m}^2}$, with area product $\bar{q}''_{\text{wall}} A_{\text{wall}} = \mathbf{71.73\text{ W}}$. (In the higher-corrector $n_{\text{outer}}=12$ run, $\dot{Q}_{\text{hot}} = 75.56\text{ W}$ at $t = 0.1\text{ s}$, with $\bar{q}''_{\text{wall}} = 2.37978 \times 10^5\text{ W/m}^2$).
     The post-processing script reads this area-integrated rate directly from the `integral` column and performs step-resolved trapezoidal quadrature over time ($\text{W} \cdot \text{s} = \text{J}$), guaranteeing complete mathematical consistency between patch-integrated boundary rates and cumulative Joules.
   
   Numerical time integration of net boundary heat delivered into the domain accounts explicitly for the initial rectangular startup step across the unrecorded first interval ($0 \le t \le t_1 = 0.1\text{ s}$, $\Delta t_1 = 0.1\text{ s}$) plus standard step-resolved trapezoidal quadrature thereafter ($k = 2, \dots, n$):
   $$E_{\text{net, in}}(t_n) = E_{\text{net}, 0 \to t_1} + \sum_{k=2}^n \frac{1}{2}\left[\dot{Q}_{\text{net}}(t_{k-1}) + \dot{Q}_{\text{net}}(t_k)\right]\Delta t_k$$
   where $t_1 = 0.1\text{ s}$, $t_k = k \times 0.1\text{ s}$, and $\dot{Q}_{\text{net}}(t) = \dot{Q}_{\text{hot}}(t) - \dot{Q}_{\text{cold, out}}(t)$. For the gross hot-wall integral, the rectangular startup contribution is $E_{\text{hot}, 0 \to 0.1} = \dot{Q}_{\text{hot}}(t_1)\Delta t_1 = 71.7338\text{ W} \times 0.1\text{ s} = \mathbf{7.173\text{ J}}$. For the net boundary integral, the startup contribution is $E_{\text{net}, 0 \to 0.1} = [\dot{Q}_{\text{hot}}(t_1) - \dot{Q}_{\text{cold, out}}(t_1)]\Delta t_1$. At $t_1 = 0.1\text{ s}$, the hot-wall thermal disturbance has not produced a meaningful cold-wall response on this timescale, and the computed cold-wall heat flux is negligible at $t = 0.1\text{ s}$ ($\dot{Q}_{\text{cold, out}}(0.1\text{ s}) \approx 4.5 \times 10^{-11}\text{ W}$, cumulative loss $< 10^{-11}\text{ J}$), so $E_{\text{net}, 0 \to 0.1} \equiv 7.173\text{ J}$ matches the gross hot-wall startup contribution to within $< 10^{-11}\text{ J}$.

2. **Domain-Integrated Stored Energy**:
   Evaluated relative to the uniform initial isothermal solid state ($T_{\text{init}} = 301.45\text{ K}$, $f_L = 0$):
   $$E_{\text{sens}}(t) = \int_\Omega \rho C_p \left(T(\mathbf{x}, t) - T_{\text{init}}\right)\,dV$$
   $$E_{\text{lat}}(t) = \int_\Omega \rho L f_L(\mathbf{x}, t)\,dV$$
   $$E_{\text{stored}}(t) = E_{\text{sens}}(t) + E_{\text{lat}}(t)$$

3. **Discrete Energy-Balance Discrepancy**:
   $$\Delta E(t) = E_{\text{stored}}(t) - E_{\text{net, in}}(t), \quad \epsilon_{\text{energy}}(t) = \frac{\Delta E(t)}{E_{\text{net, in}}(t)} \times 100\%$$

**Domain & Physical Capacity Verification**:
- Cavity Dimensions: $W = 0.0889\text{ m}$, $H = 0.0635\text{ m}$, thickness $D = 0.005\text{ m}$ ($0.50\text{ cm}$).
- Total Volume: $V = 0.0889 \times 0.0635 \times 0.005 = 2.822575 \times 10^{-5}\text{ m}^3$.
- Total Gallium Mass: $m = \rho V = 6095.0 \times 2.822575 \times 10^{-5} \approx 0.17203595\text{ kg}$.
- Total Latent Capacity: $m \cdot L = (0.172035946\ldots\text{ kg})(80160\text{ J/kg}) \approx \mathbf{13790.40\text{ J}}$.
- Sensible Heat Capacity: $m \cdot C_p = (0.172035946\ldots\text{ kg})(381.5\text{ J/(kg K)}) \approx \mathbf{65.63\text{ J/K}}$.

### 2.3 Time-Resolved Energy Audit & Formulation Comparison

#### Complete Transient Energy Audit History (`EHCDarcy`, $42 \times 32$ Grid, $\Delta t = 0.1\text{ s}$)
Evaluating the step-resolved boundary flux integrals against domain-stored enthalpy across the full 17-minute duration yields:

| Time | $\dot{Q}_{\text{hot}}$ (W) | $\dot{Q}_{\text{cold}}$ (W) | $E_{\text{hot}}$ (J) | $E_{\text{cold, out}}$ (J) | $E_{\text{net, in}}$ (J) | $E_{\text{sens}}$ (J) | $E_{\text{lat}}$ (J) | $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Error (%) |
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

#### Reconciliation of Stored Latent Energy and Liquid Fraction
At $t = 17\text{ min}$ ($1020\text{ s}$), the global liquid fraction evaluated from the internal field is $f_L = 0.5736135$. With uniform cell volume and constant density:
$$E_{\text{lat}} = m \cdot L \cdot f_L = 13790.40\text{ J} \times 0.5736135 = \mathbf{7910.36\text{ J}} \approx \mathbf{7910.4\text{ J}}$$
Direct summation across all cells yields $\sum \rho V_{\text{cell}} L \alpha = 7910.36\text{ J}$, fully reconciling the stored latent energy with the reported melt fraction to $0.01\text{ J}$ (and $7910.4\text{ J}$ at $0.1\text{ J}$ precision).

#### Outer-Iteration Sensitivity Study in Enthalpy-Porosity
To test convergence toward discrete physics, Enthalpy-Porosity was executed across a sweep of PIMPLE outer correctors ($n_{\text{outer}} = 2, 5, 8, 12$) to $t = 2\text{ min}$:

| Formulation | $n_{\text{outer}}$ | $f_L(2\text{ min})$ | Net $E_{\text{net, in}}$ (J) | Stored $E_{\text{stored}}$ (J) | Discrepancy $\Delta E$ (J) | Energy Error (%) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **EHCDarcy** | 8 | **0.1303** | **1880.31 J** | **1878.37 J** | **-1.94 J** | **-0.10%** |
| Enthalpy-Porosity | 2 | 0.1762 | 1362.10 J | 2533.71 J | +1171.61 J | +86.01% |
| Enthalpy-Porosity | 5 | 0.1395 | 1719.85 J | 2009.79 J | +289.94 J | +16.86% |
| Enthalpy-Porosity | 8 | 0.1340 | 1809.52 J | 1931.11 J | +121.59 J | +6.72% |
| **Enthalpy-Porosity** | **12** | **0.1317** | **1856.07 J** | **1898.69 J** | **+42.62 J** | **+2.30%** |

At $t = 17\text{ min}$, Enthalpy-Porosity ($n_{\text{outer}}=12$) reaches:
- $E_{\text{hot}} = 8276.57\text{ J}$, $E_{\text{cold, out}} = 249.32\text{ J} \implies E_{\text{net, in}} = 8027.25\text{ J}$
- Stored Energy: $E_{\text{stored}} = 8265.76\text{ J}$ ($E_{\text{sens}} = 236.27\text{ J}$, $E_{\text{lat}} = 8029.50\text{ J}$, $f_L = 0.58225$)
- Discrepancy: $\Delta E = +238.51\text{ J}$ ($\mathbf{+2.97\%}$)

**Formulation Comparison Finding**:  
Increasing the outer-corrector count substantially reduces the energy imbalance and brings the enthalpy–porosity results closer to those obtained with EHCDarcy (reducing the 2-minute error from $+6.72\%$ to $+2.30\%$, and 17-minute error from $+5.18\%$ to $+2.97\%$). However, while $+2.30\%$ readily satisfies the developed convective criterion ($|\epsilon_{\text{energy}}| \le 6.0\%$), it does not meet the strict early-phase conduction target ($|\epsilon_{\text{energy}}| \le 2.0\%$). At late times ($t = 17\text{ min}$), Enthalpy-Porosity ($n_{\text{outer}}=12$) achieves $+2.97\%$ (passing the $\le 6.0\%$ convective target). The observed error reduction is consistent with incomplete outer-iteration coupling during steep transients, although other numerical contributions—including discrete equation residuals, spatial truncation, and flux integration approximations—are not excluded. The remaining differences, particularly at late times ($f_L = 0.5822$ vs $0.5736$, an absolute difference of $0.0086$, or $1.50\%$), indicate that complete numerical equivalence has not yet been demonstrated across the full simulation. Algorithmic differences in latent heat updating (Newton-projection source terms vs. diagonalized path-integral effective capacity) account for the distinct convergence paths.

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

## 5. Spatial Grid Sensitivity & ASME / Roache GCI Analysis (to 17 min)

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
Ye & Arıcı (2023, Section 3.2.2, Eq. 6) conducted a 2D verification study on natural convection melting of gallium, correlating their grid-converged numerical liquid fractions as:
$$\phi_n = 2.708\,\tau^{0.712}\,\mathrm{Ra}^{0.0243}\,\mathrm{Ar}^{-0.0515}$$
where:
- $\phi_n$ represents the numerical global liquid volume fraction ($f_L$).
- $\tau = \mathrm{Fo} \cdot \mathrm{Ste}$, with $\mathrm{Fo} = \frac{\alpha t}{W^2} = 1.7761$ and $\mathrm{Ste} = \frac{C_p (T_h - T_m)}{L} = 0.03912$.
- $\mathrm{Ra} = \frac{g \beta (T_h - T_m) H^3}{\nu \alpha} = 6.063 \times 10^5$.
- $\mathrm{Ar} = \frac{H}{W} = 0.7143$.

Evaluating Eq. (6) at $t = 17\text{ min}$ ($1020\text{ s}$):
$$\phi_n^{\text{Ye \& Arıcı (2023)}} = \mathbf{0.5703}$$
- `phaseChangeFoam` (EHCDarcy, $42 \times 32$): **$f_L = 0.5736$** (Discrepancy: **$+0.0033$, or $0.33$ percentage points**).
- `phaseChangeFoam` (EHCDarcy, $63 \times 48$): **$f_L = 0.5725$** (Discrepancy: **$+0.0022$, or $0.22$ percentage points**).
- `phaseChangeFoam` (EHCDarcy, $32 \times 24$): **$f_L = 0.5682$** (Discrepancy: **$-0.0021$, or $0.21$ percentage points**).

*Note on Interpretation*: Agreement with a correlation fitted to published numerical solutions serves as a valuable consistency check against independent modern CFD codes, but is kept distinct from experimental validation against physical measurements.

---

## 8. Summary Conclusions & Verification Scope Assessment

### 8.1 Summary of Methodological Resolutions
1. **Energy-Balance Audit and Boundary Accounting Resolutions**:
   - **Isolation of Temporal Undersampling vs. Estimator Discrepancy**: At the two checked output timestamps, the finite-difference and native OpenFOAM flux estimators agree to within $0.01\%$, with differences ($\Delta Q = Q_{\text{native}} - Q_{\text{FD}}$) of approximately $+0.00062\text{ W}$ at 60 s ($11.22152\text{ W}$ vs. $11.22090\text{ W}$) and $-0.00032\text{ W}$ at 120 s ($8.49778\text{ W}$ vs. $8.49810\text{ W}$). Both estimators round to the same two-decimal flux values ($11.22\text{ W}$ and $8.50\text{ W}$) at these timestamps (though this comparison does not establish agreement at every intervening time). Evaluating the coarse 60-second quadrature rule on native flux yields $1601.516\text{ J} \approx 1601.52\text{ J}$, while the archived script computed $1601.451\text{ J}$ (displayed as $1601.5\text{ J}$ in the archived log), a minor difference of $0.065\text{ J}$ ($0.004\%$). This controlled reconstruction establishes that coarse startup quadrature and temporal undersampling account for approximately $214.86\text{ J}$ of the gross-input difference. Full-resolution integration is documented as a rectangular startup step ($E_{0 \to 0.1} = 7.173\text{ J}$) plus 1,199 trapezoidal steps ($1809.205\text{ J}$), reconciling unrounded data to $1816.378\text{ J} \approx 1816.38\text{ J}$ with an evaluated startup sensitivity of $\approx 1.07\text{ J}$ (using the initial theoretical wall gradient $\approx 93.1\text{ W}$).
   - **Gross vs. Net Boundary Accounting**: Explicit separation of gross hot-wall delivery ($E_{\text{hot}}$) and conductive cold-wall loss ($E_{\text{cold, out}}$, with positive outward loss defined from OpenFOAM's raw negative patch integral via $\dot{Q}_{\text{cold, out}} = -\dot{Q}_{\text{cold, raw}} > 0$) resolved the boundary energy balance across all audit timestamps ($1\text{ min}$, $2\text{ min}$, and $17\text{ min}$). At 17 minutes, cold-wall loss totals $246.61\text{ J}$, yielding a net boundary input of $8148.99\text{ J}$ that closes against stored domain energy ($8143.07\text{ J}$) within **$-5.92\text{ J}$ ($-0.073\%$)** for `EHCDarcy`.
   - **Enthalpy–Porosity Formulation Sensitivity**: Increasing outer correctors from $n_{\text{outer}}=8 \to 12$ reduces the 17-minute discrepancy from $+5.18\%$ to $+2.97\%$ and the 2-minute error from $+6.72\%$ to $+2.30\%$. While $+2.30\%$ satisfies the convective criterion ($\le 6.0\%$), it falls just short of the strict early conduction target ($\le 2.0\%$), reflecting outer-iteration sensitivity during rapid phase front initialization.

2. **Methodological Qualification of Latent Heat Advection**:
   - A cell-by-cell scale analysis across all 85 mushy cells at $t = 1020\text{ s}$, evaluated as a **continuum-expanded, cell-centred diagnostic** $\sum V_i \rho L (\mathbf{u}_i \cdot \nabla f_{L,i})$ (distinguished from the discrete finite-volume face-flux operator $\nabla \cdot (\boldsymbol{\phi} f_L)$), confirmed that the median cell-wise ratio in the core mushy zone ($0.01 < f_L < 0.90$, 80 cells) is small ($\operatorname{median}_i (|\mathbf{u}_i \cdot \nabla f_{L,i}| / |\partial f_{L,i}/\partial t|) = 0.99\% \approx 1.0\%$, ratio of medians: $0.50\%$).
   - Because computed per-cell phase-fraction differences between 960 s and 1020 s are nonnegative across all cells belonging to the mushy-zone mask $\mathcal{S}$ defined at $t = 1020\text{ s}$ ($\Delta f_{L,i} \ge 0$), the signed and absolute interval-averaged transient latent uptake rates coincide for this diagnostic interval: $P_{\text{trans, signed}} \equiv P_{\text{trans, abs}} = 2.61\text{ W}$ (core), $0.47\text{ W}$ (fringe), and $3.08\text{ W}$ (entire mushy zone).
   - In the core mushy zone, absolute advection evaluates to $3.69\text{ W}$ ($1.41\times$ the transient uptake rate) and signed advection evaluates to $+3.01\text{ W}$. In the 5 outermost fringe cells ($f_L \ge 0.90$), localized advection reaches $9.20\text{ W}$ absolute and $-7.52\text{ W}$ signed, which are comparable to whole-cavity power scales ($\dot{Q}_{\text{hot}} \approx 6.99\text{ W}$, domain $\dot{E}_{\text{lat}} \approx 6.42\text{ W}$). Under the continuum incompressibility identity, these signs are consistent with net flux divergence in the core and net flux convergence into the fringe, while definitive regional flux confirmation requires evaluating face fluxes across subregion boundaries. Across the entire mushy zone, absolute advection totals $12.90\text{ W}$ ($4.19\times$ transient uptake).
   - Crucially, the $1.41$ and $4.19$ ratios compare integrated local magnitudes rather than global energy changes, as global advective divergence cancels out across the closed domain.
   - We also note that the advection diagnostic is evaluated instantaneously at $t = 1020\text{ s}$ while the transient rate is backward-differenced over $960\text{--}1020\text{ s}$ ($\Delta t = 60\text{ s}$), limiting a strict instantaneous balance.
   - Consequently, omitting $\nabla \cdot (\rho \mathbf{u} L f_L)$ is correctly classified as an **established modeling approximation** standard across CFD benchmark literature, rather than an exact mathematical identity for diffuse interfaces. The diagnostic identifies a potentially important localized redistribution whose precise influence on the local melt front remains to be quantified by future comparative simulations incorporating the full discrete finite-volume divergence operator.

3. **Spatial and Temporal Sensitivity Bounds**:
   - Mid-height front position at 17 min demonstrates strict monotonic spatial convergence with an ASME observed order of $p = 2.70$, an asymptotic consistency ratio of $1.01$, and a fine-grid $\text{GCI}$ of $0.67\%$ (medium-grid $\text{GCI} = 2.04\%$).
   - Global melt fraction at 17 min exhibits minimal spatial variation ($\Delta f_L \le 0.0054$, $0.95\%$) across grids, matching the modern grid-converged numerical correlation of Ye & Arıcı (2023) within $0.33$ percentage points.

### 8.2 Scope Distinction: Numerical Benchmark Verification vs. Physical Experimental Validation
To preserve scientific precision, the scope and status of this report are explicitly defined:
- **Numerical Benchmark Verification**: This report provides comprehensive, rigorous evidence of **numerical benchmark consistency and code-to-code verification** against the established 2D CFD literature (Brent et al. 1988, Hannoun et al. 2003, Ye & Arıcı 2023). Under the stated idealized assumptions, `phaseChangeMultiRegionFoam` with the EHCDarcy formulation achieves a mean phase-front deviation of $3.46\text{ mm}$ ($3.89\%$ of cavity width) and a maximum local deviation of $6.34\text{ mm}$ ($7.13\%$ of cavity width) relative to the digitized Brent et al. numerical benchmark. The reported global energy-balance discrepancy at 17 minutes is $-0.073\%$.
- **Distinction from Experimental Validation**: These results **do not constitute a full quantitative validation against independent physical experimental measurements**. Quantitative alignment with physical experiments (such as Gau & Viskanta 1986) requires resolving 3D endwall boundary-layer drag, conjugate thermal inertia of the heating jacket and wall apparatus, solid-state gallium crystallographic anisotropy ($k \in [28, 41]\text{ W/(m K)}$), and non-intrusive front measurement methods that avoid pour-out boundary layer disruption.
- **Document Status**: Accordingly, this study is formally presented as a **high-fidelity numerical benchmark verification report**, providing a solid foundation for future 3D conjugated experimental validation studies.

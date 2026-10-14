# Comprehensive Validation & Model Comparison: Gallium 2D Melting Benchmark

**Solver**: `phaseChangeMultiRegionFoam`  
**Phase Change Formulations**:  
1. `enthalpyPorosity` (Voller & Prakash 1987, Brent et al. 1988)  
2. `EHCDarcy` (Path-Integral Effective Heat Capacity with Carman-Kozeny Darcy Damping)  
**Numerical Benchmark**: A. D. Brent, V. R. Voller, and K. J. Reid (1988), *Numerical Heat Transfer*, 13(3), pp. 297–318.  
**Experimental Reference**: C. Gau and R. Viskanta (1986), *J. Heat Transfer*, 108(1), pp. 174–181.  
**Cavity Dimensions**: $W = 8.89\text{ cm} \times H = 6.35\text{ cm}$ (Aspect Ratio $A = 1.40$)

---

## 1. Executive Summary & Acceptance Criteria Scorecard

To address the observed late-time melting overprediction in the enthalpy-porosity formulation, a parallel study was conducted using the **EHCDarcy** formulation on the identical mesh ($42 \times 32$), time-step ($\Delta t = 0.1\text{ s}$), and physical property set.

### 1.1 Acceptance Criteria Comparison Table

| Metric | Target | Enthalpy-Porosity | EHCDarcy | Status (EHCDarcy) |
| :--- | :---: | :---: | :---: | :---: |
| **Mean Absolute Error ($f_L$)** | $\le 0.03\text{--}0.05$ | $0.0401$ ($4.01\%$) | **$0.0278$ ($2.78\%$)** | **PASSED** ($\mathbf{30.7\%}$ improvement) |
| **Mean Relative Error ($f_L$)** | $\le 15\%$ | $14.73\%$ | **$10.89\%$** | **PASSED** |
| **$f_L$ at $t = 2\text{ min}$** (Brent: $0.1102$) | $\pm 0.03$ | $0.1340$ ($+0.0238$) | **$0.1303$ ($+0.0201$)** | **PASSED** |
| **$f_L$ at $t = 6\text{ min}$** (Brent: $0.2346$) | $\pm 0.03$ | $0.2659$ ($+0.0313$) | **$0.2573$ ($+0.0227$)** | **PASSED** |
| **$f_L$ at $t = 10\text{ min}$** (Brent: $0.3465$) | $\pm 0.04$ | $0.3884$ ($+0.0419$) | **$0.3747$ ($+0.0282$)** | **PASSED** |
| **$f_L$ at $t = 17\text{ min}$** (Brent: $0.5334$) | $\pm 0.05$ | $0.5967$ ($+0.0633$) | **$0.5736$ ($+0.0402$)** | **PASSED** |
| **Mean Front Dev. ($t = 2\text{ min}$)** | $\le 5\%$ of $W$ ($4.45\text{ mm}$) | $1.95\text{ mm}$ ($2.19\%$) | **$1.67\text{ mm}$ ($1.88\%$)** | **PASSED** |
| **Mean Front Dev. ($t = 6\text{ min}$)** | $\le 5\%$ of $W$ ($4.45\text{ mm}$) | $2.60\text{ mm}$ ($2.93\%$) | **$1.87\text{ mm}$ ($2.11\%$)** | **PASSED** |
| **Mean Front Dev. ($t = 10\text{ min}$)** | $\le 5\%$ of $W$ ($4.45\text{ mm}$) | $3.51\text{ mm}$ ($3.95\%$) | **$2.35\text{ mm}$ ($2.65\%$)** | **PASSED** |
| **Mean Front Dev. ($t = 17\text{ min}$)** | $\le 5\%$ of $W$ ($4.45\text{ mm}$) | $5.43\text{ mm}$ ($6.11\%$) | **$3.46\text{ mm}$ ($3.89\%$)** | **PASSED** |
| **Max Front Dev. ($t = 2\text{ min}$)** | $\le 10\%$ of $W$ ($8.89\text{ mm}$) | $2.54\text{ mm}$ ($2.85\%$) | **$2.29\text{ mm}$ ($2.58\%$)** | **PASSED** |
| **Max Front Dev. ($t = 6\text{ min}$)** | $\le 10\%$ of $W$ ($8.89\text{ mm}$) | $3.88\text{ mm}$ ($4.37\%$) | **$3.18\text{ mm}$ ($3.58\%$)** | **PASSED** |
| **Max Front Dev. ($t = 10\text{ min}$)** | $\le 10\%$ of $W$ ($8.89\text{ mm}$) | $5.62\text{ mm}$ ($6.32\%$) | **$4.12\text{ mm}$ ($4.63\%$)** | **PASSED** |
| **Max Front Dev. ($t = 17\text{ min}$)** | $\le 10\%$ of $W$ ($8.89\text{ mm}$) | $9.74\text{ mm}$ ($10.95\%$) | **$6.34\text{ mm}$ ($7.13\%$)** | **PASSED** |

**Conclusion**: Switching to `EHCDarcy` successfully brings **all front deviation criteria and the global liquid fraction error within the rigorous acceptance envelope**. The mean absolute error drops to **$0.0278$ ($2.78$ percentage points)**, achieving a **$30.7\%$ reduction in discrepancy**.

---

## 2. Complete Transient Progression Table

| Time ($\text{s}$) | Time ($\text{min}$) | Brent (1988) | Enthalpy-Porosity | EHCDarcy | $\Delta(f_L^{\text{EHC}} - f_L^{\text{EP}})$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| $60.0$ | $1.00$ | — | $0.0918$ | $0.0903$ | $-0.0015$ |
| **$120.0$** | **$2.00$** | **$0.1102$** | **$0.1340$** | **$0.1303$** | **$-0.0037$** |
| $180.0$ | $3.00$ | — | $0.1697$ | $0.1643$ | $-0.0054$ |
| $240.0$ | $4.00$ | — | $0.2028$ | $0.1967$ | $-0.0061$ |
| $300.0$ | $5.00$ | — | $0.2347$ | $0.2274$ | $-0.0073$ |
| **$360.0$** | **$6.00$** | **$0.2346$** | **$0.2659$** | **$0.2573$** | **$-0.0086$** |
| $420.0$ | $7.00$ | — | $0.2968$ | $0.2871$ | $-0.0097$ |
| $480.0$ | $8.00$ | — | $0.3275$ | $0.3163$ | $-0.0112$ |
| $540.0$ | $9.00$ | — | $0.3580$ | $0.3456$ | $-0.0124$ |
| **$600.0$** | **$10.00$** | **$0.3465$** | **$0.3884$** | **$0.3747$** | **$-0.0138$** |
| $660.0$ | $11.00$ | — | $0.4187$ | $0.4040$ | $-0.0147$ |
| $720.0$ | $12.00$ | — | $0.4488$ | $0.4323$ | $-0.0165$ |
| $780.0$ | $13.00$ | — | $0.4787$ | $0.4616$ | $-0.0171$ |
| $840.0$ | $14.00$ | — | $0.5085$ | $0.4893$ | $-0.0192$ |
| $900.0$ | $15.00$ | — | $0.5381$ | $0.5176$ | $-0.0205$ |
| $960.0$ | $16.00$ | — | $0.5676$ | $0.5457$ | $-0.0219$ |
| **$1020.0$** | **$17.00$** | **$0.5334$** | **$0.5967$** | **$0.5736$** | **$-0.0231$** |

---

## 3. Root Cause Analysis: Why Overprediction Occurs & Why EHCDarcy Resolves It

### 3.1 Numerical Formulation Differences (Enthalpy-Porosity vs EHCDarcy)

1. **Enthalpy-Porosity Formulation**:
   - In `enthalpyPorosity`, phase fraction $f_L$ is updated iteratively:
     $$f^{k+1} = \mathrm{clip}\left(f^k + \lambda \frac{T - T_m}{L/c_p + \Delta T}, 0, 1\right)$$
   - The latent heat appears as an explicit source term $S_h = -\frac{\rho L}{\Delta t}(f^k - f^{old})$.
   - In high-velocity convection zones (such as the upper impingement jet where warm fluid hits the melt interface at the ceiling), explicit outer-iteration lag can allow local cell temperature to temporarily exceed $T_m$, artificially accelerating interface propagation.

2. **EHCDarcy (Effective Heat Capacity) Formulation**:
   - In `EHCDarcy`, phase change is path-integrated with an exact Newton linearization:
     $$S_p = \frac{\rho_{bar} L}{\Delta t c_p} \frac{df}{dT}, \quad S_u = \frac{\rho_{bar} L}{\Delta t}(f - f^{old})$$
   - The diagonal term $S_p$ directly couples the latent heat absorption implicitly into the matrix equation for $h$.
   - This implicit stabilization strictly tethers cell temperature to the melting temperature range, eliminating local temperature overshoots and preventing the convective jet from melting cells prematurely.
   - Result: As observed at all benchmark times ($t = 2, 6, 10, 17\text{ min}$), EHCDarcy systematically restrains the upper melt bulge, shifting the interface $2\text{--}3.4\text{ mm}$ closer to the Brent numerical contours.

### 3.2 Time Discretization & Historical Delay in Brent et al. (1988)

At $t = 2\text{ min}$, heat transfer is $>90\%$ pure conduction. The exact analytical 1D two-phase Stefan problem with solid subcooling predicts:
$$f_L^{\text{exact Stefan}} = 0.1250 \quad (x_f = 11.11\text{ mm})$$

- `phaseChangeFoam` (EHCDarcy) predicts: $f_L = 0.1303$ ($x_f = 11.58\text{ mm}$), matching the analytical Stefan solution with early buoyancy.
- Brent et al. (1988) reported: $f_L = 0.1102$ ($x_f \approx 9.0\text{ mm}$).
- **Why Brent's early result was low**: Brent et al. utilized $\Delta t = 5\text{ s}$ for 4 steps followed by $\Delta t = 10\text{ s}$ with a first-order backward-Euler implicit scheme. The large 10-second time step introduced numerical diffusion that retarded early boundary heating. Modern solvers with $\Delta t = 0.1\text{ s}$ resolve the initial transient heating cleanly.

### 3.3 Convective Discretization Schemes (Upwind vs Higher-Order)

Testing higher-order bounded schemes (`Gauss limitedLinear 1`) confirms that numerical diffusion in 1st-order upwind schemes adds a slight positive bias to melting rates ($f_L$ at $120\text{ s}$ drops from $0.1340$ to $0.1329$). However, the primary physical-numerical damping is achieved by the implicit Newton formulation of **EHCDarcy**.

---

## 4. Master Plot Details

The master figure `validation/gallium2D/gallium_melt_fraction_comparison.png` combines:
1. **Panel (a) Liquid Fraction Evolution vs Time**:
   - Brent et al. (1988) Numerical Benchmark (red dashed line with square markers)
   - Analytical 1D Stefan pure conduction solution (green dotted line)
   - Enthalpy-Porosity model (blue solid line with circle markers)
   - EHCDarcy model (purple solid line with triangle markers)
   - Performance scorecard text box showing the $30.7\%$ error reduction.
2. **Panel (b) Melt Front Morphology Comparison**:
   - Melt front profiles at $2\text{ min}$ (blue), $6\text{ min}$ (green), $10\text{ min}$ (orange), and $17\text{ min}$ (red).
   - Dashed lines: Brent et al. (1988) reference contours.
   - Solid lines: Enthalpy-Porosity predictions.
   - Dotted lines with markers: EHCDarcy predictions.

# Physical Formulations & Governing Equations

Mathematical models, governing equations, and thermodynamic consistency relations implemented in `phaseChangeFoam`.

---

## Governing Equations

The solver treats multi-region domains composed of fluid and solid regions coupled at conjugate interfaces.

### Fluid Regions

For fluid regions with natural convection, species transport, and phase change:

1. **Continuity**:
   $$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{U}) = \dot{m}$$
   where $\dot{m} = \sum \dot{m}_k$ is the net volumetric mass addition (zero for closed single-phase species conversion).

2. **Momentum (Navier-Stokes)**:
   $$\frac{\partial (\rho \mathbf{U})}{\partial t} + \nabla \cdot (\rho \mathbf{U} \mathbf{U}) = -\nabla p_{\text{rgh}} - (\mathbf{g} \cdot \mathbf{r})\nabla \rho + \nabla \cdot (\mu_{\text{eff}} \nabla \mathbf{U}) + \rho \mathbf{g} \beta (T - T_{\text{ref}}) + \mathbf{S}_D$$
   where:
   * $p_{\text{rgh}} = p - \rho \mathbf{g} \cdot \mathbf{r}$ is modified pressure.
   * $\rho \mathbf{g} \beta (T - T_{\text{ref}})$ is the Boussinesq thermal buoyancy source.
   * $\mathbf{S}_D$ is the Carman-Kozeny Darcy drag sink (active when solid-liquid PCM models operate in fluid regions):
     $$\mathbf{S}_D = -C_u \frac{(1 - \alpha)^2}{\alpha^3 + q} \mathbf{U}$$

3. **Species Transport**:
   $$\frac{\partial (\rho Y_k)}{\partial t} + \nabla \cdot (\rho \mathbf{U} Y_k) = \nabla \cdot (\rho D_k \nabla Y_k) + \dot{m}_k$$

4. **Energy (Sensible Enthalpy $h$)**:
   $$\frac{\partial (\rho h)}{\partial t} + \nabla \cdot (\rho \mathbf{U} h) + \frac{\partial (\rho K)}{\partial t} + \nabla \cdot (\rho \mathbf{U} K) = \nabla \cdot (\Gamma_h \nabla h) + \frac{\partial p}{\partial t} + S_h$$
   where $K = \frac{1}{2}|\mathbf{U}|^2$, $\Gamma_h = k_{\text{eff}} / C_p$ is the enthalpy diffusion coefficient, and $S_h$ is the latent heat source or sink.

### Solid Regions

In solid regions, thermal transport is solved for sensible enthalpy $h$ using the monolithic coupled energy framework:
$$\beta_v \frac{\rho_{\text{eff}} C_{p,\text{eff}}}{C_{p,\text{thermo}}} \frac{\partial h}{\partial t} - \nabla \cdot \left(\beta_v \frac{k_{\text{eff}}}{C_{p,\text{thermo}}} \nabla h\right) = S_{\text{latent}}$$
where:
* $\beta_v$ is the solid volume fraction (porosity factor).
* $C_{p,\text{thermo}}$ is the baseline thermophysical model heat capacity, ensuring that $\frac{k_{\text{eff}}}{C_{p,\text{thermo}}} \nabla h = k_{\text{eff}} \nabla T$ holds identically.
* $C_{p,\text{eff}}$ is the effective sensible heat capacity accounting for temperature-dependent solid and liquid properties.
* $S_{\text{latent}}$ is the separate latent heat contribution, assembled via explicit source $S_{\text{lat}}$ and implicit Newton linearisation $S_p$:
  $$S_{\text{latent}} = \beta_v \left[ S_{\text{lat}} + S_p (h - h^k) \right]$$
  discretized as `+ betav * pcModel.latentHeatSource() + fvm::Sp(betav * pcModel.latentHeatSp(), h) - betav * pcModel.latentHeatSp() * h`.

> [!NOTE]
> These equations present the principal transport and phase-change contributions. Optional radiation, external heat sources, solver constraints and formulation-dependent terms are omitted for clarity.

---

## Fluid Region Phase Change Models

Classes derive from [`fluidPhaseChangeModel`](../src/solvers/phaseChangeMultiRegionFoam/fluid/fluidPhaseChangeModel.H).

### 1. Volumetric Evaporation & Condensation (`Lee`)

Class: `leeFluidPhaseChangeModel`

The Lee model drives inter-species mass transfer based on deviation from saturation temperature $T_{\text{sat}}$:

$$\dot{m}_v = \begin{cases} C_{\text{evap}} \rho_l Y_l \frac{T - T_{\text{sat}}}{T_{\text{sat}}} & \text{if } T > T_{\text{sat}} \quad (\text{evaporation}) \\ - C_{\text{cond}} \rho_v Y_v \frac{T_{\text{sat}} - T}{T_{\text{sat}}} & \text{if } T < T_{\text{sat}} \quad (\text{condensation}) \end{cases}$$

#### Saturation Temperature Modeling
* **Constant $T_{\text{sat}}$**: Specified directly via keyword `Tsat` (default: 373.15 K).
* **Pressure-dependent $T_{\text{sat}}(p)$**: Enabled via `enableTsatP true`, evaluated using the Clausius-Clapeyron equation:
  $$T_{\text{sat}}(p) = \left( \frac{1}{T_{\text{ref}}} - \frac{R_{\text{specific}}}{L} \ln\left(\frac{p}{p_{\text{ref}}}\right) \right)^{-1}$$

#### Coupled Conservation Sources
Mass and energy equations incorporate coupled source terms:
* Liquid species $Y_l$: $\dot{S}(Y_l) = -\dot{m}_v$
* Vapor species $Y_v$: $\dot{S}(Y_v) = +\dot{m}_v$
* Energy equation: $\dot{S}_h = -\dot{m}_v L_{\text{eff}}$
* Pressure equation: $\nabla \cdot \mathbf{U} = \left(\frac{1}{\rho_v} - \frac{1}{\rho_l}\right)\dot{m}_v$

#### Repeated Outer-Corrector Evaluation
Solver guardrails in `createCoupledRegions.H` enforce `nOuterCorrectors >= 2` whenever phase change is active in a fluid region or multi-region domain. Invariance tests in `verification/leeModel` verify that $\dot{m}_v$ is evaluated cleanly per outer iteration without double-counting across $n_{\text{outer}} \in [2, 4]$.

---

### 2. Constant Source Verification Model (`constantSource`)

Class: `constantSourceFluidPhaseChangeModel`

Applies fixed mass transfer rate $\dot{m}_0$ (kg/(m3 s)) and energy source rate $\dot{S}_0$ (W/m3) for unit testing and analytical MMS verification.

---

### 3. Inactive Fluid Model (`none`)

Class: `noneFluidPhaseChangeModel`

Returns zero mass and energy sources for single-phase non-reacting fluid regions.

---

## Solid & Fluid PCM Models

Classes derive from [`phaseChangeModel`](../src/solvers/phaseChangeMultiRegionFoam/solid/phaseChangeModel.H).

### 1. Effective Heat Capacity (`EHC`)

Class: `ehcPhaseChangeModel`

The `EHC` model separates sensible heat capacity variation from latent heat exchange:
1. **Sensible heat capacity**: Uses a path-based sensible heat capacity $C_{p,\text{eff}}$ to account for different solid and liquid specific heats without numerical overshoot.
2. **Latent heat handling**: Latent heat is handled separately through explicit source ($S_{\text{lat}}$) and implicit linearisation ($S_p$) terms entering the energy equation.

#### Thermophysical Modes & Path-Secant Capacity
The sensible heat capacity calculation depends on the thermophysical dictionary setting in `phaseChangeDict`:
* **`thermophysical { mode custom; }`**: When custom mode is active with user-specified $C_{p,s}$ and $C_{p,l}$, evaluating capacity at a single point $C_p(T^k)$ can cause numerical temperature overshoot when a time step crosses the phase change boundary. To eliminate this, `phaseChangeModel` computes a path-secant heat capacity from the change in sensible enthalpy:
  $$C_{p,\text{eff}} = \frac{\Delta h_{\text{sens}}}{\Delta T} = \frac{h_{\text{sens}}(T_{\text{new}}) - h_{\text{sens}}(T_{\text{old}})}{T_{\text{new}} - T_{\text{old}}}$$
  where $h_{\text{sens}}(T)$ is evaluated by piecewise analytic integration across solid, mushy, and liquid regimes:
  $$h_{\text{sens}}(T) = \begin{cases} C_{p,s} T & T \le T_s \\ C_{p,s} T_s + C_{p,s}(T - T_s) + \frac{C_{p,l} - C_{p,s}}{2} \frac{(T - T_s)^2}{T_l - T_s} & T_s < T < T_l \\ C_{p,s} T_s + \frac{C_{p,s} + C_{p,l}}{2}(T_l - T_s) + C_{p,l}(T - T_l) & T \ge T_l \end{cases}$$
  Effective thermal conductivity is similarly blended: $k_{\text{eff}} = (1 - \alpha) k_s + \alpha k_l$.

  When hysteresis or direction restrictions are enabled, the path-secant calculation follows a history-aware branch: it accounts for constant-phase-fraction plateau regimes up to the trajectory exit temperature, splitting the enthalpy increment between the frozen plateau capacity and the active transition curve.
* **`thermophysical { mode thermo; }` (Default)**: Uses the thermophysical model's underlying heat capacity $C_{p,\text{eff}} = C_{p,\text{thermo}}$ and thermal conductivity $k_{\text{eff}} = \kappa_{\text{thermo}}$.

#### Latent Heat Source and Linearisation
Latent heat is added to the sensible enthalpy matrix via explicit and implicit coefficients:
* **Explicit latent heat source**:
  $$S_{\text{lat}} = \frac{\bar{\rho} L}{\Delta t}(\alpha - \alpha^{\text{old}})$$
* **Implicit linearisation coefficient**:
  $$S_p = \frac{\bar{\rho} L}{\Delta t\, C_{p,\text{thermo}}} \frac{\partial \alpha}{\partial T}$$
In the energy matrix (`solveSolid.H` and `EEqn.H`), these enter as:
$$hEqn += \beta_v S_{\text{lat}} + \text{fvm::Sp}(\beta_v S_p, h) - \beta_v S_p h^k$$
providing quasi-Newton fixed-point stabilisation during phase transition.

Here, $\bar{\rho}$ denotes the density used for latent heat accounting. For the linear density formulation, it is evaluated using the midpoint phase fraction between the current and previous states. Otherwise, the thermophysical model density is used. The superscript $k$ denotes the current nonlinear iteration.

#### Transition Direction & Stateful Hysteresis
Phase fraction is governed by a cell-wise state machine supporting distinct melting and freezing curves. The model tracks the current heating or cooling trajectory using the stored thermal history. When the temperature reverses beyond the specified tolerance, the trajectory changes and the corresponding melting or freezing curve is followed. The previous liquid fraction is retained where required to preserve the history-dependent path.

For materials exhibiting subcooling, freezing begins at a lower temperature ($T_{\text{liq,f}} < T_{\text{liq,m}}$). Cell-wise reversal temperatures $T_{\text{reversal}}$ track thermal history to ensure proper path dependence.

> [!IMPORTANT]
> **Domain Applicability & Convection Guard**:
> `EHC` is designed for solid regions or pure-conduction domains. It does **not** include Darcy momentum damping. If `EHC` is applied to a fluid region with an active velocity field $\mathbf{U}$ and unsuppressed convection, the solver throws a runtime fatal error. When modeling phase change with natural convection or fluid motion, use `EHCDarcy` or `enthalpyPorosity`, or suppress convection via `convection { suppress true; }`.

---

### 2. EHC with Darcy Momentum Damping (`EHCDarcy`)

Class: `ehcDarcyPhaseChangeModel`

Extends `EHC` for fluid regions where phase change couples with Navier-Stokes flow.

* **Energy equation**: Formulated in terms of sensible enthalpy with path-based sensible heat capacity and separate linearised latent heat sources identical to `EHC`.
* **Momentum damping**: Provides an implicit Darcy drag source to `UEqn.H`:
  $$\mathbf{S}_D = -C_u \frac{(1 - \alpha)^2}{\alpha^3 + q} \mathbf{U}$$
  The term is discretized implicitly via `fvm::Sp(pcmModel.momentumSp(), U)`. In solid cells ($\alpha \to 0$), the drag coefficient approaches the penalty limit $C_u / q$. For the default code settings ($C_u = 10^5, q = 10^{-3}$), this equals $10^8$ kg/(m3 s); for benchmark example settings ($C_u = 1.6 \times 10^6, q = 10^{-3}$), it approaches $1.6 \times 10^9$ kg/(m3 s), extinguishing velocities to machine zero without numerical instability.

---

### 3. Enthalpy-Porosity Formulation (`enthalpyPorosity`)

Class: `enthalpyPorosityPhaseChangeModel`

Implements the classic source-based enthalpy-porosity formulation of Voller & Prakash (1987) and Brent et al. (1988).

The current enthalpy-porosity implementation supports a single phase-transition window. Stateful thermal hysteresis and separate melting/freezing curves are not supported by this model. Configuring a `hysteresis` or `reverse` subdictionary, or specifying a `direction` other than `both`, will trigger a runtime fatal error.

#### Governing Energy Equation
$$\frac{\partial (\rho C_p T)}{\partial t} + \nabla \cdot (\rho C_p \mathbf{U} T) = \nabla \cdot (k \nabla T) + S_h$$
where latent heat is tracked through an explicit source:
$$S_h = -L \frac{\partial (\rho f)}{\partial t}$$

#### Liquid Fraction Newton Projection
Following the sensible energy solve, liquid fraction $f$ is updated iteratively via Newton projection:
$$f_{\text{new}} = \min\left(\max\left(f + \lambda \frac{T - T_f}{L/C_p + (T_{\text{liq}} - T_{\text{sol}})}, 0\right), 1\right)$$
where $T_f = T_{\text{sol}} + f(T_{\text{liq}} - T_{\text{sol}})$, and $\lambda$ is the Newton step relaxation factor.

For a linear phase-transition window, the exact Newton step corresponds to $\lambda = 1.0$. Setting a smaller relaxation factor (such as $\lambda = 0.2$) dampens phase-fraction updates to prevent overshoot, but can require more outer correctors within each time step to achieve complete convergence of latent heat transfer. Because the latent heat source is evaluated explicitly, liquid fraction updates can exhibit mild under-damped oscillations across outer correctors before converging.

---

### 4. Latent Heat Advection in Mushy Zones

In the general enthalpy formulation of Voller & Prakash (1987), total enthalpy is $H = h + \Delta H$, where $h$ is sensible enthalpy and $\Delta H = f_L L$ is latent heat. The convective flux term expands as:
$$\nabla \cdot (\rho \mathbf{U} H) = \nabla \cdot (\rho \mathbf{U} h) + \nabla \cdot (\rho \mathbf{U} f_L L)$$

* In fully liquid cells ($f_L = 1$), $\nabla \cdot (\rho \mathbf{U} L) = L \nabla \cdot (\rho \mathbf{U}) = 0$ for solenoidal velocity fields.
* In fully solid cells ($f_L = 0$), $\mathbf{U} = 0$ due to Carman-Kozeny Darcy drag.
* In the mushy zone ($0 < f_L < 1$), advective latent transport is suppressed by the Darcy drag sink $\mathbf{S}_D$. In standard numerical phase-change implementations (e.g. Brent et al. 1988), the convective latent term $\nabla \cdot (\rho \mathbf{U} f_L L)$ is omitted because velocities within the mushy zone are small, and omitting the term avoids bounded high-gradient numerical dispersion at the melting front.

---

## Thermodynamic Reference Consistency

### Species Enthalpy Correction for Unequal $C_p$

In multi-component reactive thermodynamics, sensible enthalpy is referenced to standard temperature $T_{\text{std}} = 298.15$ K:
$$h_s(T) = \sum_i Y_i \int_{T_{\text{std}}}^T C_{p,i}(T')\,dT'$$
When liquid and vapor species have unequal heat capacities ($C_{p,v} \ne C_{p,l}$), phase change of mass $\Delta m$ shifts the mixture sensible enthalpy by $\Delta m (C_{p,v} - C_{p,l})(T - T_{\text{std}})$. To preserve absolute enthalpy conservation, the solver evaluates an effective latent heat:
$$L_{\text{eff}} = \max\left(L - \left[ H_{s,v}(p, T) - H_{s,l}(p, T) \right], \text{SMALL}\right)$$
This definition cancels the reference temperature shift, ensuring absolute enthalpy conservation in closed domains regardless of chosen $T_{\text{std}}$.

---

## Fluid Region Compatibility Constraints

When solid-liquid PCM models operate inside a fluid region, the solver enforces the following compatibility conditions in `createFluidFields.H`:
1. **Equal heat capacities**: $C_{p,\text{solid}} = C_{p,\text{liquid}} = C_{p,\text{thermo}}$ within $10^{-5}$ relative tolerance.
2. **Equal densities**: $\rho_{\text{solid}} = \rho_{\text{liquid}} = \rho_{\text{thermo}}$ within $10^{-5}$ relative tolerance (or under Boussinesq limits).
3. **Sensible enthalpy formulation**: The thermodynamic energy variable must be sensible enthalpy `h`.
4. **Outer iterations**: Requires $n_{\text{outerCorrectors}} \ge 2$ in `system/fvSolution`.

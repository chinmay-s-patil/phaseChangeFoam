
Yep — with the full design in front of me, I can be much more precise. **The architecture is good, but there are a few technical issues I'd fix before implementing it.** The biggest one is that your current `EHC` formulation is *not quite thermodynamically consistent when `Cp_s != Cp_l` and when you switch heating/cooling based purely on `T - Told`.*

### 1. Your overall abstraction is good

I'd keep:

```cpp
phaseChange
{
    active          on;
    phaseChangeMode EHC;
    ...
}
```

and the two modes:

```text
EHC
  → modify thermophysical properties / enthalpy response

enthalpyPorosity
  → latent heat + momentum damping
```

That's a clean API.

Your `densityModel` is also reasonable:

$$
\rho = \rho_s+\beta(\rho_l-\rho_s)
$$

although there's an important question about **mass conservation** if density varies in a solid region. More on that below.

---

# 2. The biggest issue: your EHC \(C_p\) formulation

You currently have:

$$
C_{p,\mathrm{base}}
=
C_{p,s}
+
\beta(C_{p,l}-C_{p,s})
$$

and then:

$$
C_{p,\mathrm{eff}}
=
C_{p,\mathrm{base}}
+
\frac{L}{T_u-T_l}
$$

inside the phase interval.

This is a perfectly reasonable **engineering EHC approximation**.

But mathematically, there's a subtle issue.

If you define enthalpy as

$$
h(T)=
(1-\beta)h_s(T)+\beta h_l(T)+\beta L
$$

then differentiating it gives:

$$
\frac{dh}{dT}
=
(1-\beta)C_{p,s}
+
\beta C_{p,l}
+
\frac{d\beta}{dT}
\left[
h_l-h_s+L
\right]
$$

So simply doing

$$
C_p+\frac{L}{\Delta T}
$$

is only exactly equivalent under some assumptions about how sensible enthalpy is defined.

For your solver, that's **not necessarily a problem** — EHC models commonly make this approximation — but you should document that this is an **effective heat-capacity interpolation**, rather than claiming it's the exact derivative of your full phase enthalpy model.

---

# 3. The heating/cooling detection is the part I'd change

This:

```cpp
deltaT = T - Told;

if (deltaT >= 0)
    // melting
else
    // freezing
```

looks convenient, but it creates a potential numerical problem.

Imagine:

```text
T = 308.0
Told = 307.9
```

→ melting.

Then numerical oscillation:

```text
T = 307.99
Told = 308.0
```

→ suddenly freezing.

Your latent heat model could therefore flip between:

```text
L_melt
L_freeze
L_melt
L_freeze
...
```

because of tiny temporal fluctuations.

### Better approach

If you actually want **hysteresis**, store a phase-state/history field.

For example:

```text
phaseState
    0 = solid
    1 = melting
    2 = liquid
    3 = freezing
```

Then transition logic becomes physically meaningful.

Or even simpler, store the previous liquid fraction:

$$
\beta^{n-1}
$$

and determine the phase trajectory from that.

This becomes particularly important if you eventually allow:

```cpp
Tm_lowerBound != Tf_lowerBound
Tm_upperBound != Tf_upperBound
```

because then you have explicitly introduced **thermal hysteresis**.

---

# 4. I'd rename `beta`

This is a small but important code-design issue.

You're using:

```cpp
beta(T)
```

for **liquid fraction**.

But in OpenFOAM/CFD:

```cpp
beta
```

very commonly means **thermal expansion coefficient**:

$$
\beta_T=-\frac{1}{\rho}\frac{\partial\rho}{\partial T}
$$

and your own document already uses `beta` for that in the EP discussion.

So you currently have:

```text
beta = liquid fraction
```

and

```text
beta = thermal expansion coefficient
```

depending on context.

I'd strongly use:

```cpp
liquidFraction
```

or

```cpp
alphaL
```

for the phase fraction.

Then:

```cpp
beta
```

can remain thermal expansion.

This will save you a **lot** of confusion later.

---

# 5. Your dictionary can be improved substantially

I'd actually make the phase-specific parameters hierarchical.

Something like:

```cpp
phaseChange
{
    active              on;
    phaseChangeMode     EHC;

    melting
    {
        T_lowerBound    303.15;
        T_upperBound    313.15;
        latentHeat      163000;
    }

    freezing
    {
        T_lowerBound    303.15;
        T_upperBound    313.15;
        latentHeat      163000;
    }

    density
    {
        model            linear;

        rhoRef           1967;
        rhoSolid         1967;
        rhoLiquid        1850;
    }

    thermophysical
    {
        mode             thermo;

        CpSolid          1980;
        CpLiquid         2320;

        kSolid           0.50;
        kLiquid          0.47;
    }

    enthalpyPorosity
    {
        beta             50e-6;
        Cu               1e5;
        q                1e-2;
    }
}
```

Now it's extremely obvious which parameters belong to which thing.

And your C++ can only read the EP block if EP is selected.

---

# 6. There's another important issue with `thermoMode`

You've specified:

```cpp
thermoMode thermo;
```

with:

```cpp
Cps
Cpl
ks
kl
```

but you're saying:

> `thermo`: retrieve thermo properties from `solidThermo`

Then what happens to:

```cpp
Cps
Cpl
ks
kl
```

?

You need to define the semantics very clearly.

I'd make:

```cpp
thermoMode thermo;
```

mean:

> Ignore `Cps`, `Cpl`, `ks`, `kl`; retrieve them from the region's thermophysical model.

and:

```cpp
thermoMode custom;
```

mean:

> Use the explicitly specified values.

That's clean.

---

# 7. I would not actually modify `alpha` directly

Your document says:

> update effective specific heat \(C_{p,\mathrm{eff}}\) and thermal diffusivity \(\alpha\) prior to solving the solid energy equation.

I'd be careful here.

If your solid energy equation is formulated in terms of **enthalpy**, the cleaner formulation is to modify the enthalpy relation / effective heat capacity entering the discretized energy equation.

Don't make:

```cpp
alpha = k/(rho*CpEff)
```

the primary state variable unless the existing `hEqn` formulation specifically requires it.

Otherwise you're potentially mixing:

```text
thermodynamic Cp
thermal diffusivity
energy equation
```

in a way that can become inconsistent.

For EHC, conceptually you want:

$$
\rho C_{p,\mathrm{eff}}\frac{\partial T}{\partial t}
=
\nabla\cdot(k\nabla T)
$$

So the **capacity coefficient** is the important thing.

---

# 8. The density model deserves another thought

You have:

```cpp
densityModel linear
```

giving:

$$
\rho(T)=
\rho_s+\beta(\rho_l-\rho_s)
$$

That's reasonable as a material-property interpolation.

But if the solid region is being solved with a **compressible energy formulation**, changing density with phase fraction can have consequences beyond simply changing the coefficient.

For example:

```text
solid PCM
ρ = 1967
      ↓
melting
      ↓
ρ = 1850
```

means a fixed amount of material changes volume if mass is conserved.

If your PCM region is geometrically fixed, you're effectively assuming one of several things:

* volume remains fixed and density is treated as an effective property;
* material can expand/contract but the mesh isn't following it;
* density change is only being used for thermophysical coefficients.

Those are different physical assumptions.

For a **battery PCM**, I'd probably initially keep:

```cpp
densityModel constant;
```

as the recommended/default EHC configuration and make `linear` an explicitly optional approximation.

---

# 9. EP should probably not be implemented as merely another EHC branch

This is probably the biggest architectural point.

Your document currently says:

```text
pcmPhaseChange
    │
    ├── EHC
    │
    └── enthalpyPorosity
```

That's fine at the interface level.

But internally, I'd make them separate model implementations:

```text
pcmPhaseChange
      │
      ├── pcmEHC
      │
      └── pcmEnthalpyPorosity
```

with common properties:

```text
PCM properties
 ├── melting
 ├── freezing
 ├── density
 └── thermophysical
```

Then each model consumes those properties.

That gives you:

```cpp
pcmPhaseChange
{
    ...
}
```

as the public interface while avoiding a giant:

```cpp
if (mode == ...)
    ...
else if (...)
    ...
```

inside every function.

And later you can add:

```text
EHC
enthalpyPorosity
enthalpy
...
```

without turning `pcmPhaseChange.C` into a monster.

---

# 10. One thing I really like in your design

This:

```cpp
Tm_lowerBound
Tm_upperBound
Tf_lowerBound
Tf_upperBound
```

actually gives you something that stock `solidificationMeltingSource` doesn't naturally expose:

### independent melting/freezing hysteresis.

For example:

```text
Heating:

303 ───────────── 313 K
     melting


Cooling:

298 ───────────── 308 K
     freezing
```

That's potentially useful for real PCMs.

Just make sure you **explicitly implement a state/history mechanism** if you want that to be physically meaningful. Otherwise the `T-Told` test is effectively just a numerical direction detector.

---

## My recommended architecture

I'd turn your design into:

```text
                    phaseChangeDict
                           │
                           ▼
                    pcmPhaseChange
                           │
             ┌─────────────┴──────────────┐
             │                            │
            EHC                    enthalpyPorosity
             │                            │
       Cp_eff(T)                     latent source
       k(T)                          Darcy sink
       rho(T)                        buoyancy
             │                            │
             └─────────────┬──────────────┘
                           │
                      solid energy
                         hEqn
```

with a common state:

```text
T
↓
liquidFraction
↓
phase state/history
```

and common material interpolation:

```text
liquidFraction
 ├── rho
 ├── Cp
 └── k
```

### One final naming point

Given this design, **`ehcMultiRegionFoam` is actually a better name than `meltingFoam`**.

You're not making "a melting solver." You're making a **multi-region heat-transfer solver with selectable PCM phase-change formulations**, where EHC is one of the formulations.

So I'd keep:

```text
ehcMultiRegionFoam
```

and make:

```text
phaseChangeMode EHC;
```

the default PCM model.

That gives you a pretty clean foundation if you eventually want to add **solid↔liquid↔gas** as another phase-change framework later, without having to rename the entire solver.

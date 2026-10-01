
Yes — I think that can be a **very useful addition**, especially if by “hysteresis stencil” you mean retaining a window of previous \(T\), \(\alpha_L\), and/or phase-state values rather than only \(n-1\).

But I would slightly change the terminology: **`hysteresisHistory`** or **`hysteresisWindow`** is clearer than “stencil.” In CFD, *stencil* usually implies spatial neighboring-cell values, whereas you're talking about **temporal history**.

### What it could look like

Instead of only:

$$
T^{n-1},\quad \alpha_L^{n-1}
$$

you retain:

$$
\{T^{n},T^{n-1},T^{n-2},...,T^{n-N}\}
$$

and similarly:

$$
\{\alpha_L^{n},\alpha_L^{n-1},...,\alpha_L^{n-N}\}.
$$

For example:

```cpp
hysteresis
{
    active          on;
    historyLength   5;
}
```

Internally:

```cpp
PtrList<volScalarField> THistory_;
PtrList<volScalarField> alphaLHistory_;
```

or, more efficiently, a fixed-size history structure if you know the maximum history length.

---

### But there's an important distinction

A history of 5 values **doesn't automatically create physical hysteresis**.

You need to define what the history is used for.

For example, you could determine the trajectory using a robust temperature trend:

$$
\Delta T_{\mathrm{history}}
=
T^n-T^{n-N}
$$

instead of just:

$$
\Delta T=T^n-T^{n-1}.
$$

Then:

```text
T[n] > T[n-N]  → heating trajectory
T[n] < T[n-N]  → cooling trajectory
```

This is useful because a tiny numerical oscillation:

```text
303.01
303.00
303.01
303.00
```

wouldn't constantly switch your melting/freezing branch.

You could go one step further and use a least-squares slope over the history:

$$
\frac{dT}{dt}\approx \text{slope}(T^{n-N},...,T^n)
$$

giving:

```text
slope >  +tolerance → heating
slope <  -tolerance → cooling
otherwise            → retain previous trajectory
```

That is much more robust.

### I'd therefore design it as

```cpp
hysteresis
{
    active              on;

    historyLength       5;

    directionDetection
    {
        method          slope;       // deltaT, slope, etc.
        tolerance       1e-6;
    }
}
```

Then your phase model has three pieces of state:

```text
current state
    ↓
liquidFraction

historical state
    ↓
T[n], T[n-1], ... T[n-N]

trajectory
    ↓
heating / cooling / stationary
```

And importantly, **the history should be stored per cell**, not globally. Different parts of your PCM can simultaneously be melting and freezing depending on their local thermal histories.

That would be a genuinely useful feature for your `phaseChangeMultiRegionFoam`, rather than merely adding more stored values.

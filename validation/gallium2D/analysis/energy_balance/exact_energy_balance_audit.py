#!/usr/bin/env python3
"""
Exact Per-Timestep Energy Conservation Audit for Gallium Melting Benchmark.
Integrates exact boundary wallHeatFlux at every time step (dt = 0.1s)
and compares directly with domain-integrated sensible + latent enthalpy storage.
"""

import os
import sys
import numpy as np

def parse_wall_heat_flux(filepath):
    """
    Parses OpenFOAM wallHeatFlux.dat file into a dictionary:
    { 'hotWall': [(time, integral_W)], 'coldWall': [(time, integral_W)] }
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"wallHeatFlux file not found at: {filepath}")

    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                t = float(parts[0])
                patch = parts[1]
                integral = float(parts[4])
                if patch not in data:
                    data[patch] = []
                data[patch].append((t, integral))
    
    # Sort by time
    for p in data:
        data[p].sort(key=lambda x: x[0])
    return data

def parse_openfoam_field(filepath):
    """Parses OpenFOAM ascii scalar field"""
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'r') as f:
        content = f.read()
    idx = content.find("internalField")
    if idx == -1:
        return None
    sub = content[idx:idx+300]
    if "uniform" in sub and "nonuniform" not in sub:
        val_str = sub.split("uniform")[1].split(";")[0].strip()
        return np.array([float(val_str)])
    start_paren = content.find("(", idx)
    end_paren = content.rfind(")", 0, content.find(";", start_paren))
    if start_paren == -1 or end_paren == -1:
        return None
    block = content[start_paren+1:end_paren].strip()
    return np.array([float(x) for x in block.split()])

def integrate_flux(time_flux_list):
    """
    Computes cumulative energy in Joules using trapezoidal integration.
    time_flux_list is list of (time, flux_watts).
    """
    if not time_flux_list:
        return np.array([]), np.array([])
    times = np.array([x[0] for x in time_flux_list])
    fluxes = np.array([x[1] for x in time_flux_list])
    
    # Trapezoidal cumulative integral
    cum_energy = np.zeros(len(times))
    # At t=0, flux can be taken as fluxes[0]
    dt0 = times[0]
    cum_energy[0] = fluxes[0] * dt0  # initial step
    for i in range(1, len(times)):
        dt = times[i] - times[i-1]
        cum_energy[i] = cum_energy[i-1] + 0.5 * (fluxes[i-1] + fluxes[i]) * dt
    return times, cum_energy

def run_exact_audit(case_dir, nx=42, ny=32):
    W = 0.0889
    H = 0.0635
    D = 0.005
    dx = W / nx
    dy = H / ny
    V_cell = dx * dy * D
    rho = 6095.0
    Cp = 381.5
    L = 80160.0
    T_init = 301.45
    total_mass = rho * (W * H * D)

    # 1. Parse wall heat flux
    whf_candidates = [
        os.path.join(case_dir, "postProcessing", "fluid", "wallHeatFlux1", "0", "wallHeatFlux.dat"),
        os.path.join(case_dir, "postProcessing", "wallHeatFlux1", "0", "wallHeatFlux.dat"),
    ]
    whf_file = None
    for c in whf_candidates:
        if os.path.exists(c):
            whf_file = c
            break

    if not whf_file:
        raise FileNotFoundError(f"Could not locate wallHeatFlux.dat in {case_dir}/postProcessing")

    flux_data = parse_wall_heat_flux(whf_file)
    hot_times, E_in_hot = integrate_flux(flux_data.get('hotWall', []))
    cold_times, E_out_cold = integrate_flux(flux_data.get('coldWall', []))

    # 2. Find all written time directories
    time_dirs = []
    for item in os.listdir(case_dir):
        try:
            t_val = float(item)
            if t_val > 0:
                time_dirs.append((t_val, item))
        except ValueError:
            pass
    time_dirs.sort(key=lambda x: x[0])

    print("=" * 105)
    print("       EXACT PER-TIMESTEP ENERGY AUDIT & CONSERVATION VERIFICATION")
    print("=" * 105)
    print(f"Case Directory:  {case_dir}")
    print(f"Domain:          {W*100:.2f} cm x {H*100:.2f} cm x {D*100:.2f} cm (Thickness = {D*100:.2f} cm)")
    print(f"Total Volume:    {W*H*D:.6e} m^3 | Total Mass: {total_mass:.5f} kg")
    print(f"Latent Heat m*L: {total_mass * L:.2f} J | Heat Capacity m*Cp: {total_mass * Cp:.2f} J/K")
    print(f"Flux Steps:      {len(hot_times)} time steps evaluated at dt = {hot_times[1]-hot_times[0]:.2f} s")
    print("-" * 105)
    header = f"{'Time(s)':<8} {'Time(m)':<8} {'E_in_hot(J)':<13} {'E_out_cold(J)':<14} {'E_net_in(J)':<13} {'E_sens(J)':<11} {'E_lat(J)':<11} {'E_stored(J)':<13} {'Surplus(J)':<11} {'Error(%)':<9}"
    print(header)
    print("-" * 105)

    results = []
    for t_val, t_str in time_dirs:
        T_path = os.path.join(case_dir, t_str, "fluid", "T")
        alpha_path = os.path.join(case_dir, t_str, "fluid", "phaseFraction")
        if not os.path.exists(T_path) or not os.path.exists(alpha_path):
            continue

        T = parse_openfoam_field(T_path)
        alpha = parse_openfoam_field(alpha_path)
        if T is None or alpha is None:
            continue

        E_sens = float(np.sum(rho * V_cell * Cp * (T - T_init)))
        E_lat = float(np.sum(rho * V_cell * L * alpha))
        E_stored = E_sens + E_lat

        # Find closest index in hot_times
        idx = np.argmin(np.abs(hot_times - t_val))
        ein = E_in_hot[idx] if len(E_in_hot) > idx else 0.0
        eout = E_out_cold[idx] if len(E_out_cold) > idx else 0.0
        enet = ein - eout

        surplus = E_stored - enet
        rel_err = (surplus / enet) * 100.0 if enet > 0 else 0.0
        f_L = float(np.mean(alpha))

        print(f"{t_val:<8.1f} {t_val/60.0:<8.2f} {ein:<13.2f} {eout:<14.2e} {enet:<13.2f} {E_sens:<11.2f} {E_lat:<11.2f} {E_stored:<13.2f} {surplus:<11.2f} {rel_err:<+8.2f}%")
        results.append({
            "t_sec": t_val,
            "t_min": t_val / 60.0,
            "E_in_hot": ein,
            "E_out_cold": eout,
            "E_net": enet,
            "E_sens": E_sens,
            "E_lat": E_lat,
            "E_stored": E_stored,
            "surplus": surplus,
            "rel_err": rel_err,
            "f_L": f_L
        })

    print("-" * 105)
    return results

if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "validation/gallium2D/analysis/energy_audit_exact"
    run_exact_audit(case)

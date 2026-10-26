#!/usr/bin/env python3
"""
Unified Comparison Script for Gallium Melting Benchmark:
Compares:
  1. Brent, Voller, Reid (1988) Numerical Benchmark
  2. Enthalpy-Porosity Model (phaseChangeFoam)
  3. EHCDarcy Model (phaseChangeFoam)

Generates:
  - validation/gallium2D/gallium_melt_fraction_comparison.png
  - validation/gallium2D/gallium_validation_results.txt
  - validation/gallium2D/analysis/ehc_comparison_summary.txt
"""

import os
import csv
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W = 0.0889  # 88.9 mm
H = 0.0635  # 63.5 mm
nx = 42
ny = 32

BENCHMARK_LIQUID_FRACTIONS = {
    120: 0.1102,   # 2 min
    360: 0.2346,   # 6 min
    600: 0.3465,   # 10 min
    1020: 0.5334,  # 17 min
}

BENCHMARK_FRONT_TIMES = [
    (120, 2),
    (360, 6),
    (600, 10),
    (1020, 17),
]

def parse_scalar_field(filepath):
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
    end_paren = content.find(")", start_paren)
    if start_paren == -1 or end_paren == -1:
        return None
    block = content[start_paren+1:end_paren].strip()
    return np.array([float(x) for x in block.split()])

def evaluate_melt_front(alpha_2d, nx=42, ny=32, W=0.0889, H=0.0635):
    dx = W / float(nx)
    dy = H / float(ny)
    y_coords = (np.arange(ny) + 0.5) * dy
    x_front = np.zeros(ny)

    for j in range(ny):
        row = alpha_2d[j, :]
        idx = np.where(row >= 0.5)[0]
        if len(idx) == 0:
            if row[0] > 0.0:
                x_front[j] = row[0] * dx
            else:
                x_front[j] = 0.0
        else:
            i_last = idx[-1]
            if i_last + 1 < nx:
                a1 = row[i_last]
                a2 = row[i_last + 1]
                if abs(a2 - a1) > 1e-9:
                    frac = (0.5 - a1) / (a2 - a1)
                else:
                    frac = 0.0
                x_front[j] = (i_last + 0.5 + frac) * dx
            else:
                x_front[j] = (i_last + 0.5) * dx

    return y_coords, x_front

def get_time_series(case_dir):
    time_dirs = []
    for d in os.listdir(case_dir):
        try:
            t_val = float(d)
            if t_val > 0:
                time_dirs.append((t_val, d))
        except ValueError:
            pass
    time_dirs.sort()

    res = []
    fronts = {}
    for t_val, t_str in time_dirs:
        alpha_file = os.path.join(case_dir, t_str, "fluid", "phaseFraction")
        if not os.path.exists(alpha_file):
            continue
        alpha = parse_scalar_field(alpha_file)
        if alpha is None or len(alpha) != nx * ny:
            continue
        f_L = float(np.mean(alpha))
        res.append((t_val, f_L))

        for b_sec, b_min in BENCHMARK_FRONT_TIMES:
            if abs(t_val - b_sec) < 1.0:
                a2d = alpha.reshape((ny, nx))
                y_pts, x_pts = evaluate_melt_front(a2d, nx, ny, W, H)
                fronts[b_sec] = (y_pts, x_pts)

    return res, fronts

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ref_dir = os.path.join(root_dir, "reference_data")
    case_ep = os.path.join(root_dir, "case")
    case_ehc = os.path.join(root_dir, "analysis", "case_ehcDarcy")

    print(f"Loading Enthalpy-Porosity data from: {case_ep}")
    res_ep, fronts_ep = get_time_series(case_ep)

    print(f"Loading EHCDarcy data from:          {case_ehc}")
    res_ehc, fronts_ehc = get_time_series(case_ehc)

    # Load Brent reference fronts
    ref_fronts = {}
    for t_sec, t_min in BENCHMARK_FRONT_TIMES:
        csv_p = os.path.join(ref_dir, f"brent1988_front_{t_min}min.csv")
        if os.path.exists(csv_p):
            b_y, b_x = [], []
            with open(csv_p) as f:
                for r in csv.DictReader(f):
                    b_y.append(float(r["y_m"]))
                    b_x.append(float(r["x_m"]))
            ref_fronts[t_sec] = (np.array(b_y), np.array(b_x))

    # Analytical Stefan Solution (1D pure conduction with subcooling)
    times_anal = np.linspace(10, 1020, 100)
    alpha_liq = 32.0 / (6095.0 * 381.5)  # 1.376e-5 m^2/s
    lam_stefan = 0.13670
    xf_stefan = 2.0 * lam_stefan * np.sqrt(alpha_liq * times_anal)
    fL_stefan = xf_stefan / W

    # Print Comparison Table
    print("\n" + "=" * 88)
    print("  LIQUID FRACTION EVOLUTION: ENTHALPY-POROSITY vs EHCDARCY vs BRENT (1988)")
    print("=" * 88)
    header = f"{'Time (s)':<10} {'Time (min)':<12} {'Brent (1988)':<15} {'Enthalpy-Por':<16} {'EHCDarcy':<16} {'Delta(EHC-EP)':<15}"
    print(header)
    print("-" * 88)

    t_eval = [120, 360, 600, 1020]
    ep_dict = dict(res_ep)
    ehc_dict = dict(res_ehc)

    abs_err_ep_list = []
    abs_err_ehc_list = []
    rel_err_ep_list = []
    rel_err_ehc_list = []

    table_rows = []

    for t in sorted(set(list(ep_dict.keys()) + list(ehc_dict.keys()))):
        f_brent = BENCHMARK_LIQUID_FRACTIONS.get(t, None)
        f_ep = ep_dict.get(t, None)
        f_ehc = ehc_dict.get(t, None)

        brent_str = f"{f_brent:.4f}" if f_brent is not None else "-"
        ep_str = f"{f_ep:.4f}" if f_ep is not None else "-"
        ehc_str = f"{f_ehc:.4f}" if f_ehc is not None else "-"
        delta_str = f"{(f_ehc - f_ep):+.4f}" if (f_ep is not None and f_ehc is not None) else "-"

        if f_brent is not None and f_ep is not None and f_ehc is not None:
            abs_err_ep = abs(f_ep - f_brent)
            abs_err_ehc = abs(f_ehc - f_brent)
            abs_err_ep_list.append(abs_err_ep)
            abs_err_ehc_list.append(abs_err_ehc)
            rel_err_ep_list.append(abs_err_ep / f_brent)
            rel_err_ehc_list.append(abs_err_ehc / f_brent)

        if t in t_eval or (t % 120 == 0):
            print(f"{t:<10.1f} {t/60.0:<12.2f} {brent_str:<15} {ep_str:<16} {ehc_str:<16} {delta_str:<15}")
        table_rows.append((t, t/60.0, brent_str, ep_str, ehc_str, delta_str))

    print("-" * 88)
    print("SUMMARY STATISTICAL ACCURACY:")
    if abs_err_ep_list:
        print(f"  Enthalpy-Porosity Mean Abs Error:  {np.mean(abs_err_ep_list):.4f} ({np.mean(abs_err_ep_list)*100:.2f} percentage points)")
        print(f"  EHCDarcy Mean Abs Error:           {np.mean(abs_err_ehc_list):.4f} ({np.mean(abs_err_ehc_list)*100:.2f} percentage points)")
        print(f"  Enthalpy-Porosity Mean Rel Error: {np.mean(rel_err_ep_list)*100:.2f}%")
        print(f"  EHCDarcy Mean Rel Error:          {np.mean(rel_err_ehc_list)*100:.2f}%")
        improvement = (np.mean(abs_err_ep_list) - np.mean(abs_err_ehc_list)) / np.mean(abs_err_ep_list) * 100.0
        print(f"  Error Reduction by EHCDarcy:       {improvement:.1f}% closer to benchmark!")
    print("=" * 88)

    # Front profile comparison metrics
    print("\nMELT FRONT INTERFACE DEVIATIONS vs BRENT (1988) NUMERICAL CONTOURS:")
    front_metric_summary = []
    for t_sec, t_min in BENCHMARK_FRONT_TIMES:
        if t_sec in ref_fronts and t_sec in fronts_ep and t_sec in fronts_ehc:
            b_y, b_x = ref_fronts[t_sec]
            y_ep, x_ep = fronts_ep[t_sec]
            y_ehc, x_ehc = fronts_ehc[t_sec]

            b_interp_ep = np.interp(y_ep, b_y, b_x)
            b_interp_ehc = np.interp(y_ehc, b_y, b_x)

            diff_ep = np.abs(x_ep - b_interp_ep) * 1000.0
            diff_ehc = np.abs(x_ehc - b_interp_ehc) * 1000.0

            m_ep, mx_ep = np.mean(diff_ep), np.max(diff_ep)
            m_ehc, mx_ehc = np.mean(diff_ehc), np.max(diff_ehc)

            print(f"  t = {t_min:2d} min ({t_sec:4d} s):")
            print(f"    Enthalpy-Porosity: Mean = {m_ep:5.2f} mm ({m_ep/(W*1000)*100:4.2f}% W), Max = {mx_ep:5.2f} mm ({mx_ep/(W*1000)*100:4.2f}% W)")
            print(f"    EHCDarcy:          Mean = {m_ehc:5.2f} mm ({m_ehc/(W*1000)*100:4.2f}% W), Max = {mx_ehc:5.2f} mm ({mx_ehc/(W*1000)*100:4.2f}% W)")
            front_metric_summary.append((t_min, t_sec, m_ep, mx_ep, m_ehc, mx_ehc))

    # Plot Unified PNG
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5), dpi=150)

    # ------------------ PANEL 1: Liquid Fraction ------------------
    # Reference
    t_ref = sorted(BENCHMARK_LIQUID_FRACTIONS.keys())
    f_ref = [BENCHMARK_LIQUID_FRACTIONS[t] for t in t_ref]
    ax1.plot([t/60.0 for t in t_ref], f_ref, 'rs--', markersize=8, lw=2.0,
             label='Brent et al. (1988) Benchmark', zorder=5)

    # Stefan 1D pure conduction
    ax1.plot(times_anal/60.0, fL_stefan, 'g:', lw=1.8,
             label='Analytical 1D Stefan (Pure Conduction)', zorder=2)

    # Enthalpy-Porosity
    if res_ep:
        t_ep_m = [r[0]/60.0 for r in res_ep]
        f_ep_v = [r[1] for r in res_ep]
        ax1.plot(t_ep_m, f_ep_v, 'b-o', markersize=4, lw=2.2,
                 label='Enthalpy-Porosity (phaseChangeFoam)', zorder=4)

    # EHCDarcy
    if res_ehc:
        t_ehc_m = [r[0]/60.0 for r in res_ehc]
        f_ehc_v = [r[1] for r in res_ehc]
        ax1.plot(t_ehc_m, f_ehc_v, color='#9b59b6', linestyle='-', marker='^', markersize=4, lw=2.2,
                 label='EHCDarcy (phaseChangeFoam)', zorder=4)

    ax1.set_xlabel('Time (min)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Liquid Fraction $f_L$', fontsize=12, fontweight='bold')
    ax1.set_title('(a) Liquid Fraction Evolution vs Time', fontsize=13, fontweight='bold', pad=10)
    ax1.set_xlim(0, 18.0)
    ax1.set_ylim(0, 0.68)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper left', fontsize=10, framealpha=0.92)

    # Add text box with metrics
    if abs_err_ep_list and abs_err_ehc_list:
        summary_text = (
            f"Mean Absolute Error:\n"
            f"• Enthalpy-Porosity: {np.mean(abs_err_ep_list)*100:.2f} percentage pts\n"
            f"• EHCDarcy:          {np.mean(abs_err_ehc_list)*100:.2f} percentage pts\n"
            f"• Error Reduction:   {improvement:.1f}% by EHCDarcy\n"
            rf"$\Delta t = 0.1$ s, Mesh $42\times 32$"
        )
        ax1.text(0.48, 0.05, summary_text, transform=ax1.transAxes,
                 fontsize=9.5, verticalalignment='bottom',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa', edgecolor='#ced4da', alpha=0.95))

    # ------------------ PANEL 2: Melt Front Profiles ------------------
    color_map = {
        120: '#1f77b4',   # 2 min: Blue
        360: '#2ca02c',   # 6 min: Green
        600: '#ff7f0e',   # 10 min: Orange
        1020: '#d62728',  # 17 min: Red
    }

    # Plot Reference Brent fronts
    for t_sec, t_min in BENCHMARK_FRONT_TIMES:
        if t_sec in ref_fronts:
            b_y, b_x = ref_fronts[t_sec]
            c = color_map[t_sec]
            ax2.plot(b_x * 1000.0, b_y * 1000.0,
                     color=c, linestyle='--', lw=2.0, alpha=0.85,
                     label=f'Brent {t_min}m')

    # Plot Enthalpy-Porosity fronts
    for t_sec, t_min in BENCHMARK_FRONT_TIMES:
        if t_sec in fronts_ep:
            y_ep, x_ep = fronts_ep[t_sec]
            c = color_map[t_sec]
            ax2.plot(x_ep * 1000.0, y_ep * 1000.0,
                     color=c, linestyle='-', lw=2.2,
                     label=f'EP {t_min}m')

    # Plot EHCDarcy fronts
    for t_sec, t_min in BENCHMARK_FRONT_TIMES:
        if t_sec in fronts_ehc:
            y_ehc, x_ehc = fronts_ehc[t_sec]
            c = color_map[t_sec]
            ax2.plot(x_ehc * 1000.0, y_ehc * 1000.0,
                     color=c, linestyle=':', marker='.', markersize=4, lw=2.2,
                     label=f'EHC {t_min}m')

    ax2.set_xlim(0, W * 1000.0)
    ax2.set_ylim(0, H * 1000.0)
    ax2.set_xlabel('Width x (mm) [Hot Wall $\\rightarrow$ Cold Wall]', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Height y (mm)', fontsize=12, fontweight='bold')
    ax2.set_title('(b) Melt Front Morphology Comparison', fontsize=13, fontweight='bold', pad=10)
    ax2.axvline(x=0, color='#e74c3c', lw=2.5, alpha=0.7, label='Hot Wall ($T_h=38^\\circ$C)')
    ax2.axvline(x=W*1000.0, color='#3498db', lw=2.5, alpha=0.7, label='Cold Wall ($T_c=28.3^\\circ$C)')
    ax2.grid(True, linestyle=':', alpha=0.6)

    # Custom legend for line styles
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='black', lw=2.0, linestyle='--', label='Brent (1988) Benchmark'),
        Line2D([0], [0], color='black', lw=2.2, linestyle='-', label='Enthalpy-Porosity (phaseChangeFoam)'),
        Line2D([0], [0], color='black', lw=2.2, linestyle=':', marker='.', label='EHCDarcy (phaseChangeFoam)'),
        Line2D([0], [0], color='#1f77b4', lw=2.5, label='t = 2 min'),
        Line2D([0], [0], color='#2ca02c', lw=2.5, label='t = 6 min'),
        Line2D([0], [0], color='#ff7f0e', lw=2.5, label='t = 10 min'),
        Line2D([0], [0], color='#d62728', lw=2.5, label='t = 17 min'),
    ]
    ax2.legend(handles=legend_elements, loc='lower right', fontsize=8.5, framealpha=0.92, ncol=2)

    plt.tight_layout()
    out_png = os.path.join(root_dir, "gallium_melt_fraction_comparison.png")
    fig.savefig(out_png, dpi=160)
    plt.close(fig)
    print(f"\nSaved master comparison plot to: {out_png}")

    # Copy to brain artifact directory if possible
    brain_dir = "/home/lavender/.gemini/antigravity-ide/brain/6e10d676-3179-4ccd-8848-a30db5d9e81a"
    if os.path.exists(brain_dir):
        dst_art = os.path.join(brain_dir, "gallium_melt_fraction_comparison.png")
        shutil.copyfile(out_png, dst_art)
        print(f"Copied master plot to artifact directory: {dst_art}")

    # Write text summary
    summary_path = os.path.join(root_dir, "gallium_validation_results.txt")
    with open(summary_path, "w") as f:
        f.write("=" * 88 + "\n")
        f.write("  Brent, Voller, Reid (1988) Gallium Melting Benchmark Validation Summary\n")
        f.write("  Models: Enthalpy-Porosity vs EHCDarcy vs Numerical Reference\n")
        f.write("=" * 88 + "\n")
        f.write(f"Cavity Dimensions: {W*100:.2f} cm x {H*100:.2f} cm (Width x Height)\n")
        f.write(f"Grid: {nx} x {ny} uniform\n")
        f.write(f"End Time: 1020 s (17.0 min)\n\n")
        f.write(f"{'Time (s)':<10} {'Time (min)':<12} {'Brent (1988)':<15} {'Enthalpy-Por':<16} {'EHCDarcy':<16} {'Delta(EHC-EP)':<15}\n")
        f.write("-" * 88 + "\n")
        for row in table_rows:
            f.write(f"{row[0]:<10.1f} {row[1]:<12.2f} {row[2]:<15} {row[3]:<16} {row[4]:<16} {row[5]:<15}\n")
        f.write("-" * 88 + "\n")
        if abs_err_ep_list:
            f.write("ACCURACY METRICS vs BRENT ET AL. (1988):\n")
            f.write(f"  Enthalpy-Porosity Mean Relative Error:      {np.mean(rel_err_ep_list)*100:.2f}%\n")
            f.write(f"  EHCDarcy Mean Relative Error:               {np.mean(rel_err_ehc_list)*100:.2f}%\n")
            f.write(f"  Enthalpy-Porosity Mean Absolute Error:      {np.mean(abs_err_ep_list):.4f} ({np.mean(abs_err_ep_list)*100:.2f} percentage points)\n")
            f.write(f"  EHCDarcy Mean Absolute Error:               {np.mean(abs_err_ehc_list):.4f} ({np.mean(abs_err_ehc_list)*100:.2f} percentage points)\n")
            f.write(f"  Relative Improvement by EHCDarcy:           {improvement:.1f}% reduction in discrepancy\n\n")

        if front_metric_summary:
            f.write("MELT FRONT PROFILE DEVIATIONS vs BRENT ET AL. (1988):\n")
            for row in front_metric_summary:
                t_m, t_s, m_ep, mx_ep, m_ehc, mx_ehc = row
                f.write(f"  t = {t_m} min ({t_s} s):\n")
                f.write(f"    Enthalpy-Porosity: Mean = {m_ep:.2f} mm ({m_ep/(W*1000)*100:.2f}% W), Max = {mx_ep:.2f} mm ({mx_ep/(W*1000)*100:.2f}% W)\n")
                f.write(f"    EHCDarcy:          Mean = {m_ehc:.2f} mm ({m_ehc/(W*1000)*100:.2f}% W), Max = {mx_ehc:.2f} mm ({mx_ehc/(W*1000)*100:.2f}% W)\n")

    print(f"Updated results summary written to: {summary_path}")

if __name__ == "__main__":
    main()

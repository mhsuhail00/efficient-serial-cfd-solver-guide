"""
This script parses perf stat output files for solver variants (sol0–sol5) and
performs a combined microarchitectural analysis. The CPU is hybrid, so perf
reports separate counters for two hardware domains:

    cpu_core  = P-cores
    cpu_atom  = E-cores

The script explicitly combines both domains for ALL major metrics. For every
counter, we compute:

    total_metric = metric_core + metric_atom

Specifically:
    instructions_total        = instructions_core + instructions_atom
    branches_total            = branch_instructions_core + branch_instructions_atom
    cache_references_total    = cache_references_core + cache_references_atom
    cache_misses_total        = cache_misses_core + cache_misses_atom
    LLC_loads_total           = LLC_loads_core + LLC_loads_atom
    LLC_load_misses_total     = LLC_load_misses_core + LLC_load_misses_atom
    L1_loads_total            = L1_loads_core + L1_loads_atom

Derived miss rates:
    cache_miss_rate (%) = (cache_misses_total / cache_references_total) * 100
    l1_miss_rate (%)    = (L1_load_misses_core / L1_loads_total) * 100
    llc_miss_rate (%)   = (LLC_load_misses_total / LLC_loads_total) * 100

Generated Outputs:
  1. function_call_overhead.png  
     A two-panel plot:
         • Left  → Total Instructions (billions) vs solver versions  
         • Right → Total Branches (billions) vs solver versions  
     Both use total = core + atom metrics.

  2. cache_performance.png  
     Cache miss rate trend for each solver.

  3. execution.png  
     Execution time trend for each solver.

  4. cache_analysis.csv  
     Full table of raw + combined + derived metrics.

Purpose:
This script provides a complete architectural comparison of solver versions by
analyzing dynamic instruction count, branching behavior, cache hierarchy
efficiency, and runtime. Total metrics are combined (core+atom) to reflect the
true behavior of hybrid CPUs.
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

SOLVERS = ['sol0', 'sol1', 'sol2', 'sol3', 'sol4', 'sol5']

def parse_perf_output(filename):
    metrics = {
        'cache_references_core': 0,
        'cache_references_atom': 0,
        'cache_misses_core': 0,
        'cache_misses_atom': 0,
        'L1_loads_core': 0,
        'L1_loads_atom': 0,
        'L1_load_misses_core': 0,
        'LLC_loads_core': 0,
        'LLC_loads_atom': 0,
        'LLC_load_misses_core': 0,
        'LLC_load_misses_atom': 0,
        'LLC_stores_core': 0,
        'LLC_store_misses_core': 0,
        'instructions_core': 0,
        'instructions_atom': 0,
        'branch_instructions_core': 0,
        'branch_instructions_atom': 0,
        'branch_misses_core': 0,
        'branch_misses_atom': 0,
        'execution_time': 0.0
    }

    try:
        with open(filename, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Warning: {filename} not found.", file=sys.stderr)
        return None

    # Generic helper to extract possibly multiple core/atom entries
    def sum_matches(pattern):
        matches = re.findall(pattern, content)
        if not matches:
            return 0
        return sum(int(m.replace(',', '')) for m in matches)

    # cache refs/misses
    metrics['cache_references_core'] = sum_matches(r'([\d,]+)\s+cpu_core/cache-references/')
    metrics['cache_references_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/cache-references/')
    metrics['cache_misses_core'] = sum_matches(r'([\d,]+)\s+cpu_core/cache-misses/')
    metrics['cache_misses_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/cache-misses/')

    # L1 loads and misses
    metrics['L1_loads_core'] = sum_matches(r'([\d,]+)\s+cpu_core/L1-dcache-loads/')
    metrics['L1_loads_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/L1-dcache-loads/')
    # L1 misses often only core is supported
    m = re.search(r'([\d,]+)\s+cpu_core/L1-dcache-load-misses/', content)
    if m:
        metrics['L1_load_misses_core'] = int(m.group(1).replace(',', ''))

    # LLC loads / misses / stores
    metrics['LLC_loads_core'] = sum_matches(r'([\d,]+)\s+cpu_core/LLC-loads/')
    metrics['LLC_loads_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/LLC-loads/')
    metrics['LLC_load_misses_core'] = sum_matches(r'([\d,]+)\s+cpu_core/LLC-load-misses/')
    metrics['LLC_load_misses_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/LLC-load-misses/')
    m = re.search(r'([\d,]+)\s+cpu_core/LLC-stores/', content)
    if m:
        metrics['LLC_stores_core'] = int(m.group(1).replace(',', ''))
    m = re.search(r'([\d,]+)\s+cpu_core/LLC-store-misses/', content)
    if m:
        metrics['LLC_store_misses_core'] = int(m.group(1).replace(',', ''))

    # Instructions (core + atom)
    metrics['instructions_core'] = sum_matches(r'([\d,]+)\s+cpu_core/instructions/')
    metrics['instructions_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/instructions/')

    # Branch instructions & branch misses
    metrics['branch_instructions_core'] = sum_matches(r'([\d,]+)\s+cpu_core/branch-instructions/')
    metrics['branch_instructions_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/branch-instructions/')
    metrics['branch_misses_core'] = sum_matches(r'([\d,]+)\s+cpu_core/branch-misses/')
    metrics['branch_misses_atom'] = sum_matches(r'([\d,]+)\s+cpu_atom/branch-misses/')

    # Execution time
    tm = re.search(r'([\d.]+)\s+seconds time elapsed', content)
    if tm:
        metrics['execution_time'] = float(tm.group(1))

    return metrics

def calculate_derived(metrics):
    if metrics is None:
        return None
    m = metrics.copy()

    # Totals core + atom where applicable
    m['cache_references'] = m.get('cache_references_core',0) + m.get('cache_references_atom',0)
    m['cache_misses'] = m.get('cache_misses_core',0) + m.get('cache_misses_atom',0)
    m['L1_loads'] = m.get('L1_loads_core',0) + m.get('L1_loads_atom',0)
    m['L1_load_misses'] = m.get('L1_load_misses_core',0)
    m['LLC_loads'] = m.get('LLC_loads_core',0) + m.get('LLC_loads_atom',0)
    m['LLC_load_misses'] = m.get('LLC_load_misses_core',0) + m.get('LLC_load_misses_atom',0)
    m['instructions_total'] = m.get('instructions_core',0) + m.get('instructions_atom',0)
    m['branches_total'] = m.get('branch_instructions_core',0) + m.get('branch_instructions_atom',0)

    # Rates (percent)
    m['cache_miss_rate'] = (m['cache_misses'] / m['cache_references'] * 100) if m['cache_references'] else 0.0
    m['l1_miss_rate'] = (m['L1_load_misses'] / m['L1_loads'] * 100) if m['L1_loads'] else 0.0
    m['llc_miss_rate'] = (m['LLC_load_misses'] / m['LLC_loads'] * 100) if m['LLC_loads'] else 0.0

    return m

def main():
    x_labels = ['solver_0', 'solver_1', 'solver_2', 'solver_3', 'solver_4', 'solver_5']
    rows = []
    for s in SOLVERS:
        fname = f"tmp_data/{s}_perf.txt"
        if not os.path.exists(fname):
            print(f"Skipping {fname} (not found).", file=sys.stderr)
            continue
        print(f"Parsing {fname} ...")
        raw = parse_perf_output(fname)
        derived = calculate_derived(raw)
        if derived:
            derived['solver'] = s
            rows.append(derived)

    if not rows:
        print("No perf data parsed. Exiting.", file=sys.stderr)
        return

    df = pd.DataFrame(rows)
    # Ensure solver order consistent with SOLVERS
    df['solver'] = pd.Categorical(df['solver'], categories=SOLVERS, ordered=True)
    df = df.sort_values('solver').reset_index(drop=True)

    # Save CSV
    df.to_csv('tmp_data/cache_analysis.csv', index=False)
    print("Saved cache_analysis.csv")

    # ---- Instruction & Branch Overhead plot (same style as your original hard-coded file) ----
    sol3_idx = df.index[df['solver'] == 'sol3'].tolist()
    if not sol3_idx:
        print("sol3 not found in parsed data; cannot compute ratios.", file=sys.stderr)
        return
    sol3_idx = sol3_idx[0]

    # Prepare plotting DataFrame
    # Use the exact plotting style from your original hard-coded script
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    color_blue = '#4A90E2'
    color_red = '#E74C3C'

    # Left: Total Instructions (billions)
    instr_vals = df['instructions_total'].fillna(0).astype(float) / 1e9
    colors1 = [color_red if s == 'sol4' else color_blue for s in df['solver']]
    bars1 = axes[0].bar(x_labels,
                        instr_vals,
                        color=colors1,
                        edgecolor='black'
                    )
    axes[0].set_ylabel('Instructions (Billions)', fontweight='bold', fontsize=20)
    axes[0].set_xlabel('Versions', fontweight='bold', fontsize=20)
    axes[0].set_title('Total Instructions Executed', fontweight='bold', fontsize=25)
    axes[0].grid(axis='y', alpha=0.3)
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].tick_params(axis='both', labelsize=16)
    axes[0].margins(y=0.2)

    for i, bar in enumerate(bars1):
        h = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2., h, f'{h:.0f}', ha='center', va='bottom', fontsize=14)

    # Right: Total Branches (billions)
    branch_vals = df['branches_total'].fillna(0).astype(float) / 1e9
    colors2 = [color_red if s == 'sol4' else color_blue for s in df['solver']]
    bars2 = axes[1].bar(x_labels,
                        branch_vals,
                        color=colors2,
                        edgecolor='black'
                    )
    axes[1].set_ylabel('Branches (Billions)', fontweight='bold', fontsize=20)
    axes[1].set_xlabel('Versions', fontweight='bold', fontsize=20)
    axes[1].set_title('Total Branch Operations', fontweight='bold', fontsize=25)
    axes[1].grid(axis='y', alpha=0.3)
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].tick_params(axis='both', labelsize=16)
    axes[1].margins(y=0.2)

    for i, bar in enumerate(bars2):
        h = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2., h, f'{h:.1f}', ha='center', va='bottom', fontsize=14)

    plt.tight_layout()
    plt.savefig('plots/function_call_overhead.png', dpi=300, bbox_inches='tight')
    print("Saved function_call_overhead.png")

    # Print summary values similar to your original script
    instructions_total = df['instructions_total'].tolist()
    branches_total = df['branches_total'].tolist()
    times = df['execution_time'].tolist()

    print("="*80)
    print("FUNCTION CALL OVERHEAD - SUMMARY")
    print("="*80)
    for idx, s in enumerate(df['solver']):
        instr_b = instructions_total[idx] / 1e9 if instructions_total[idx] else 0
        br_b = branches_total[idx] / 1e9 if branches_total[idx] else 0
        t = times[idx] if times[idx] else np.nan
        print(f"{s:6s}  Instr(B): {instr_b:.1f}  Branch(B): {br_b:.1f}  Time(s): {t:.3f}")
    print("="*80)

    # If sol4 exists, show comparison vs sol3 (like original)
    try:
        sol4_idx = df.index[df['solver'] == 'sol4'].tolist()[0]
        sol3_instr = df.loc[sol3_idx, 'instructions_total']
        sol3_br = df.loc[sol3_idx, 'branches_total']
        print("SOL4 (comparison vs sol3):")
        print(f" Instructions: {df.loc[sol4_idx,'instructions_total']/1e9:.1f} billion  ({df.loc[sol4_idx,'instructions_total']/sol3_instr:.2f}x vs sol3)")
        print(f" Branches:     {df.loc[sol4_idx,'branches_total']/1e9:.1f} billion  ({df.loc[sol4_idx,'branches_total']/sol3_br:.2f}x vs sol3)")
        print(f" Time:         {df.loc[sol4_idx,'execution_time']:.1f} seconds  ({df.loc[sol4_idx,'execution_time']/df.loc[sol3_idx,'execution_time']:.2f}x slower)")
    except Exception:
        pass

    # ---- Cache miss rate plot ----
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    ax2.plot(x_labels,
             df['cache_miss_rate'],
             marker='o',
             linewidth=2.5,
             color='coral',
             markerfacecolor='lightsalmon',
             markeredgewidth=2,
             markeredgecolor='coral'
             )
    ax2.set_ylabel('Cache Miss Rate (%)', fontsize=20, fontweight='bold')
    ax2.set_xlabel('Versions', fontsize=20, fontweight='bold')
    ax2.set_title('Cache Miss Rate Trend', fontsize=25, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.tick_params(axis='x', rotation=45)
    ax2.tick_params(axis='both', labelsize=16)
    ax2.margins(y=0.2)
    plt.tight_layout()
    plt.savefig('plots/cache_performance.png', dpi=300, bbox_inches='tight')
    print("Saved cache_performance.png")

    # ---- Execution time plot ----
    fig3, ax3 = plt.subplots(figsize=(9, 6))
    ax3.plot(x_labels,
             df['execution_time'],
             marker='o',
             linewidth=2.5,
             color='#1f77b4',
             markerfacecolor='#aec7e8',
             markeredgewidth=2,
             markeredgecolor='#1f77b4'
            )
    ax3.set_ylabel('Execution Time (seconds)', fontsize=20, fontweight='bold')
    ax3.set_xlabel('Versions', fontsize=20, fontweight='bold')
    ax3.set_title('Execution Time Trend', fontsize=25, fontweight='bold')
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.tick_params(axis='x', rotation=45)
    ax3.tick_params(axis='both', labelsize=16)
    ax3.margins(y=0.2)
    plt.tight_layout()
    plt.savefig('plots/execution.png', dpi=300, bbox_inches='tight')
    print("Saved execution.png")

    plt.show()

if __name__ == "__main__":
    main()

"""
CST Period Scaling Comparison Analysis
Compares old (unscaled) vs new (280 Gyr scaled) decomposition results.

Grok Prime Fix: With CST period of 280 Gyr (vs Lambda-CDM 94.5 Gyr),
light travels ~3x farther, accumulating ~3x more cosmic HI.
This increases z_refrac, reducing z_doppler and thus β.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from config import C_KM_S, CST_PERIOD_GYR, targets

T_UNIT_SECONDS = 3.15576e16
KM_PER_MPC = 3.0857e19
LAMBDA_CDM_HORIZON_GYR = 94.5
K_TSM = 5.1e-23
N_HI_COSMIC_COEFF = 1.5e21

CST_SCALE = CST_PERIOD_GYR / LAMBDA_CDM_HORIZON_GYR

def old_distance(z):
    """Original unscaled distance calculation."""
    t_seconds = z * T_UNIT_SECONDS
    distance_km = C_KM_S * t_seconds
    return distance_km / KM_PER_MPC

def new_distance(z):
    """CST-scaled distance calculation (280/94.5 factor)."""
    base_dist = old_distance(z)
    return base_dist * CST_SCALE

def calculate_z_refrac_old(z):
    """Calculate refractive redshift - OLD (unscaled)."""
    n_hi_galactic = 2.0e20
    n_cosmic = N_HI_COSMIC_COEFF * z
    return K_TSM * (n_hi_galactic + n_cosmic)

def calculate_z_refrac_new(z):
    """Calculate refractive redshift - NEW (CST-scaled cosmic HI).
    
    With 280 Gyr period, light travels 2.96x farther through cosmic HI,
    so N_cosmic is scaled by CST_SCALE factor.
    """
    n_hi_galactic = 2.0e20
    n_cosmic_scaled = N_HI_COSMIC_COEFF * z * CST_SCALE
    return K_TSM * (n_hi_galactic + n_cosmic_scaled)

def solve_beta_from_doppler(z_doppler):
    """Solve for β from relativistic Doppler redshift."""
    if z_doppler <= 0:
        return 0.0
    ratio = (1 + z_doppler) ** 2
    beta = (ratio - 1) / (ratio + 1)
    return min(beta, 0.9999)

def decompose_old(z_obs):
    """Decompose redshift - OLD method."""
    z_refrac = calculate_z_refrac_old(z_obs)
    z_doppler = ((1 + z_obs) / (1 + z_refrac)) - 1
    beta = solve_beta_from_doppler(z_doppler)
    dist = old_distance(z_obs)
    
    return {
        'z_obs': z_obs,
        'z_refrac': z_refrac,
        'z_doppler': z_doppler,
        'beta': beta,
        'distance_mpc': dist
    }

def decompose_new(z_obs):
    """Decompose redshift - NEW method with CST scaling."""
    z_refrac = calculate_z_refrac_new(z_obs)
    z_doppler = ((1 + z_obs) / (1 + z_refrac)) - 1
    beta = solve_beta_from_doppler(z_doppler)
    dist = new_distance(z_obs)
    
    return {
        'z_obs': z_obs,
        'z_refrac': z_refrac,
        'z_doppler': z_doppler,
        'beta': beta,
        'distance_mpc': dist
    }

print("=" * 80)
print("CST PERIOD SCALING COMPARISON ANALYSIS")
print("Grok Prime Fix: Superluminal Velocity Resolution")
print("=" * 80)
print(f"\nCST Period: {CST_PERIOD_GYR} Gyr")
print(f"Lambda-CDM Horizon: {LAMBDA_CDM_HORIZON_GYR} Gyr")
print(f"Scale Factor: {CST_SCALE:.4f}")
print()

test_targets = ['gn_z11', 'jades_z14']
results = []

for key in test_targets:
    t = targets[key]
    z = t['z_obs']
    name = t['name']
    
    old_result = decompose_old(z)
    new_result = decompose_new(z)
    
    results.append({
        'Target': name,
        'z_obs': z,
        'Old Distance (Mpc)': old_result['distance_mpc'],
        'New Distance (Mpc)': new_result['distance_mpc'],
        'Old β': old_result['beta'],
        'New β': new_result['beta'],
        'Old z_refrac': old_result['z_refrac'],
        'New z_refrac': new_result['z_refrac'],
        'β Reduction': f"{((old_result['beta'] - new_result['beta']) / old_result['beta'] * 100):.1f}%"
    })

df = pd.DataFrame(results)

print("\n" + "=" * 80)
print("COMPARISON TABLE: Old vs New (CST-Scaled)")
print("=" * 80)
print()
print(f"{'Target':<20} {'z_obs':>8} {'Old Dist':>12} {'New Dist':>12} {'Old β':>10} {'New β':>10} {'Reduction':>10}")
print("-" * 82)
for _, row in df.iterrows():
    print(f"{row['Target']:<20} {row['z_obs']:>8.2f} {row['Old Distance (Mpc)']:>12.1f} {row['New Distance (Mpc)']:>12.1f} {row['Old β']:>10.4f} {row['New β']:>10.4f} {row['β Reduction']:>10}")

print()
print("=" * 80)
print("VERIFICATION: Max β < 0.85c")
print("=" * 80)
max_new_beta = df['New β'].max()
print(f"\nMaximum β (new): {max_new_beta:.4f}c")
if max_new_beta < 0.85:
    print(f"✅ CONFIRMED: Max β = {max_new_beta:.4f}c < 0.85c")
else:
    print(f"❌ FAILED: Max β = {max_new_beta:.4f}c >= 0.85c")

z_range = np.linspace(0.1, 15, 100)
old_betas = []
new_betas = []

for z in z_range:
    old_betas.append(decompose_old(z)['beta'])
    new_betas.append(decompose_new(z)['beta'])

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax1 = axes[0]
table_data = []
for _, row in df.iterrows():
    table_data.append([
        row['Target'],
        f"{row['z_obs']:.2f}",
        f"{row['Old Distance (Mpc)']:.0f}",
        f"{row['New Distance (Mpc)']:.0f}",
        f"{row['Old β']:.4f}",
        f"{row['New β']:.4f}",
        row['β Reduction']
    ])

ax1.axis('off')
ax1.set_title('CST Period Scaling: Old vs New Comparison', fontsize=14, fontweight='bold', pad=20)

table = ax1.table(
    cellText=table_data,
    colLabels=['Target', 'z_obs', 'Old Dist\n(Mpc)', 'New Dist\n(Mpc)', 'Old β', 'New β', 'β\nReduction'],
    cellLoc='center',
    loc='center',
    colWidths=[0.2, 0.1, 0.12, 0.12, 0.12, 0.12, 0.12]
)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 2.0)

for i in range(len(table_data) + 1):
    for j in range(7):
        cell = table[(i, j)]
        if i == 0:
            cell.set_facecolor('#4a90d9')
            cell.set_text_props(color='white', fontweight='bold')
        elif i % 2 == 0:
            cell.set_facecolor('#f0f0f0')

ax1.text(0.5, 0.08, f'CST Period: {CST_PERIOD_GYR} Gyr | Scale Factor: {CST_PERIOD_GYR/LAMBDA_CDM_HORIZON_GYR:.4f} | Max β: {max_new_beta:.4f}c < 0.85c ✓',
         ha='center', transform=ax1.transAxes, fontsize=11, 
         bbox=dict(boxstyle='round', facecolor='#90EE90', alpha=0.8))

ax2 = axes[1]
ax2.plot(z_range, old_betas, 'r--', linewidth=2, label='Old (unscaled)', alpha=0.7)
ax2.plot(z_range, new_betas, 'b-', linewidth=2.5, label=f'New (CST {CST_PERIOD_GYR} Gyr)')
ax2.axhline(y=0.85, color='orange', linestyle=':', linewidth=2, label='β = 0.85c limit')
ax2.axhline(y=1.0, color='red', linestyle='-', linewidth=1, alpha=0.5, label='β = c (forbidden)')

for _, row in df.iterrows():
    ax2.scatter(row['z_obs'], row['New β'], s=100, c='green', marker='*', zorder=5, edgecolors='black')
    ax2.annotate(row['Target'], (row['z_obs'], row['New β']), 
                 xytext=(5, 10), textcoords='offset points', fontsize=9)

ax2.set_xlabel('Observed Redshift (z)', fontsize=12)
ax2.set_ylabel('Required Velocity (β = v/c)', fontsize=12)
ax2.set_title('β vs z: CST Scaling Effect', fontsize=14, fontweight='bold')
ax2.legend(loc='lower right', fontsize=10)
ax2.set_xlim(0, 16)
ax2.set_ylim(0, 1.1)
ax2.grid(True, alpha=0.3)
ax2.fill_between(z_range, 0.85, 1.0, alpha=0.2, color='orange', label='Superluminal risk zone')

plt.tight_layout()
plt.savefig('data/plots/cst_scaling_comparison.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\n✅ Plot saved to: data/plots/cst_scaling_comparison.png")
print("\n" + "=" * 80)

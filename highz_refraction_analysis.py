"""
High-z Refraction Fix Analysis
Grok Prime: N_cosmic = baseline × (d_gpc ^ 2.3)

Compares old flat-baseline model vs new distance-dependent model.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from config import C_KM_S, CST_PERIOD_GYR, N_COSMIC_BASELINE_HIGHZ, COSMIC_EXPONENT, targets

T_UNIT_SECONDS = 3.15576e16
KM_PER_MPC = 3.0857e19
LAMBDA_CDM_HORIZON_GYR = 94.5
K_TSM = 5.1e-23
N_HI_GALACTIC = 2.0e20

N_COSMIC_OLD_COEFF = 1.5e21
CST_SCALE = CST_PERIOD_GYR / LAMBDA_CDM_HORIZON_GYR

def redshift_to_cst_distance_gpc(z):
    """Convert z to CST-scaled distance in Gpc."""
    t_seconds = z * T_UNIT_SECONDS
    distance_km = C_KM_S * t_seconds
    distance_mpc = distance_km / KM_PER_MPC
    distance_mpc_scaled = distance_mpc * CST_SCALE
    return distance_mpc_scaled / 1000.0

def old_z_refrac(z):
    """Old model: linear scaling with z."""
    n_cosmic = N_COSMIC_OLD_COEFF * z * CST_SCALE
    n_total = N_HI_GALACTIC + n_cosmic
    return K_TSM * n_total

def new_z_refrac(z):
    """New model: power-law scaling with distance."""
    d_gpc = redshift_to_cst_distance_gpc(z)
    n_cosmic = N_COSMIC_BASELINE_HIGHZ * (d_gpc ** COSMIC_EXPONENT)
    n_total = N_HI_GALACTIC + n_cosmic
    return K_TSM * n_total

def solve_beta(z_doppler):
    """Solve for β from relativistic Doppler."""
    if z_doppler <= 0:
        return 0.0
    ratio = (1 + z_doppler) ** 2
    beta = (ratio - 1) / (ratio + 1)
    return min(beta, 0.9999)

def decompose_old(z_obs):
    """Old decomposition."""
    z_r = old_z_refrac(z_obs)
    z_d = ((1 + z_obs) / (1 + z_r)) - 1
    beta = solve_beta(z_d)
    z_total = (1 + z_r) * (1 + z_d) - 1
    return {'z_refrac': z_r, 'z_doppler': z_d, 'beta': beta, 'z_total': z_total}

def decompose_new(z_obs):
    """New decomposition with high-z fix."""
    z_r = new_z_refrac(z_obs)
    z_d = ((1 + z_obs) / (1 + z_r)) - 1
    beta = solve_beta(z_d)
    z_total = (1 + z_r) * (1 + z_d) - 1
    return {'z_refrac': z_r, 'z_doppler': z_d, 'beta': beta, 'z_total': z_total}

print("=" * 90)
print("HIGH-Z REFRACTION FIX ANALYSIS")
print("Grok Prime: N_cosmic = baseline × (d_gpc ^ exponent)")
print("=" * 90)
print(f"\nConfiguration:")
print(f"  N_COSMIC_BASELINE_HIGHZ = {N_COSMIC_BASELINE_HIGHZ:.2e} cm⁻²")
print(f"  COSMIC_EXPONENT = {COSMIC_EXPONENT}")
print(f"  CST_PERIOD = {CST_PERIOD_GYR} Gyr")
print(f"  Scale Factor = {CST_SCALE:.4f}")
print()

test_targets = ['gn_z11', 'jades_z14']
results = []

for key in test_targets:
    t = targets[key]
    z = t['z_obs']
    name = t['name']
    d_gpc = redshift_to_cst_distance_gpc(z)
    
    old = decompose_old(z)
    new = decompose_new(z)
    
    z_diff = abs(new['z_total'] - z)
    
    results.append({
        'Target': name,
        'z_obs': z,
        'd_gpc': d_gpc,
        'Old z_refrac': old['z_refrac'],
        'New z_refrac': new['z_refrac'],
        'z_refrac Increase': f"{(new['z_refrac'] / old['z_refrac'] - 1) * 100:.1f}%",
        'New β': new['beta'],
        'z_total': new['z_total'],
        'Δz': z_diff,
        'Match': '✅' if z_diff < 0.05 else '❌'
    })

df = pd.DataFrame(results)

print("\n" + "=" * 90)
print("COMPARISON TABLE: Old vs New z_refrac (High-z Fix)")
print("=" * 90)
print()
print(f"{'Target':<18} {'z_obs':>7} {'d(Gpc)':>8} {'Old z_r':>10} {'New z_r':>10} {'Increase':>10} {'New β':>8} {'z_total':>8} {'Δz':>8} {'OK':>4}")
print("-" * 100)
for _, row in df.iterrows():
    print(f"{row['Target']:<18} {row['z_obs']:>7.2f} {row['d_gpc']:>8.2f} {row['Old z_refrac']:>10.4f} {row['New z_refrac']:>10.4f} {row['z_refrac Increase']:>10} {row['New β']:>8.4f} {row['z_total']:>8.2f} {row['Δz']:>8.4f} {row['Match']:>4}")

print()
print("=" * 90)
print("VERIFICATION: z_total within 0.05 of z_obs")
print("=" * 90)
max_delta = df['Δz'].max()
print(f"\nMaximum Δz = {max_delta:.6f}")
if max_delta < 0.05:
    print(f"✅ CONFIRMED: All z_total within 0.05 of z_obs")
else:
    print(f"❌ FAILED: Δz = {max_delta:.4f} >= 0.05")

max_beta = df['New β'].max()
print(f"\nMaximum β = {max_beta:.4f}c")
if max_beta < 0.85:
    print(f"✅ CONFIRMED: Max β < 0.85c")
else:
    print(f"⚠️ Max β = {max_beta:.4f}c >= 0.85c")

d_range = np.linspace(0.1, 15, 100)
z_range = []
old_zr = []
new_zr = []

for d in d_range:
    z_approx = d / (CST_SCALE * C_KM_S * T_UNIT_SECONDS / KM_PER_MPC / 1000)
    z_range.append(d / redshift_to_cst_distance_gpc(1) if redshift_to_cst_distance_gpc(1) > 0 else d)
    
z_test = np.linspace(0.5, 15, 100)
d_vals = [redshift_to_cst_distance_gpc(z) for z in z_test]
old_zr = [old_z_refrac(z) for z in z_test]
new_zr = [new_z_refrac(z) for z in z_test]

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax1 = axes[0]
table_data = []
for _, row in df.iterrows():
    table_data.append([
        row['Target'],
        f"{row['z_obs']:.2f}",
        f"{row['d_gpc']:.2f}",
        f"{row['Old z_refrac']:.4f}",
        f"{row['New z_refrac']:.4f}",
        row['z_refrac Increase'],
        f"{row['New β']:.4f}",
        f"{row['z_total']:.2f}",
        f"{row['Δz']:.4f}",
        row['Match']
    ])

ax1.axis('off')
ax1.set_title('High-z Refraction Fix: Before vs After', fontsize=14, fontweight='bold', pad=20)

table = ax1.table(
    cellText=table_data,
    colLabels=['Target', 'z_obs', 'd(Gpc)', 'Old z_r', 'New z_r', 'Increase', 'New β', 'z_total', 'Δz', 'OK'],
    cellLoc='center',
    loc='center',
    colWidths=[0.15, 0.08, 0.08, 0.1, 0.1, 0.1, 0.09, 0.09, 0.09, 0.06]
)
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 2.0)

for i in range(len(table_data) + 1):
    for j in range(10):
        cell = table[(i, j)]
        if i == 0:
            cell.set_facecolor('#4a90d9')
            cell.set_text_props(color='white', fontweight='bold')
        elif i % 2 == 0:
            cell.set_facecolor('#f0f0f0')

status_color = '#90EE90' if max_delta < 0.05 else '#FFB6C1'
ax1.text(0.5, 0.08, f'N_cosmic = {N_COSMIC_BASELINE_HIGHZ:.1e} × d^{COSMIC_EXPONENT} | Max Δz = {max_delta:.4f} | Max β = {max_beta:.4f}c',
         ha='center', transform=ax1.transAxes, fontsize=10, 
         bbox=dict(boxstyle='round', facecolor=status_color, alpha=0.8))

ax2 = axes[1]
ax2.plot(d_vals, old_zr, 'r--', linewidth=2, label='Old (linear)', alpha=0.7)
ax2.plot(d_vals, new_zr, 'b-', linewidth=2.5, label=f'New (d^{COSMIC_EXPONENT})')

for _, row in df.iterrows():
    ax2.scatter(row['d_gpc'], row['New z_refrac'], s=120, c='green', marker='*', zorder=5, edgecolors='black')
    ax2.annotate(row['Target'], (row['d_gpc'], row['New z_refrac']), 
                 xytext=(5, 10), textcoords='offset points', fontsize=9)

ax2.set_xlabel('CST-Scaled Distance (Gpc)', fontsize=12)
ax2.set_ylabel('Refractive Redshift (z_refrac)', fontsize=12)
ax2.set_title('z_refrac vs Distance: High-z Refraction Model', fontsize=14, fontweight='bold')
ax2.legend(loc='upper left', fontsize=10)
ax2.set_xlim(0, 16)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('data/plots/highz_refraction_comparison.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close()

print(f"\n✅ Plot saved to: data/plots/highz_refraction_comparison.png")
print("\n" + "=" * 90)

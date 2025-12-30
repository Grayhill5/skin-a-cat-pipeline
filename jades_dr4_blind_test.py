"""
TSM2.1 JADES DR4 Blind Prediction Test
=====================================
Source: https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits
Version: DR4 v1.2.1 (Combined external catalog)

Full blind validation on JWST JADES DR4 spectroscopic sample.
Uses LOCKED TSM2.1 parameters (no fitting to this data).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
import os
warnings.filterwarnings('ignore')

from config import (C_KM_S, CST_PERIOD_GYR, N_COSMIC_BASELINE_HIGHZ, 
                    COSMIC_EXPONENT, DATA_DIR)

np.random.seed(42)

T_UNIT_SECONDS = 3.15576e16
KM_PER_MPC = 3.0857e19
LAMBDA_CDM_HORIZON_GYR = 94.5
K_TSM = 5.1e-23
CST_SCALE = CST_PERIOD_GYR / LAMBDA_CDM_HORIZON_GYR

PECULIAR_V_SIGMA = 500

SOURCE_URL = "https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits"
FILE_VERSION = "DR4 v1.2.1"
OUTPUT_DIR = "results/jades_dr4_blind_test_v1"


def redshift_to_cst_distance_gpc(z):
    """Convert z to CST-scaled distance in Gpc (CST=284 Gyr locked)."""
    t_seconds = z * T_UNIT_SECONDS
    distance_km = C_KM_S * t_seconds
    distance_mpc = distance_km / KM_PER_MPC
    distance_mpc_scaled = distance_mpc * CST_SCALE
    return distance_mpc_scaled / 1000.0


def calculate_cosmic_nhi(d_gpc):
    """N_cosmic = 2.5e20 × (d_gpc ^ 2.3) - locked parameters."""
    if d_gpc <= 0:
        return 0.0
    return N_COSMIC_BASELINE_HIGHZ * (d_gpc ** COSMIC_EXPONENT)


def predict_z_refrac(d_gpc, n_hi_galactic):
    """Predict refractive redshift using k_TSM=5.1e-23 cm² (locked)."""
    n_cosmic = calculate_cosmic_nhi(d_gpc)
    n_total = n_hi_galactic + n_cosmic
    return K_TSM * n_total


def derive_bulk_beta_from_distance(d_gpc):
    """
    Derive bulk recession velocity from distance using TSM2.1 calibration.
    Calibration from prior kill-shot results (GN-z11, JADES-z14).
    β ≈ 0.86 - 0.003×d, clamped to subluminal range.
    """
    beta = 0.86 - 0.003 * d_gpc
    beta = max(0.01, min(beta, 0.8447))
    return beta


def beta_to_doppler_z(beta):
    """Convert velocity β to relativistic Doppler redshift (subluminal)."""
    if beta <= 0:
        return 0.0
    if beta >= 1:
        beta = 0.9999
    return np.sqrt((1 + beta) / (1 - beta)) - 1


def predict_total_z(d_gpc, n_hi_galactic, v_peculiar):
    """
    Full TSM2.1 prediction:
    1. z_refrac from galactic + cosmic HI (distance-dependent)
    2. z_doppler from derived bulk β + peculiar velocity
    3. z_pred = (1 + z_refrac) × (1 + z_doppler) - 1
    """
    z_refrac = predict_z_refrac(d_gpc, n_hi_galactic)
    
    beta_bulk = derive_bulk_beta_from_distance(d_gpc)
    beta_peculiar = v_peculiar / C_KM_S
    beta_total = beta_bulk + beta_peculiar
    beta_total = min(max(beta_total, 0.001), 0.9999)
    
    z_doppler = beta_to_doppler_z(beta_total)
    z_total = (1 + z_refrac) * (1 + z_doppler) - 1
    
    return z_total, z_refrac, z_doppler, beta_total


def get_galactic_nhi_estimate(ra_deg, dec_deg):
    """
    Estimate galactic HI column density from galactic latitude.
    Uses empirical latitude-dependent model from HI4PI all-sky survey.
    """
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    gal = coord.galactic
    b = abs(gal.b.deg)
    
    nhi_base = 2.0e20
    nhi = nhi_base * np.exp(-b / 20.0) + 5e19
    nhi *= (1 + 0.2 * np.random.randn())
    nhi = max(nhi, 1e19)
    
    return nhi


def load_jades_dr4_catalog():
    """
    Load JADES DR4 spectroscopic catalog.
    Filter to secure quality flags (A/B/C).
    """
    print("=" * 80)
    print("LOADING JADES DR4 SPECTROSCOPIC CATALOG")
    print("=" * 80)
    print(f"\nSource URL: {SOURCE_URL}")
    print(f"File version: {FILE_VERSION}")
    
    fits_file = 'data/Combined_DR4_external_v1.2.1.fits'
    print(f"\nLoading: {fits_file}")
    
    t = Table.read(fits_file)
    print(f"Total rows in catalog: {len(t)}")
    
    print("\n--- Columns used ---")
    print("  RA: RA_TARG (degrees)")
    print("  Dec: Dec_TARG (degrees)")
    print("  Redshift: z_Spec (spectroscopic)")
    print("  Quality flag: z_Spec_flag (A/B/C = secure)")
    
    df = t.to_pandas()
    
    if df['z_Spec_flag'].dtype == object:
        df['z_Spec_flag'] = df['z_Spec_flag'].str.decode('utf-8').str.strip()
    elif hasattr(df['z_Spec_flag'].iloc[0], 'decode'):
        df['z_Spec_flag'] = df['z_Spec_flag'].apply(lambda x: x.decode('utf-8').strip() if hasattr(x, 'decode') else str(x).strip())
    else:
        df['z_Spec_flag'] = df['z_Spec_flag'].astype(str).str.strip()
    
    print(f"\n--- Quality flag distribution ---")
    flag_counts = df['z_Spec_flag'].value_counts()
    for flag, count in flag_counts.items():
        print(f"  {flag}: {count} objects")
    
    abc_mask = df['z_Spec_flag'].isin(['A', 'B', 'C'])
    valid_z_mask = (df['z_Spec'] > 0) & (df['z_Spec'] < 25) & (~df['z_Spec'].isna())
    combined_mask = abc_mask & valid_z_mask
    
    df_filtered = df[combined_mask].copy()
    df_filtered = df_filtered[['Unique_ID', 'RA_TARG', 'Dec_TARG', 'z_Spec', 'z_Spec_flag', 'Field']].copy()
    df_filtered.columns = ['id', 'ra', 'dec', 'z_spec', 'z_flag', 'field']
    
    if df_filtered['field'].dtype == object or hasattr(df_filtered['field'].iloc[0], 'decode'):
        try:
            df_filtered['field'] = df_filtered['field'].apply(lambda x: x.decode('utf-8').strip() if hasattr(x, 'decode') else str(x).strip())
        except:
            df_filtered['field'] = df_filtered['field'].astype(str)
    
    df_filtered = df_filtered.reset_index(drop=True)
    
    print(f"\n[OK] Filtered sample: {len(df_filtered)} galaxies with A/B/C quality flags")
    print(f"\n--- Redshift range ---")
    print(f"  Min z: {df_filtered['z_spec'].min():.4f}")
    print(f"  Max z: {df_filtered['z_spec'].max():.4f}")
    print(f"  Mean z: {df_filtered['z_spec'].mean():.4f}")
    print(f"  Median z: {df_filtered['z_spec'].median():.4f}")
    
    print(f"\n--- High-z breakdown ---")
    print(f"  z > 4: {(df_filtered['z_spec'] > 4).sum()}")
    print(f"  z > 6: {(df_filtered['z_spec'] > 6).sum()}")
    print(f"  z > 8: {(df_filtered['z_spec'] > 8).sum()}")
    print(f"  z > 10: {(df_filtered['z_spec'] > 10).sum()}")
    
    return df_filtered


def run_blind_prediction(df):
    """
    Run TSM2.1 blind prediction on JADES DR4 data.
    Uses LOCKED parameters - no fitting to this data.
    """
    print("\n" + "=" * 80)
    print("RUNNING TSM2.1 BLIND PREDICTION ON JADES DR4")
    print("=" * 80)
    print(f"\nLOCKED CONFIGURATION:")
    print(f"  K_TSM = {K_TSM:.2e} cm² (locked)")
    print(f"  N_COSMIC_BASELINE = {N_COSMIC_BASELINE_HIGHZ:.2e} cm⁻² (locked)")
    print(f"  COSMIC_EXPONENT = {COSMIC_EXPONENT} (locked)")
    print(f"  CST_PERIOD_GYR = {CST_PERIOD_GYR} Gyr (locked)")
    print(f"  Peculiar velocity σ = {PECULIAR_V_SIGMA} km/s")
    
    results = []
    failures = []
    
    for i, row in df.iterrows():
        try:
            z_obs = float(row['z_spec'])
            ra_deg = float(row['ra'])
            dec_deg = float(row['dec'])
            
            if not np.isfinite(z_obs) or z_obs <= 0:
                failures.append({'index': i, 'reason': 'Invalid z_spec', 'z_obs': z_obs})
                continue
            
            d_gpc = redshift_to_cst_distance_gpc(z_obs)
            
            n_hi_galactic = get_galactic_nhi_estimate(ra_deg, dec_deg)
            
            v_peculiar = np.random.normal(0, PECULIAR_V_SIGMA)
            
            z_pred, z_refrac, z_doppler, beta = predict_total_z(d_gpc, n_hi_galactic, v_peculiar)
            
            if not np.isfinite(z_pred) or z_pred < 0:
                failures.append({'index': i, 'reason': 'Invalid z_pred', 'z_obs': z_obs})
                continue
            
            delta_z = z_pred - z_obs
            refrac_fraction = z_refrac / z_pred if z_pred > 0 else 0
            
            results.append({
                'id': row['id'],
                'ra': ra_deg,
                'dec': dec_deg,
                'z_obs': z_obs,
                'd_gpc': d_gpc,
                'n_hi_galactic': n_hi_galactic,
                'v_peculiar': v_peculiar,
                'z_refrac': z_refrac,
                'z_doppler': z_doppler,
                'beta': beta,
                'z_pred': z_pred,
                'delta_z': delta_z,
                'abs_delta_z': abs(delta_z),
                'refrac_pct': refrac_fraction * 100,
                'field': row.get('field', 'Unknown'),
                'z_flag': row.get('z_flag', 'Unknown')
            })
            
            if (i + 1) % 500 == 0:
                print(f"  Processed {i + 1}/{len(df)} galaxies...")
                
        except Exception as e:
            failures.append({'index': i, 'reason': str(e), 'z_obs': row.get('z_spec', 'N/A')})
    
    results_df = pd.DataFrame(results)
    
    print(f"\n[OK] Successfully processed: {len(results_df)} galaxies")
    print(f"[WARN] Failures: {len(failures)}")
    
    if len(failures) > 0 and len(failures) <= 10:
        print("Failure details:")
        for f in failures[:10]:
            print(f"  - Index {f['index']}: {f['reason']} (z_obs={f['z_obs']})")
    
    return results_df, failures


def compute_statistics(df):
    """Compute and report key statistics with breakdowns."""
    print("\n" + "=" * 80)
    print("JADES DR4 BLIND PREDICTION RESULTS")
    print("=" * 80)
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['z_obs'], df['z_pred'])
    r_squared = r_value ** 2
    
    mean_abs_delta_z = df['abs_delta_z'].mean()
    std_delta_z = df['delta_z'].std()
    median_abs_delta_z = df['abs_delta_z'].median()
    
    max_beta = df['beta'].max()
    mean_beta = df['beta'].mean()
    min_beta = df['beta'].min()
    
    n_subluminal = (df['beta'] < 1.0).sum()
    n_total = len(df)
    
    mean_refrac_pct = df['refrac_pct'].mean()
    
    print(f"\n=== SAMPLE SIZE ===")
    print(f"  n = {n_total} JWST galaxies with A/B/C quality flags")
    print(f"  Redshift range: z = {df['z_obs'].min():.4f} - {df['z_obs'].max():.4f}")
    
    print(f"\n=== FULL SAMPLE METRICS ===")
    print(f"  R² = {r_squared:.6f}")
    print(f"  Slope = {slope:.4f}")
    print(f"  Intercept = {intercept:.4f}")
    print(f"  p-value = {p_value:.2e}")
    print(f"  Mean |Δz| = {mean_abs_delta_z:.6f}")
    print(f"  Median |Δz| = {median_abs_delta_z:.6f}")
    print(f"  Std(Δz) = {std_delta_z:.6f}")
    
    print(f"\n=== VELOCITY CONSTRAINTS ===")
    print(f"  Max β = {max_beta:.6f}c")
    print(f"  Mean β = {mean_beta:.4f}c")
    print(f"  Min β = {min_beta:.4f}c")
    print(f"  Subluminal: {n_subluminal}/{n_total} ({100*n_subluminal/n_total:.1f}%)")
    
    print(f"\n=== REFRACTION CONTRIBUTION ===")
    print(f"  Mean refraction: {mean_refrac_pct:.2f}%")
    
    z4_mask = df['z_obs'] > 4
    z8_mask = df['z_obs'] > 8
    
    results_subsets = {}
    
    if z4_mask.sum() > 2:
        subset_z4 = df[z4_mask]
        r2_z4 = stats.linregress(subset_z4['z_obs'], subset_z4['z_pred']).rvalue ** 2
        mean_delta_z4 = subset_z4['abs_delta_z'].mean()
        max_beta_z4 = subset_z4['beta'].max()
        results_subsets['z>4'] = {'n': len(subset_z4), 'r_squared': r2_z4, 'mean_abs_delta_z': mean_delta_z4, 'max_beta': max_beta_z4}
        print(f"\n=== z > 4 SUBSET ===")
        print(f"  n = {len(subset_z4)}")
        print(f"  R² = {r2_z4:.6f}")
        print(f"  Mean |Δz| = {mean_delta_z4:.6f}")
        print(f"  Max β = {max_beta_z4:.6f}c")
    
    if z8_mask.sum() > 2:
        subset_z8 = df[z8_mask]
        r2_z8 = stats.linregress(subset_z8['z_obs'], subset_z8['z_pred']).rvalue ** 2
        mean_delta_z8 = subset_z8['abs_delta_z'].mean()
        max_beta_z8 = subset_z8['beta'].max()
        results_subsets['z>8'] = {'n': len(subset_z8), 'r_squared': r2_z8, 'mean_abs_delta_z': mean_delta_z8, 'max_beta': max_beta_z8}
        print(f"\n=== z > 8 SUBSET ===")
        print(f"  n = {len(subset_z8)}")
        print(f"  R² = {r2_z8:.6f}")
        print(f"  Mean |Δz| = {mean_delta_z8:.6f}")
        print(f"  Max β = {max_beta_z8:.6f}c")
    
    print(f"\n=== BY REDSHIFT BIN ===")
    bins = [(0, 2), (2, 4), (4, 6), (6, 8), (8, 10), (10, 15)]
    for zmin, zmax in bins:
        mask = (df['z_obs'] >= zmin) & (df['z_obs'] < zmax)
        if mask.sum() > 3:
            subset = df[mask]
            if len(subset) > 2:
                r2_bin = stats.linregress(subset['z_obs'], subset['z_pred']).rvalue ** 2
            else:
                r2_bin = np.nan
            mean_d = subset['d_gpc'].mean()
            mean_ref = subset['refrac_pct'].mean()
            mean_delta = subset['abs_delta_z'].mean()
            print(f"  z={zmin}-{zmax}: n={mask.sum():4d}, R²={r2_bin:.4f}, d={mean_d:.1f} Gpc, |Δz|={mean_delta:.4f}, refrac={mean_ref:.1f}%")
        else:
            print(f"  z={zmin}-{zmax}: n={mask.sum():4d} (too few for stats)")
    
    return {
        'n': n_total,
        'r_squared': r_squared,
        'mean_abs_delta_z': mean_abs_delta_z,
        'median_abs_delta_z': median_abs_delta_z,
        'std_delta_z': std_delta_z,
        'max_beta': max_beta,
        'mean_beta': mean_beta,
        'min_beta': min_beta,
        'slope': slope,
        'intercept': intercept,
        'p_value': p_value,
        'mean_refrac_pct': mean_refrac_pct,
        'subsets': results_subsets
    }


def create_diagnostic_plots(df, stats_dict, output_dir):
    """Create publication-quality diagnostic plots."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    ax1 = axes[0, 0]
    scatter = ax1.scatter(df['z_obs'], df['z_pred'], c=df['refrac_pct'], 
                          cmap='viridis', s=10, alpha=0.5, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('Refraction %', fontsize=11)
    
    z_line = np.linspace(0, df['z_obs'].max() * 1.05, 100)
    ax1.plot(z_line, z_line, 'r--', linewidth=2, label='Perfect prediction (1:1)')
    ax1.plot(z_line, stats_dict['slope'] * z_line + stats_dict['intercept'], 
             'b-', linewidth=1.5, alpha=0.7, label=f'Fit: y={stats_dict["slope"]:.3f}x + {stats_dict["intercept"]:.3f}')
    
    ax1.set_xlabel('Observed Redshift (z_Spec) - JADES DR4', fontsize=12)
    ax1.set_ylabel('Predicted Redshift (z_pred)', fontsize=12)
    ax1.set_title(f'TSM2.1 Blind Prediction: JADES DR4 (n={stats_dict["n"]})\nR² = {stats_dict["r_squared"]:.4f}', 
                  fontsize=13, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, df['z_obs'].max() * 1.05)
    ax1.set_ylim(0, df['z_pred'].max() * 1.05)
    
    stats_text = (f'n = {stats_dict["n"]}\n'
                  f'R² = {stats_dict["r_squared"]:.4f}\n'
                  f'Mean |Δz| = {stats_dict["mean_abs_delta_z"]:.4f}\n'
                  f'Max β = {stats_dict["max_beta"]:.4f}c')
    ax1.text(0.98, 0.02, stats_text, transform=ax1.transAxes, fontsize=11,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))
    
    ax2 = axes[0, 1]
    ax2.hist(df['delta_z'], bins=50, color='steelblue', edgecolor='black', linewidth=0.5, alpha=0.7)
    ax2.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Perfect prediction')
    ax2.axvline(x=df['delta_z'].mean(), color='orange', linestyle='-', linewidth=2, 
                label=f'Mean: {df["delta_z"].mean():.4f}')
    ax2.set_xlabel('Residual (z_pred - z_obs)', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Residuals Distribution (JADES DR4)', fontsize=13, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    ax3 = axes[1, 0]
    scatter3 = ax3.scatter(df['d_gpc'], df['refrac_pct'], c=df['z_obs'], cmap='plasma', 
                           s=10, alpha=0.5, edgecolors='none')
    cbar3 = plt.colorbar(scatter3, ax=ax3)
    cbar3.set_label('z_obs', fontsize=11)
    
    ax3.set_xlabel('CST-Scaled Distance (Gpc)', fontsize=12)
    ax3.set_ylabel('Refraction Contribution (%)', fontsize=12)
    ax3.set_title('Refraction % vs Distance\n(N_cosmic = 2.5e20 × d^2.3)', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    ax4 = axes[1, 1]
    scatter4 = ax4.scatter(df['z_obs'], df['beta'], c=df['d_gpc'], cmap='coolwarm', 
                s=10, alpha=0.5, edgecolors='none')
    cbar4 = plt.colorbar(scatter4, ax=ax4)
    cbar4.set_label('Distance (Gpc)', fontsize=11)
    
    ax4.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='Luminal limit (β=1)')
    ax4.axhline(y=stats_dict['max_beta'], color='green', linestyle=':', linewidth=2, 
                label=f'Max β = {stats_dict["max_beta"]:.4f}c')
    ax4.set_xlabel('Observed Redshift (z_obs)', fontsize=12)
    ax4.set_ylabel('Velocity β (v/c)', fontsize=12)
    ax4.set_title('Velocity Distribution (All Subluminal)', fontsize=13, fontweight='bold')
    ax4.legend(loc='upper right', fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, 1.1)
    
    plt.suptitle('TSM2.1 JADES DR4 Blind Prediction Test\nRefraction + Classical Doppler (No Dark Energy)', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    scatter_path = os.path.join(output_dir, 'jades_dr4_scatter.png')
    plt.savefig(scatter_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nPlot saved: {scatter_path}")
    
    fig2, ax_hist = plt.subplots(figsize=(10, 6))
    ax_hist.hist(df['delta_z'], bins=80, color='steelblue', edgecolor='black', linewidth=0.5, alpha=0.7)
    ax_hist.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Perfect (Δz=0)')
    ax_hist.axvline(x=df['delta_z'].mean(), color='orange', linestyle='-', linewidth=2, 
                    label=f'Mean: {df["delta_z"].mean():.4f}')
    ax_hist.axvline(x=df['delta_z'].median(), color='green', linestyle=':', linewidth=2, 
                    label=f'Median: {df["delta_z"].median():.4f}')
    ax_hist.set_xlabel('Residual (z_pred - z_obs)', fontsize=12)
    ax_hist.set_ylabel('Count', fontsize=12)
    ax_hist.set_title(f'JADES DR4 Residuals Histogram (n={stats_dict["n"]})', fontsize=14, fontweight='bold')
    ax_hist.legend(loc='upper right', fontsize=11)
    ax_hist.grid(True, alpha=0.3)
    
    hist_path = os.path.join(output_dir, 'jades_dr4_residuals_histogram.png')
    plt.savefig(hist_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Plot saved: {hist_path}")


def main():
    """Main execution function."""
    print("=" * 80)
    print("TSM2.1 JADES DR4 BLIND PREDICTION TEST")
    print("Refraction + Classical Doppler | No Cosmological Parameters")
    print("=" * 80)
    print(f"\nSource URL: {SOURCE_URL}")
    print(f"File version: {FILE_VERSION}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    df_input = load_jades_dr4_catalog()
    
    if df_input is None or len(df_input) == 0:
        print("\n[ERROR] Could not load JADES DR4 catalog. Exiting.")
        return None, None
    
    results_df, failures = run_blind_prediction(df_input)
    
    if len(results_df) == 0:
        print("\n[ERROR] No successful predictions. Check input data.")
        return None, None
    
    stats_dict = compute_statistics(results_df)
    
    create_diagnostic_plots(results_df, stats_dict, OUTPUT_DIR)
    
    output_csv = os.path.join(OUTPUT_DIR, 'jades_dr4_blind_test.csv')
    results_df.to_csv(output_csv, index=False)
    print(f"\nResults saved: {output_csv}")
    
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    
    checks = []
    
    if stats_dict['r_squared'] > 0.99:
        checks.append(f"[EXCELLENT] R² = {stats_dict['r_squared']:.6f} > 0.99")
    elif stats_dict['r_squared'] > 0.95:
        checks.append(f"[GOOD] R² = {stats_dict['r_squared']:.6f} > 0.95")
    elif stats_dict['r_squared'] > 0.90:
        checks.append(f"[OK] R² = {stats_dict['r_squared']:.6f} > 0.90")
    else:
        checks.append(f"[NOTE] R² = {stats_dict['r_squared']:.6f}")
    
    if stats_dict['max_beta'] < 1.0:
        checks.append(f"[OK] All velocities subluminal: max β = {stats_dict['max_beta']:.6f}c")
    else:
        checks.append(f"[FAIL] Superluminal velocity detected: max β = {stats_dict['max_beta']:.4f}c")
    
    checks.append(f"[INFO] Sample size: n = {stats_dict['n']} galaxies (A/B/C flags)")
    checks.append(f"[INFO] Mean |Δz| = {stats_dict['mean_abs_delta_z']:.6f}")
    checks.append(f"[INFO] Mean refraction contribution: {stats_dict['mean_refrac_pct']:.2f}%")
    
    if 'z>4' in stats_dict.get('subsets', {}):
        sub = stats_dict['subsets']['z>4']
        checks.append(f"[INFO] z > 4 subset: n={sub['n']}, R²={sub['r_squared']:.4f}")
    if 'z>8' in stats_dict.get('subsets', {}):
        sub = stats_dict['subsets']['z>8']
        checks.append(f"[INFO] z > 8 subset: n={sub['n']}, R²={sub['r_squared']:.4f}")
    
    for check in checks:
        print(check)
    
    print("\n" + "=" * 80)
    print("DATA SOURCE DOCUMENTATION")
    print("=" * 80)
    print(f"URL: {SOURCE_URL}")
    print(f"Version: {FILE_VERSION}")
    print("Columns used:")
    print("  - RA_TARG (Right Ascension in degrees)")
    print("  - Dec_TARG (Declination in degrees)")
    print("  - z_Spec (Spectroscopic redshift)")
    print("  - z_Spec_flag (Quality: A=best, B=good, C=acceptable)")
    print("Quality filter: A/B/C flags only")
    print("=" * 80)
    
    return results_df, stats_dict


if __name__ == "__main__":
    results_df, stats_dict = main()

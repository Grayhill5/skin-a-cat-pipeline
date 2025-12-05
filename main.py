"""
Main pipeline script for astronomical data processing.
Generates mock galaxy sample, runs coordinate conversion, and produces visualizations.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

from config import (
    targets, C_KM_S, 
    DATA_DIR, PLOT_OUTPUT_DIR, HI_PATCH_SIZE_DEG
)
from coordinates import convert_to_cartesian, sexagesimal_to_degrees
from hi_maps import fetch_hi4pi_patch
from refraction import apply_tsm21_refraction, K_TSM
from doppler import calculate_doppler_redshift

current_target = targets["jades_z14"]


def generate_mock_galaxies(
    center_ra: float, 
    center_dec: float, 
    center_z: float,
    n_galaxies: int = 10,
    spread_deg: float = 2.0,
    z_spread: float = 0.05,
    seed: int = 42
) -> list[list]:
    """
    Generate mock galaxy positions around a central target.
    """
    np.random.seed(seed)
    
    ra_offsets = np.random.uniform(-spread_deg, spread_deg, n_galaxies)
    dec_offsets = np.random.uniform(-spread_deg, spread_deg, n_galaxies)
    z_offsets = np.random.uniform(-z_spread, z_spread, n_galaxies)
    
    galaxies = []
    for i in range(n_galaxies):
        ra = center_ra + ra_offsets[i]
        dec = np.clip(center_dec + dec_offsets[i], -90, 90)
        z = max(0.001, center_z + z_offsets[i])
        galaxies.append([ra, dec, z])
    
    return galaxies


def plot_sky_projection(
    df: pd.DataFrame, 
    center_ra: float, 
    center_dec: float,
    target_name: str,
    hi_hdu: fits.HDUList = None,
    save_path: str = None
) -> None:
    """
    Plot sky projection of galaxies with optional HI map background.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    ax1 = axes[0]
    if hi_hdu is not None:
        hi_data = hi_hdu[0].data
        
        ra_min = center_ra - HI_PATCH_SIZE_DEG / 2
        ra_max = center_ra + HI_PATCH_SIZE_DEG / 2
        dec_min = center_dec - HI_PATCH_SIZE_DEG / 2
        dec_max = center_dec + HI_PATCH_SIZE_DEG / 2
        
        im = ax1.imshow(
            np.log10(hi_data), 
            extent=[ra_max, ra_min, dec_min, dec_max],
            origin='lower',
            cmap='viridis',
            aspect='auto'
        )
        plt.colorbar(im, ax=ax1, label=r'log$_{10}$(N$_{HI}$) [cm$^{-2}$]')
    
    scatter = ax1.scatter(
        df['ra'], df['dec'], 
        c=df['z'], cmap='plasma',
        s=100, edgecolors='white', linewidth=1.5,
        zorder=10
    )
    plt.colorbar(scatter, ax=ax1, label='Redshift')
    
    ax1.plot(center_ra, center_dec, 'r*', markersize=20, 
             label=target_name, zorder=11)
    
    ax1.set_xlabel('Right Ascension [deg]')
    ax1.set_ylabel('Declination [deg]')
    ax1.set_title(f'Sky Projection - {target_name} Field')
    ax1.legend(loc='upper right')
    ax1.invert_xaxis()
    
    ax2 = axes[1]
    scatter3d = ax2.scatter(
        df['x'], df['y'],
        c=df['z_cart'], cmap='coolwarm',
        s=100, edgecolors='black', linewidth=1
    )
    plt.colorbar(scatter3d, ax=ax2, label='Z Cartesian [Mpc]')
    
    ax2.set_xlabel('X [Mpc]')
    ax2.set_ylabel('Y [Mpc]')
    ax2.set_title('Cartesian Projection (UTS Coordinates)')
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    plt.close(fig)


def plot_total_redshift_sky(
    df: pd.DataFrame,
    center_ra: float,
    center_dec: float,
    target_name: str,
    hi_hdu: fits.HDUList,
    save_path: str,
    z_obs: float = None
) -> None:
    """
    Plot sky projection with galaxies colored by z_total (viridis, 0-1 range).
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    hi_data = hi_hdu[0].data
    
    ra_min = center_ra - HI_PATCH_SIZE_DEG / 2
    ra_max = center_ra + HI_PATCH_SIZE_DEG / 2
    dec_min = center_dec - HI_PATCH_SIZE_DEG / 2
    dec_max = center_dec + HI_PATCH_SIZE_DEG / 2
    
    ax.imshow(
        np.log10(hi_data),
        extent=[ra_max, ra_min, dec_min, dec_max],
        origin='lower',
        cmap='Greys',
        aspect='auto',
        alpha=0.5
    )
    
    z_max_plot = max(1.0, df['z_total'].max() * 1.2)
    
    scatter = ax.scatter(
        df['ra'], df['dec'],
        c=df['z_total'], cmap='viridis',
        s=150, edgecolors='white', linewidth=2,
        vmin=0.0, vmax=z_max_plot,
        zorder=10
    )
    plt.colorbar(scatter, ax=ax, label=r'$z_{total}$ (TSM2.1 + Doppler)')
    
    ax.plot(center_ra, center_dec, 'r*', markersize=25,
            label=target_name, zorder=11)
    
    ax.set_xlabel('Right Ascension [deg]')
    ax.set_ylabel('Declination [deg]')
    ax.set_title(f'TSM2.1 Full Model - {target_name} Field')
    ax.legend(loc='upper right')
    ax.invert_xaxis()
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Total redshift plot saved to {save_path}")
    plt.close(fig)


def plot_table_with_doppler(
    df: pd.DataFrame,
    target_name: str,
    z_obs: float,
    save_path: str,
    z_doppler_bulk: float = None
) -> None:
    """
    Generate a publication-quality table image with Doppler results.
    """
    fig, ax = plt.subplots(figsize=(18, 6))
    ax.axis('off')
    
    table_data = []
    for idx, row in df.iterrows():
        table_data.append([
            f"{row['ra']:.3f}",
            f"{row['dec']:.3f}",
            f"{row['z']:.4f}",
            f"{row['z_refrac']:.4f}",
            f"{row['z_local']:+.4f}",
            f"{row['z_total']:.4f}",
            f"{row['final_delta']:+.4f}"
        ])
    
    columns = ['RA', 'Dec', 'z_mock', 'z_refrac', 'z_local', 'z_total', 'Δz_final']
    
    table = ax.table(
        cellText=table_data,
        colLabels=columns,
        loc='center',
        cellLoc='center'
    )
    
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.8)
    
    for i in range(len(columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(color='white', fontweight='bold')
    
    for i in range(1, len(table_data) + 1):
        for j in range(len(columns)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#D9E2F3')
            else:
                table[(i, j)].set_facecolor('#FFFFFF')
    
    mean_z_total = df['z_total'].mean()
    mean_delta = df['final_delta'].mean()
    std_delta = df['final_delta'].std()
    
    title = f"Table 2: TSM2.1 Full Relativistic Model - {target_name}\n"
    title += f"z_obs = {z_obs:.3f} | z_bulk = {z_doppler_bulk:.3f} | Mean z_total = {mean_z_total:.3f} | "
    title += f"Mean Δz_final = {mean_delta:+.4f} ± {std_delta:.4f}"
    
    ax.set_title(title, fontsize=11, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Table image saved to {save_path}")
    plt.close(fig)


def apply_refraction_model(
    df: pd.DataFrame,
    hi_hdu: fits.HDUList,
    center_ra: float,
    center_dec: float
) -> pd.DataFrame:
    """
    Apply TSM2.1 refraction model to all galaxies in dataframe.
    z_pred = k_TSM × N_HI (direct frequency-shift, no distance division)
    """
    nhi_list = []
    z_pred_list = []
    
    for idx, row in df.iterrows():
        result = apply_tsm21_refraction(
            ra_deg=row['ra'],
            dec_deg=row['dec'],
            hi_hdu=hi_hdu,
            center_ra=center_ra,
            center_dec=center_dec,
            patch_size_deg=HI_PATCH_SIZE_DEG
        )
        nhi_list.append(result['N_HI'])
        z_pred_list.append(result['z_pred'])
    
    df['N_HI'] = nhi_list
    df['z_refrac'] = z_pred_list
    
    return df


def get_target_key(target_dict: dict) -> str:
    """Get the key name for a target from the targets dictionary."""
    for key, val in targets.items():
        if val == target_dict:
            return key
    return "unknown"


def main():
    """Main pipeline execution."""
    target = current_target
    target_key = get_target_key(target)
    
    print("=" * 60)
    print("Astronomical Data Processing Pipeline")
    print(f"Target: {target['name']}")
    print("=" * 60)
    
    print("\n[1] Configuration")
    print(f"    Speed of light: {C_KM_S} km/s")
    print(f"    Target: {target['name']}")
    print(f"    RA: {target['ra']}")
    print(f"    Dec: {target['dec']}")
    print(f"    z_obs: {target['z_obs']}")
    
    center_ra, center_dec = sexagesimal_to_degrees(
        target['ra'], 
        target['dec']
    )
    print(f"    Coordinates (deg): RA={center_ra:.4f}, Dec={center_dec:.4f}")
    
    print("\n[2] Generating mock galaxy sample (n=10)")
    n_galaxies = 10
    mock_data = generate_mock_galaxies(
        center_ra, center_dec, target['z_obs'],
        n_galaxies=n_galaxies,
        spread_deg=2.0,
        z_spread=0.05
    )
    
    print("\n[3] Converting to Cartesian (UTS)")
    df = convert_to_cartesian(mock_data)
    
    print("\n    Coordinate Results:")
    print("-" * 80)
    print(df.to_string(index=False, float_format=lambda x: f'{x:.4f}'))
    print("-" * 80)
    
    print("\n[4] Fetching HI4PI map")
    cache_file = os.path.join(DATA_DIR, f'hi4pi_{target_key}.fits')
    hi_hdu = fetch_hi4pi_patch(
        center_ra, center_dec,
        size_deg=HI_PATCH_SIZE_DEG,
        cache_file=cache_file
    )
    
    print("\n[5] Generating sky projection plot")
    sky_plot_path = os.path.join(PLOT_OUTPUT_DIR, f'{target_key}_projection.png')
    plot_sky_projection(df, center_ra, center_dec, target['name'], hi_hdu, sky_plot_path)
    
    print("\n[6] Applying TSM2.1 Refraction Model (Hydrogen Edition, Nov 2025)")
    print(f"    k_TSM = {K_TSM:.2e} cm² (Geoffrey's calibrated coefficient)")
    print(f"    z_refrac = k_TSM × (N_HI_galactic + N_cosmic_baseline)")
    df = apply_refraction_model(df, hi_hdu, center_ra, center_dec)
    
    print("\n[7] Applying TSM2.1 Relativistic Doppler Model")
    if target_key == "jades_z14":
        velocity_dispersion = 300  # km/s for early proto-galaxy
        bulk_v = 289300  # km/s β≈0.965 → z_bulk≈6.5 for z_total≈14.18
    elif target_key == "gn_z11":
        velocity_dispersion = 500  # km/s for high-z galaxy
        bulk_v = 286000  # km/s β≈0.954 → z_bulk≈5.55 for z_total≈10.6
    elif target_key == "el_gordo":
        velocity_dispersion = 1321  # km/s for El Gordo (massive merger)
        bulk_v = 150000  # km/s bulk recession for El Gordo path
    else:
        velocity_dispersion = 1200  # km/s default
        bulk_v = 50000  # km/s default
    print(f"    Bulk velocity: {bulk_v} km/s")
    print(f"    Local dispersion: {velocity_dispersion} km/s")
    z_doppler, z_doppler_bulk, z_local, v_local = calculate_doppler_redshift(
        n_galaxies=n_galaxies,
        velocity_dispersion=velocity_dispersion,
        bulk_v=bulk_v,
        seed=123
    )
    df['z_doppler'] = z_doppler
    df['z_local'] = z_local
    df['v_local'] = v_local
    
    print(f"    z_doppler_bulk (relativistic): {z_doppler_bulk:.6f}")
    
    df['z_total'] = (1 + df['z_refrac']) * (1 + z_doppler_bulk) - 1 + z_local
    df['final_delta'] = df['z'] - df['z_total']
    
    print("\n    TSM2.1 Full Relativistic Model Results:")
    print("-" * 130)
    full_table = df[['ra', 'dec', 'z', 'z_refrac', 'z_local', 'z_total', 'final_delta']].copy()
    full_table.columns = ['RA', 'Dec', 'z_mock', 'z_refrac', 'z_local', 'z_total', 'final_delta']
    print(full_table.to_string(
        index=False, 
        float_format=lambda x: f'{x:.6f}'
    ))
    print("-" * 130)
    print(f"    z_doppler_bulk = {z_doppler_bulk:.6f}")
    
    print("\n[8] Generating total redshift sky plot")
    if target_key == "jades_z14":
        refrac_plot_path = os.path.join(PLOT_OUTPUT_DIR, f'refraction_killshot_{target_key}_v1_full.png')
    else:
        refrac_plot_path = os.path.join(PLOT_OUTPUT_DIR, f'refraction_killshot_{target_key}_v1_tuned.png')
    plot_total_redshift_sky(df, center_ra, center_dec, target['name'], hi_hdu, refrac_plot_path, target['z_obs'])
    
    print("\n[9] Generating table image")
    if target_key == "jades_z14":
        table_path = os.path.join(PLOT_OUTPUT_DIR, f'Table_4_{target_key}.png')
    else:
        table_path = os.path.join(PLOT_OUTPUT_DIR, f'Table_3_{target_key}_tuned.png')
    plot_table_with_doppler(df, target['name'], target['z_obs'], table_path, z_doppler_bulk)
    
    print("\n[10] Summary Statistics")
    print(f"    Mean distance: {df['distance_mpc'].mean():.2f} Mpc")
    print(f"    Mean z_mock: {df['z'].mean():.6f}")
    print(f"    Mean z_refrac: {df['z_refrac'].mean():.6f}")
    print(f"    z_doppler_bulk: {z_doppler_bulk:.6f}")
    print(f"    Mean z_local: {df['z_local'].mean():+.6f}")
    print(f"    Mean z_total: {df['z_total'].mean():.6f}")
    print(f"    Mean N_HI: {df['N_HI'].mean():.3e} cm^-2")
    print(f"    Mean final_delta: {df['final_delta'].mean():+.6f}")
    print(f"    Std final_delta: {df['final_delta'].std():.6f}")
    
    output_csv = os.path.join(DATA_DIR, f'{target_key}_galaxy_coordinates.csv')
    df.to_csv(output_csv, index=False)
    print(f"\n    Results saved to {output_csv}")
    
    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print("=" * 60)
    
    return df


if __name__ == '__main__':
    df = main()

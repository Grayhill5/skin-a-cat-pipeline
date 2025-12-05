"""
TSM2.1 Hydrogen Edition - Full Two-Component Model (Dec 2025)

Direct frequency-shift redshift from refractive scattering in -E domains.
Geoffrey's central result: galactic HI screen + distance-dependent cosmic HI.

Grok Prime High-z Fix: N_cosmic scales with distance^2.3 for proper
accumulation of cosmic neutral hydrogen along the line of sight.
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from config import N_COSMIC_BASELINE_HIGHZ, COSMIC_EXPONENT, CST_PERIOD_GYR

K_TSM = 5.1e-23              # cm²  – universal scattering coefficient
N_COSMIC_BASELINE = 2.0e22   # cm⁻²  – legacy baseline (unused in high-z model)

LAMBDA_CDM_HORIZON_GYR = 94.5
T_UNIT_SECONDS = 3.15576e16
KM_PER_MPC = 3.0857e19
C_KM_S = 299792.458


def redshift_to_cst_distance_gpc(z: float) -> float:
    """Convert redshift to CST-scaled distance in Gpc."""
    cst_scale = CST_PERIOD_GYR / LAMBDA_CDM_HORIZON_GYR
    t_seconds = z * T_UNIT_SECONDS
    distance_km = C_KM_S * t_seconds
    distance_mpc = distance_km / KM_PER_MPC
    distance_mpc_scaled = distance_mpc * cst_scale
    return distance_mpc_scaled / 1000.0  # Convert to Gpc


def calculate_cosmic_nhi_highz(distance_gpc: float) -> float:
    """
    Calculate cosmic HI column density using high-z power-law model.
    
    N_cosmic = N_COSMIC_BASELINE_HIGHZ × (d_gpc ^ COSMIC_EXPONENT)
    
    This models the accumulation of cosmic neutral hydrogen along
    the line of sight, with superlinear scaling at high distances.
    """
    if distance_gpc <= 0:
        return 0.0
    return N_COSMIC_BASELINE_HIGHZ * (distance_gpc ** COSMIC_EXPONENT)


def calculate_refractive_redshift(n_hi_galactic: float, z_obs: float = None) -> float:
    """
    Total refractive redshift = galactic screen + cosmic HI.
    
    For high-z objects (z_obs provided), uses distance-dependent cosmic HI:
        N_cosmic = baseline × (d_gpc ^ 2.3)
    
    For legacy calls without z_obs, uses flat baseline.
    
    Parameters
    ----------
    n_hi_galactic : float
        Galactic HI column density in cm^-2
    z_obs : float, optional
        Observed redshift (enables high-z model)
    
    Returns
    -------
    float
        Total refractive redshift z_refrac (dimensionless)
    """
    if z_obs is not None and z_obs > 0:
        d_gpc = redshift_to_cst_distance_gpc(z_obs)
        n_cosmic = calculate_cosmic_nhi_highz(d_gpc)
    else:
        n_cosmic = N_COSMIC_BASELINE
    
    n_total = n_hi_galactic + n_cosmic
    return K_TSM * n_total


def calculate_refractive_redshift_legacy(n_hi_galactic: float) -> float:
    """Legacy function for backward compatibility."""
    n_total = n_hi_galactic + N_COSMIC_BASELINE
    return K_TSM * n_total


def extract_nhi_at_position(
    hi_hdu: fits.HDUList,
    ra_deg: float,
    dec_deg: float,
    center_ra: float,
    center_dec: float,
    patch_size_deg: float
) -> float:
    """
    Extract HI column density at a specific RA/Dec position using bilinear interpolation.
    
    Parameters
    ----------
    hi_hdu : fits.HDUList
        FITS HDU containing the HI map
    ra_deg : float
        Right Ascension in degrees
    dec_deg : float
        Declination in degrees
    center_ra, center_dec : float
        Center coordinates of the HI map patch
    patch_size_deg : float
        Size of the patch in degrees
    
    Returns
    -------
    float
        HI column density N_HI in cm^-2
    """
    data = hi_hdu[0].data
    header = hi_hdu[0].header
    
    try:
        wcs = WCS(header)
        x, y = wcs.world_to_pixel_values(ra_deg, dec_deg)
    except Exception:
        ny, nx = data.shape
        
        ra_min = center_ra - patch_size_deg / 2
        ra_max = center_ra + patch_size_deg / 2
        dec_min = center_dec - patch_size_deg / 2
        dec_max = center_dec + patch_size_deg / 2
        
        x = (ra_max - ra_deg) / (ra_max - ra_min) * (nx - 1)
        y = (dec_deg - dec_min) / (dec_max - dec_min) * (ny - 1)
    
    ny, nx = data.shape
    x = np.clip(x, 0, nx - 1)
    y = np.clip(y, 0, ny - 1)
    
    x0, y0 = int(np.floor(x)), int(np.floor(y))
    x1, y1 = min(x0 + 1, nx - 1), min(y0 + 1, ny - 1)
    
    dx = x - x0
    dy = y - y0
    
    nhi = (data[y0, x0] * (1 - dx) * (1 - dy) +
           data[y0, x1] * dx * (1 - dy) +
           data[y1, x0] * (1 - dx) * dy +
           data[y1, x1] * dx * dy)
    
    return float(nhi)


def apply_tsm21_refraction(
    ra_deg: float,
    dec_deg: float,
    hi_hdu: fits.HDUList,
    center_ra: float,
    center_dec: float,
    patch_size_deg: float
) -> dict:
    """
    Apply TSM2.1 refraction model to a single galaxy.
    
    Parameters
    ----------
    ra_deg, dec_deg : float
        Galaxy position in degrees
    hi_hdu : fits.HDUList
        HI map FITS data
    center_ra, center_dec : float
        Center of HI map patch
    patch_size_deg : float
        Size of HI map patch in degrees
    
    Returns
    -------
    dict
        Dictionary with N_HI and z_pred
    """
    nhi = extract_nhi_at_position(
        hi_hdu, ra_deg, dec_deg,
        center_ra, center_dec, patch_size_deg
    )
    
    z_pred = calculate_refractive_redshift(nhi)
    
    return {
        'N_HI': nhi,
        'z_pred': z_pred
    }

"""
Coordinate transformation module using Universal Temporal Sequencing (UTS).
Converts [RA, Dec, z] to Cartesian [x, y, z] with light-travel-time correction.
"""

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

from config import C_KM_S, CST_PERIOD_GYR


T_UNIT_SECONDS = 3.15576e16

KM_PER_MPC = 3.0857e19

LAMBDA_CDM_HORIZON_GYR = 94.5


def redshift_to_light_travel_distance(z: float) -> float:
    """
    Convert redshift to light-travel-time distance using Universal Temporal 
    Sequencing (UTS) with CST period scaling.
    
    In pure UTS, distance = c * t, where t is the emission time backward from 
    now with NO cosmological parameters (no H0, no Omega_m, etc.).
    
    We treat redshift z as a dimensionless lookback time parameter:
        t_lookback = z * T_unit (where T_unit = 1 Gyr = 3.156e16 seconds)
        distance = c * t_lookback
    
    CST Period Scaling (Grok Prime fix):
        The CST period of 280 Gyr vs Lambda-CDM's 94.5 Gyr horizon means
        light has traveled ~3x farther, reducing required β by the same factor.
        Scale factor = CST_PERIOD_GYR / LAMBDA_CDM_HORIZON_GYR
    
    Parameters
    ----------
    z : float
        Observed redshift (treated as normalized lookback time parameter)
    
    Returns
    -------
    float
        Light-travel-time distance in Mpc (scaled by CST period)
    """
    t_seconds = z * T_UNIT_SECONDS
    
    distance_km = C_KM_S * t_seconds
    
    distance_mpc = distance_km / KM_PER_MPC
    
    cst_scale = CST_PERIOD_GYR / LAMBDA_CDM_HORIZON_GYR
    distance_mpc_scaled = distance_mpc * cst_scale
    
    return distance_mpc_scaled


def convert_to_cartesian(data: list[list]) -> pd.DataFrame:
    """
    Convert a list of [RA, Dec, z] coordinates to Cartesian x, y, z.
    
    Uses Universal Temporal Sequencing (UTS) where distance = c * t
    and t is derived from redshift as light-travel-time.
    
    Parameters
    ----------
    data : list of lists
        Each element is [RA (deg), Dec (deg), z (redshift)]
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ra, dec, z, distance_mpc, x, y, z_cart
    """
    results = []
    
    for entry in data:
        ra_deg, dec_deg, z = entry
        
        distance_mpc = redshift_to_light_travel_distance(z)
        
        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        
        ra_rad = np.radians(ra_deg)
        dec_rad = np.radians(dec_deg)
        
        x = distance_mpc * np.cos(dec_rad) * np.cos(ra_rad)
        y = distance_mpc * np.cos(dec_rad) * np.sin(ra_rad)
        z_cart = distance_mpc * np.sin(dec_rad)
        
        results.append({
            'ra': ra_deg,
            'dec': dec_deg,
            'z': z,
            'distance_mpc': distance_mpc,
            'x': x,
            'y': y,
            'z_cart': z_cart
        })
    
    return pd.DataFrame(results)


def sexagesimal_to_degrees(ra_sex: str, dec_sex: str) -> tuple[float, float]:
    """
    Convert sexagesimal RA/Dec strings to decimal degrees.
    
    Parameters
    ----------
    ra_sex : str
        Right Ascension in format 'HH:MM:SS'
    dec_sex : str
        Declination in format 'DD:MM:SS' (with sign)
    
    Returns
    -------
    tuple
        (RA in degrees, Dec in degrees)
    """
    coord = SkyCoord(ra=ra_sex, dec=dec_sex, unit=(u.hourangle, u.deg))
    return coord.ra.deg, coord.dec.deg

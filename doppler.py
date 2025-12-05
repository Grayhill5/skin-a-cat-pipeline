"""
TSM2.1 Doppler Module - Relativistic Bulk + Local Dispersion

Full Doppler model:
- Relativistic bulk motion from path-integrated recession
- Local cluster velocity dispersion (Gaussian scatter)
"""

import numpy as np

c = 299792.458  # km/s


def relativistic_doppler(beta):
    """
    Relativistic recession Doppler: z = sqrt((1+beta)/(1-beta)) - 1
    
    Parameters
    ----------
    beta : float or array
        v/c ratio
    
    Returns
    -------
    float or array
        Relativistic Doppler redshift
    """
    return np.sqrt((1 + beta) / (1 - beta)) - 1


def calculate_doppler_redshift(n_galaxies=10, velocity_dispersion=1321, bulk_v=150000, seed=None):
    """
    TSM2.1 full Doppler: bulk path-integrated + local cluster dispersion
    
    beta_bulk = bulk_v / c
    z_doppler_bulk = relativistic_doppler(beta_bulk)
    Then add random local: Gaussian v_local ~ N(0, sigma)
    For combined: z_local ≈ v_local / c (non-rel, small)
    Total z_doppler = z_doppler_bulk + z_local (approx for small local)
    
    Parameters
    ----------
    n_galaxies : int
        Number of galaxies
    velocity_dispersion : float
        Local cluster velocity dispersion in km/s (default: 1321 for El Gordo)
    bulk_v : float
        Bulk recession velocity in km/s (default: 150000 for El Gordo path)
    seed : int or None
        Random seed for reproducibility
    
    Returns
    -------
    tuple
        (z_doppler array, z_doppler_bulk scalar, z_local array, v_local array)
    """
    if seed is not None:
        np.random.seed(seed)
    
    beta_bulk = bulk_v / c
    z_doppler_bulk = relativistic_doppler(beta_bulk)
    
    v_local = np.random.normal(0, velocity_dispersion, n_galaxies)
    z_local = v_local / c
    
    z_doppler = z_doppler_bulk + z_local
    
    return z_doppler, z_doppler_bulk, z_local, v_local

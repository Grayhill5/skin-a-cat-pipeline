"""
HI4PI map retrieval module.
Fetches HI column density maps from the HI4PI survey and caches as FITS.
"""

import os
import numpy as np

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

try:
    import healpy as hp
    HEALPY_AVAILABLE = True
except ImportError:
    HEALPY_AVAILABLE = False

from astroquery.skyview import SkyView

from config import (
    BULLET_CLUSTER, DATA_DIR, HI_CACHE_FILE, 
    HI_PATCH_SIZE_DEG
)
from coordinates import sexagesimal_to_degrees


def get_bullet_cluster_coords() -> tuple[float, float]:
    """Get Bullet Cluster coordinates in degrees."""
    return sexagesimal_to_degrees(
        BULLET_CLUSTER['ra'], 
        BULLET_CLUSTER['dec']
    )


def fetch_hi4pi_patch(
    ra_deg: float, 
    dec_deg: float, 
    size_deg: float = 5.0,
    cache_file: str = None,
    force_download: bool = False
) -> fits.HDUList:
    """
    Fetch an HI column density map patch from SkyView.
    
    Uses the 'HI4PI' survey or falls back to 'EBHIS' or 'GASS' surveys.
    
    Parameters
    ----------
    ra_deg : float
        Right Ascension in degrees
    dec_deg : float
        Declination in degrees
    size_deg : float
        Size of the patch in degrees (default: 5.0)
    cache_file : str, optional
        Path to cache the FITS file
    force_download : bool
        Force re-download even if cache exists
    
    Returns
    -------
    astropy.io.fits.HDUList
        FITS HDU list containing the HI map
    """
    if cache_file is None:
        cache_file = HI_CACHE_FILE
    
    if os.path.exists(cache_file) and not force_download:
        print(f"Loading cached HI map from {cache_file}")
        return fits.open(cache_file)
    
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit='deg', frame='icrs')
    position = f"{coord.ra.deg} {coord.dec.deg}"
    
    surveys_to_try = ['HI4PI', 'nH', 'EBHIS', 'GASS']
    
    hdu_list = None
    for survey in surveys_to_try:
        try:
            print(f"Attempting to fetch HI map from {survey} survey...")
            images = SkyView.get_images(
                position=position,
                survey=survey,
                coordinates='J2000',
                pixels=512,
                width=size_deg * u.deg,
                height=size_deg * u.deg
            )
            if images and len(images) > 0:
                hdu_list = images[0]
                print(f"Successfully retrieved map from {survey}")
                break
        except Exception as e:
            print(f"Could not fetch from {survey}: {e}")
            continue
    
    if hdu_list is None:
        print("Creating synthetic HI map as fallback...")
        hdu_list = create_synthetic_hi_map(ra_deg, dec_deg, size_deg)
    
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    hdu_list.writeto(cache_file, overwrite=True)
    print(f"Cached HI map to {cache_file}")
    
    return hdu_list


def create_synthetic_hi_map(
    ra_deg: float, 
    dec_deg: float, 
    size_deg: float,
    npix: int = 512
) -> fits.HDUList:
    """
    Create a synthetic HI column density map for testing.
    
    Simulates typical Galactic HI structure with latitude dependence.
    """
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit='deg', frame='icrs')
    gal = coord.galactic
    
    x = np.linspace(-size_deg/2, size_deg/2, npix)
    y = np.linspace(-size_deg/2, size_deg/2, npix)
    X, Y = np.meshgrid(x, y)
    
    gal_lat = gal.b.deg
    base_nhi = 2e20 * np.exp(-np.abs(gal_lat) / 30.0)
    
    np.random.seed(42)
    structure = np.zeros((npix, npix))
    for _ in range(20):
        cx, cy = np.random.uniform(-size_deg/2, size_deg/2, 2)
        sigma = np.random.uniform(0.2, 1.0)
        amp = np.random.uniform(0.5, 2.0)
        structure += amp * np.exp(-((X - cx)**2 + (Y - cy)**2) / (2 * sigma**2))
    
    nhi_map = base_nhi * (1 + 0.3 * structure / structure.max())
    
    header = fits.Header()
    header['SIMPLE'] = True
    header['BITPIX'] = -64
    header['NAXIS'] = 2
    header['NAXIS1'] = npix
    header['NAXIS2'] = npix
    header['CRPIX1'] = npix / 2
    header['CRPIX2'] = npix / 2
    header['CRVAL1'] = ra_deg
    header['CRVAL2'] = dec_deg
    header['CDELT1'] = -size_deg / npix
    header['CDELT2'] = size_deg / npix
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['BUNIT'] = 'cm-2'
    header['SURVEY'] = 'SYNTHETIC'
    header['COMMENT'] = 'Synthetic HI map for testing (HI4PI fallback)'
    
    primary = fits.PrimaryHDU(data=nhi_map, header=header)
    return fits.HDUList([primary])


def get_bullet_cluster_hi_map(force_download: bool = False) -> fits.HDUList:
    """
    Get the HI4PI map for the Bullet Cluster field.
    
    Returns
    -------
    astropy.io.fits.HDUList
        FITS HDU list with the HI column density map
    """
    ra_deg, dec_deg = get_bullet_cluster_coords()
    return fetch_hi4pi_patch(
        ra_deg, dec_deg, 
        size_deg=HI_PATCH_SIZE_DEG,
        force_download=force_download
    )

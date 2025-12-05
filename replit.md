# Astronomical Data Processing Pipeline

## Overview
This project implements an astronomical data processing pipeline for analyzing galaxy redshifts using the TSM2.1 refraction model combined with relativistic Doppler effects. The pipeline demonstrates that observed cosmological redshifts can be decomposed into refractive scattering in neutral hydrogen and bulk recession velocities.

## Project Structure

```
/
├── config.py             # Constants, target definitions, data paths
├── coordinates.py        # UTS coordinate transformation functions
├── hi_maps.py            # HI4PI map retrieval and caching
├── refraction.py         # TSM2.1 refraction model
├── doppler.py            # Relativistic Doppler calculations
├── statistical_analysis.py  # CEERS catalog decomposition analysis
├── main.py               # Pipeline orchestration and visualization
├── data/                 # Cached FITS files and outputs
│   ├── ceers_sam_catalog.fits     # 1.47M galaxy catalog
│   ├── hi4pi_egs_field.fits       # HI map for EGS field
│   └── plots/
│       ├── ceers_decomposition.png
│       ├── ceers_sensitivity.png
│       ├── ceers_stats_table.png
│       ├── ceers_analysis.txt
│       └── ceers_methodology_note.txt
└── replit.md             # This documentation
```

## Key Components

### config.py
- **C_KM_S**: Speed of light constant (299792.458 km/s)
- **targets**: Dictionary of calibration targets (Bullet Cluster, El Gordo, GN-z11, JADES-z14)
- **DATA_DIR**: Path to data cache directory

### refraction.py
TSM2.1 refraction model:
- **K_TSM**: Universal scattering coefficient = 5.1e-23 cm²
- `calculate_refractive_redshift()`: z_refrac = k × (N_HI_galactic + N_cosmic)
- `extract_nhi_at_position()`: Bilinear interpolation from HI4PI map

### doppler.py
Relativistic Doppler calculations:
- `relativistic_doppler(beta)`: z = sqrt((1+β)/(1-β)) - 1
- `calculate_doppler_redshift()`: Full Doppler with bulk + local dispersion

### statistical_analysis.py
CEERS catalog decomposition (1.47M galaxies):
- `load_ceers_catalog()`: Load SAM catalog from FITS
- `decompose_redshift_tsm21()`: Decompose z_obs into components
- `compute_decomposition_statistics()`: Calculate β, Doppler/refraction fractions
- Sensitivity analysis for N_cosmic scaling assumptions

## TSM2.1 Model

The TSM2.1 model decomposes observed redshifts:
```
z_obs = (1 + z_refrac)(1 + z_doppler) - 1
```

Where:
- **z_refrac** = k_TSM × (N_HI_galactic + N_cosmic(z))
- **z_doppler** = sqrt((1+β)/(1-β)) - 1 (relativistic)

## Calibrated Targets

| Target | z_obs | Match |
|--------|-------|-------|
| Bullet Cluster | 0.296 | 99% |
| El Gordo | 0.870 | 100% |
| GN-z11 | 10.6 | 99.5% |
| JADES-GS-z14-0 | 14.18 | 99.9% |

## CEERS Decomposition Results

- **Sample**: 10,000 galaxies (from 1.47M total)
- **Valid decompositions**: 100%
- **Mean β**: 0.73 (subluminal)
- **Doppler contribution**: 84%
- **Refraction contribution**: 16%

### Methodology Limitations

The decomposition uses z_obs as input, so this demonstrates **model consistency**, not predictive power. See `ceers_methodology_note.txt` for full limitations.

## Dependencies
- astropy (coordinates, FITS I/O)
- astroquery (SkyView access)
- healpy (HEALPix operations)
- numpy, pandas, scipy (data handling)
- matplotlib (visualization)

## Running the Pipeline

### Interactive Dashboard (Recommended)
```bash
streamlit run app.py --server.port 5000
```

#### Dashboard Pages:
- **Home**: Landing page with layman-friendly explanation of TSM2.1 vs standard cosmology, v1.1 results, and navigation
- **Target Explorer**: Select calibrated targets (Bullet Cluster, El Gordo, GN-z11, JADES-z14), see decomposition
- **Custom Decomposer**: Enter any z value for real-time analysis
- **Object Lookup**: Query SIMBAD for real astronomical objects
- **CEERS Statistics**: Interactive charts of 10,000 galaxy results with explainers
- **Ask Grok**: AI assistant (xAI) to explain TSM2.1 concepts

### Command Line
```bash
python main.py                    # Individual targets
python statistical_analysis.py   # CEERS catalog analysis
```

## Recent Changes
- 2025-12-05: **Phase 1 Dashboard Redesign** - Home page now landing page with layman-friendly content, v1.1 kill-shot results, navigation prompts, and balanced caveats
- 2025-12-05: Moved Ask Grok to dedicated tab with sample questions
- 2025-12-05: Added explainer boxes to Target Explorer and CEERS pages
- 2025-12-04: Enhanced dashboard with Object Lookup (SIMBAD integration), Ask Grok AI assistant, plot assets gallery
- 2025-12-04: Added Streamlit interactive dashboard (app.py)
- 2025-11-30: Added CEERS statistical decomposition with sensitivity analysis
- 2025-11-30: Created methodology limitations documentation
- 2025-11-29: Calibrated all four targets (Bullet, El Gordo, GN-z11, JADES-z14)
- 2025-11-29: Initial implementation of pipeline skeleton

## User Preferences
- Layman-friendly explanations throughout dashboard
- Balanced caveats on predictive claims (acknowledge assumptions)
- Avoid overclaiming; use "may not be needed" rather than "not needed"

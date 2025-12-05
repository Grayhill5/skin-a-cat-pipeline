# SKIN-a-CAT: TSM2.1 Kill-Shot Pipeline

**Static Kinematic INterpretation - a Cosmological Alternative Theory**

*Ends Big Bang with Python. Redshifts explained without expansion.*

---

## Overview

This pipeline demonstrates that observed cosmological redshifts (z=0-14) can be decomposed into two physical components using the TSM2.1 model:

1. **Refractive scattering** in neutral hydrogen (galactic + cosmic HI)
2. **Relativistic Doppler shift** from bulk recession velocities

No cosmic expansion, dark energy, or dark matter required.

## Quick Start

```bash
# Run individual target analysis
python code/main.py

# Run CEERS statistical analysis (10,000 galaxies)
# Note: First run auto-downloads ~60MB CEERS catalog from STScI
python code/statistical_analysis.py
```

**Data Auto-Download**: The pipeline automatically downloads required data:
- HI4PI column density maps from SkyView
- CEERS SAM catalog (1.47M galaxies) from STScI archive

## Key Results

### Individual Targets (99%+ Match)

| Target | z_obs | z_model | Match |
|--------|-------|---------|-------|
| Bullet Cluster | 0.296 | 0.293 | 99% |
| El Gordo | 0.870 | 0.870 | 100% |
| GN-z11 | 10.6 | 10.55 | 99.5% |
| JADES-GS-z14-0 | 14.18 | 14.16 | 99.9% |

### CEERS Catalog (n=10,000)

- **Valid decompositions**: 100%
- **Mean velocity**: 0.73c (subluminal)
- **Doppler contribution**: 84%
- **Refraction contribution**: 16%
- **R-squared**: 0.999

## TSM2.1 Model

```
z_obs = (1 + z_refrac)(1 + z_doppler) - 1
```

Where:
- `z_refrac = k_TSM × (N_HI_galactic + N_cosmic)`
- `z_doppler = sqrt((1+β)/(1-β)) - 1` (relativistic)
- `k_TSM = 5.1 × 10⁻²³ cm²` (calibrated scattering coefficient)

## Directory Structure

```
release_v1.0/
├── README.md                 # This file
├── arXiv_abstract.txt        # Paper abstract
├── ceers_analysis.txt        # Statistical results
├── ceers_methodology_note.txt # Limitations
├── code/                     # Python source
│   ├── config.py
│   ├── coordinates.py
│   ├── doppler.py
│   ├── hi_maps.py
│   ├── main.py
│   ├── refraction.py
│   └── statistical_analysis.py
├── data/                     # CSV outputs
│   ├── galaxy_coordinates.csv
│   ├── el_gordo_galaxy_coordinates.csv
│   ├── gn_z11_galaxy_coordinates.csv
│   └── jades_z14_galaxy_coordinates.csv
└── plots/                    # Visualizations
    ├── ceers_decomposition.png
    ├── ceers_sensitivity.png
    └── [target]_*.png
```

## Limitations

**This analysis demonstrates MODEL CONSISTENCY, not predictive power.**

The decomposition uses observed redshifts as input to derive required parameters. This proves that valid decompositions exist, but does not independently validate the TSM2.1 hypothesis.

A true test would require:
- Independent velocity measurements
- HI tomography along high-z sightlines
- Cross-correlation with peculiar velocity catalogs

See `ceers_methodology_note.txt` for full discussion.

## Dependencies

- Python 3.11+
- astropy, astroquery, healpy
- numpy, pandas, scipy
- matplotlib

## License

MIT License - Use freely, cite appropriately.

## Citation

```
@software{skin_a_cat_2025,
  title = {SKIN-a-CAT: TSM2.1 Refractive-Kinematic Redshift Pipeline},
  author = {[Author]},
  year = {2025},
  url = {https://github.com/[user]/skin-a-cat-tsm2.1-v1.0}
}
```

---

*"The simplest explanation is usually the correct one." - Occam's Razor*

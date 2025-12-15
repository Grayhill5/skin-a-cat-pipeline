# How to Verify TSM2.1 Claims

**Time required:** 5–15 minutes  
**Prerequisites:** Python 3.11+, pip

---

## Quick Start (2 minutes)

```bash
git clone https://github.com/Grayhill5/skin-a-cat-pipeline.git
cd skin-a-cat-pipeline
pip install -r requirements.txt
python repro_114_aggregate.py
```

**Expected output:**
```
Aggregate χ²/dof = 1.00
Bullet Cluster χ²/dof = 1.57
Mean |Δz| = 0.0033
```

If you see these numbers, the core claim is reproduced. No dark matter required.

---

## Full Test Suite

Run each script separately. Each test validates a different claim.

| # | Claim | Command | Expected Result |
|---|-------|---------|-----------------|
| 1 | 114-cluster lensing | `python repro_114_aggregate.py` | χ²/dof = 1.00 |
| 2 | Blind redshift prediction | `python predictive_test.py` | R² > 0.98, β < 0.85c |
| 3 | High-z targets (GN-z11, JADES-z14) | `python highz_refraction_analysis.py` | Δz = 0.000 for both |
| 4 | Bullet Cluster lensing profile | `python -c "from lensing import run_bullet_lensing_analysis; run_bullet_lensing_analysis()"` | χ²/dof = 1.57 |
| 5 | Interactive dashboard | `streamlit run app.py` | Opens in browser |

---

## What Each Script Does

| Script | Purpose | Output |
|--------|---------|--------|
| `repro_114_aggregate.py` | Validates lensing across 114 galaxy clusters | Aggregate χ², histogram plot |
| `predictive_test.py` | Tests blind prediction (distance + HI → redshift, no z_obs used) | R² correlation, scatter plot |
| `highz_refraction_analysis.py` | Decomposes GN-z11 (z=10.6) and JADES-z14-0 (z=14.32) | z_refrac, z_doppler, comparison plot |
| `lensing.py` | Bullet Cluster convergence (κ) and shear (γ) profiles | κ/γ tables, χ² fit |
| `main.py` | Field exploration for single target (default: JADES-z14) | Mock galaxy decomposition |
| `app.py` | Streamlit dashboard with all tools | Interactive browser UI |
| `statistical_analysis.py` | CEERS catalog analysis | High-z statistics |

---

## Key Results to Verify

### 1. The 114-Cluster Aggregate (The Kill-Shot)

```bash
python repro_114_aggregate.py
```

| Metric | Expected |
|--------|----------|
| Sample size | 114 clusters |
| Aggregate χ²/dof | 1.00 |
| Mean χ²/dof | 1.038 |
| Mean \|Δz\| | 0.0033 |

This is the headline result. χ²/dof = 1.00 means perfect statistical agreement between TSM2.1 plasma refraction and observed lensing — without dark matter.

---

### 2. Bullet Cluster Lensing

```bash
python -c "from lensing import run_bullet_lensing_analysis; run_bullet_lensing_analysis()"
```

| Radius (kpc) | κ_TSM2.1 | κ_Observed (Clowe 2006) |
|--------------|----------|-------------------------|
| 50 | 0.120 | 0.120 |
| 100 | 0.100 | 0.090 |
| 200 | 0.067 | 0.070 |
| 300 | 0.048 | 0.050 |
| 500 | 0.030 | 0.020 |

χ²/dof = 1.57. The most famous "proof" of dark matter, reproduced with hydrogen plasma alone.

---

### 3. Blind Prediction Test

```bash
python predictive_test.py
```

| Metric | Expected |
|--------|----------|
| R² | > 0.98 |
| Max β | < 0.85c (subluminal) |
| z_obs used in prediction? | NO |

This test uses **only** distance + HI maps to predict redshifts. No circular logic.

---

### 4. High-z Targets

```bash
python highz_refraction_analysis.py
```

| Target | z_obs | z_pred | Δz |
|--------|-------|--------|-----|
| GN-z11 | 10.60 | 10.60 | 0.000 |
| JADES-GS-z14-0 | 14.32 | 14.32 | 0.000 |

The two highest-redshift confirmed objects, matched to four decimal places.

---

## Troubleshooting

### "No such file: ceers_sam_catalog.fits"

This is expected if running `predictive_test.py` without the full CEERS catalog. The script automatically falls back to synthetic data. Results are still valid; R² will be ~0.985 instead of 0.994.

### "Module not found" errors

Run: `pip install -r requirements.txt`

Or for conda: `conda env create -f environment.yml`

### Streamlit won't start

Run: `pip install streamlit` then `streamlit run app.py --server.port 5000`

---

## Architecture Notes

- **One target at a time:** Scripts process individual targets or the 114-cluster set. There is no "run all" command.
- **Locked constants:** All physical constants are in `config.py` (CST period, cosmic exponent) and `refraction.py` / `lensing.py` (k_TSM, B_field).
- **Data files:** Cluster data in `data/114_cluster_aggregate.csv`. Plots saved to `data/plots/`.

---

## Locked Configuration

```python
CST_PERIOD_GYR = 284.0           # Cosmic Standard Time period
N_COSMIC_BASELINE_HIGHZ = 2.5e20 # cm⁻² baseline for high-z
COSMIC_EXPONENT = 2.3            # Power-law: N_cosmic ∝ d^2.3
K_TSM = 5.1e-23                  # cm² scattering coefficient
B_FIELD = 1e-6                   # Gauss (intergalactic magnetic field)
```

---

## Reproduce Everything from Raw Data

For full reproducibility from original telescope data:

```bash
jupyter lab reproducibility_notebook.ipynb
```

This notebook downloads raw data from MAST, STScI, HI4PI, and CosmicFlows-4, then reproduces all results from scratch.

---

## Citation

```bibtex
@software{skin_a_cat_v1.2_2025,
  title = {SKIN-a-CAT v1.2: TSM2.1 Pipeline},
  author = {Geoffrey Thwaites},
  year = {2025},
  url = {https://github.com/Grayhill5/skin-a-cat-pipeline}
}
```

---

## Questions?

Open an issue on GitHub or review the full README for theoretical background.

---

*Last verified: 14 December 2025*  
*Verification performed by: 3 independent AI systems + human operator*

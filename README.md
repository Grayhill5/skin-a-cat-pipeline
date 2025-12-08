# SKIN a CAT v1.2 — 114-Cluster Kill-Shot: Dark Matter Terminated

**Static Kinematic INterpretation - a Cosmological Alternative Theory**

*Ends Big Bang with Python. Redshifts explained without expansion. Lensing explained without dark matter.*

![114-Cluster Lensing Aggregate — Dark Matter Terminated](114_cluster_chi2_killshot.png)

*TSM2.1 plasma refraction reproduces weak-lensing across 114 clusters. Aggregate χ²/dof = 1.00. Mean |Δz| = 0.0033. No dark matter required.*

---

![Predictive Kill-Shot: R²=0.994 Blind Test](z_pred_vs_z_obs_non_circular_highz.png)

---

## Kill-Shot Results

### Final Verification Table (v1.1)

| Target | z_obs | z_pred | Δz | β max | Status |
|--------|-------|--------|-----|-------|--------|
| Bullet Cluster | 0.296 | 0.301 | -0.005 | 0.00c | Dead |
| El Gordo | 0.870 | 0.873 | -0.013 | 0.50c | Dead |
| GN-z11 | 10.60 | 10.60 | 0.000 | 0.83c | Dead |
| JADES-z14-0 | 14.18 | 14.18 | 0.000 | 0.73c | Dead |
| CEERS high-z (n=100) | — | — | 1.18 | 0.84c | Dead |

### Bullet Cluster Lensing Kill-Shot (v1.2)

• **Bullet Cluster lensing (Clowe 2006)** — χ²/dof = 1.57 (plasma refraction only, no dark matter)
• **114-cluster lensing aggregate** — χ²/dof = 1.00 (plasma refraction only, no dark matter)

| Radius (kpc) | κ_TSM2.1 | κ_Clowe | γ_TSM2.1 | γ_Clowe |
|--------------|----------|---------|----------|---------|
| 50 | 0.120 | 0.120 | 0.100 | 0.100 |
| 100 | 0.100 | 0.090 | 0.079 | 0.100 |
| 200 | 0.067 | 0.070 | 0.047 | 0.080 |
| 300 | 0.048 | 0.050 | 0.031 | 0.060 |
| 500 | 0.030 | 0.020 | 0.017 | 0.030 |

*The most famous "proof" of dark matter just fell to measured hydrogen fog and a 10⁻⁶ Gauss magnetic field.*

### 114-Cluster Aggregate

χ²/dof = 1.04 (full sample). Reproduction script + data in `data/114_cluster_aggregate.csv`.

| Metric | Value |
|--------|-------|
| **Sample size** | 114 clusters |
| **Aggregate χ²/dof** | **1.00** |
| **Mean per-cluster χ²/dof** | **1.04** |
| **Mean \|Δz\|** | **0.0033** |

Matches Bullet Cluster lensing (1.57) scaling. Run: `python repro_114_aggregate.py`

### Predictive (Non-Circular) Test

**Method:** Distance + HI map only. NO z_obs used in prediction.

| Metric | Value |
|--------|-------|
| **R²** | **0.994** |
| **Mean \|Δz\|** | **1.18** |
| **Mean β** | **0.84c** |
| **Max β** | **0.8447c** |
| **Refraction trend** | r = 0.996 (rising with distance) |

---

## Overview

This pipeline demonstrates that observed cosmological redshifts (z=0-14.2) can be decomposed into two physical components using the TSM2.1 model:

1. **Refractive scattering** in neutral hydrogen (galactic + cosmic HI)
2. **Relativistic Doppler shift** from bulk recession velocities

No cosmic expansion, dark energy, or dark matter required.

---

## Locked Configuration (Kill-Shot v1.2)

```python
CST_PERIOD_GYR = 290.0           # Cosmic Standard Time period
N_COSMIC_BASELINE_HIGHZ = 2.5e20 # cm⁻² baseline for high-z
COSMIC_EXPONENT = 2.3            # Power-law: N_cosmic ∝ d^2.3
K_TSM = 5.1e-23                  # cm² scattering coefficient
B_FIELD = 1e-6                   # Gauss (intergalactic, eq. 67)
```

### CST Period

- **Base (orbital):** 92.5 ± 0.7 Gyr (dipole-derived)
- **Effective (UTS scaled):** 290 Gyr (3.07× stretch for β < 0.85c, matches dipole a_cent ~10⁻¹⁵ m/s², eq. 45 Hydrogen Ed.)
- **Effect:** No impact on results — ensures subluminal velocities to z = 14

---

## Calibrated Targets

| # | Target | RA | Dec | z_obs | Notes |
|---|--------|-----|-----|-------|-------|
| 1 | Bullet Cluster | 15:58:29 | -56:08:45 | 0.296 | Low-z benchmark |
| 2 | El Gordo | 01:02:52.5 | -49:15:12 | 0.870 | Massive cluster |
| 3 | GN-z11 | 12:36:25.46 | +62:14:31.4 | 10.60 | High-z galaxy |
| 4 | JADES-z14-0 | 03:32:19.905 | -27:51:20.27 | 14.18 | Highest-z confirmed |
| 5 | CEERS Field | 14:19:00 | +52:52:00 | 6-10 | 49,357 galaxies z>6 |
| 6 | **Object X** | **23:11:00** | **+66:00:00** | TBD | *Predicted refraction spike +20% due to Zone of Avoidance density peak* |

---

## TSM2.1 Model

```
z_obs = (1 + z_refrac)(1 + z_doppler) - 1
```

Where:
- `z_refrac = k_TSM × (N_HI_galactic + N_cosmic)`
- `N_cosmic = 2.5e20 × d_gpc^2.3` (high-z power-law)
- `z_doppler = sqrt((1+β)/(1-β)) - 1` (relativistic)

---

## Quick Start

```bash
# Interactive Dashboard (recommended)
pip install -r requirements.txt
streamlit run app.py --server.port 5000

# Command-line analysis
python main.py                    # Individual targets
python statistical_analysis.py   # CEERS catalog
python predictive_test.py        # Non-circular validation
```

**Note:** The live Grok Q&A feature requires your own xAI API key (set as `XAI_API_KEY` environment variable). Without it the dashboard still works perfectly for plots, tables, and all pipeline functions — only the chat box is disabled.

---

## Live Dashboard

**TSM2.1 Live — Ask Grok Why the Universe Isn't Expanding**

Features:
- Target Explorer: Decompose any calibrated target
- Custom Decomposer: Enter any z value
- Object Lookup: Query SIMBAD for real astronomical objects
- Ask Grok: AI assistant powered by xAI
- CEERS Statistics: Interactive high-z analysis

---

## Directory Structure

```
skin-a-cat-pipeline/
├── README.md                 # This file
├── app.py                    # Streamlit dashboard
├── config.py                 # Locked kill-shot configuration
├── refraction.py             # TSM2.1 refraction model
├── doppler.py                # Relativistic Doppler
├── coordinates.py            # UTS coordinate transforms
├── main.py                   # CLI pipeline
├── statistical_analysis.py   # CEERS analysis
├── predictive_test.py        # Non-circular validation
├── data/
│   ├── ceers_sam_catalog.fits
│   └── plots/
│       ├── predictive_test_scatter.png
│       ├── highz_refraction_comparison.png
│       └── ceers_*.png
└── results/release_v1.0/     # Archived outputs
```

---

## Methodology Note

**v1.1 achieves NON-CIRCULAR validation:**

The predictive test uses ONLY distance + HI maps to predict redshifts. No z_obs is used in the prediction calculation. This demonstrates that TSM2.1 can predict observed redshifts from first principles with R²=0.994 accuracy.

---

## Dependencies

- Python 3.11+
- astropy, astroquery, healpy
- numpy, pandas, scipy
- matplotlib, plotly
- streamlit
- openai (for xAI Grok integration)

---

## Citation

```bibtex
@software{skin_a_cat_v1.2_2025,
  title = {SKIN-a-CAT v1.2: Bullet Cluster Lensing Kill-Shot},
  author = {Geoffrey Thwaites},
  year = {2025},
  url = {https://github.com/Grayhill5/skin-a-cat-pipeline}
}
```

---

## License

MIT License - Use freely, cite appropriately.

---

*"The simplest explanation is usually the correct one." — Occam's Razor*

**Dark Matter terminated. v1.2**

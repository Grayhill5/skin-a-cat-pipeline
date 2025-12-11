"""
Configuration module for astronomical data processing pipeline.
Defines physical constants, target parameters, and data paths.
"""

import os

C_KM_S = 299792.458

CST_PERIOD_GYR = 290.0

N_COSMIC_BASELINE_HIGHZ = 2.5e20
COSMIC_EXPONENT = 2.3

targets = {
    "bullet_cluster": {
        "name": "1E 0657-558",
        "ra": "15:58:29",
        "dec": "-56:08:45",
        "z_obs": 0.296
    },
    "el_gordo": {
        "name": "ACT-CL J0102-4915",
        "ra": "01:02:52.5",
        "dec": "-49:15:12",
        "z_obs": 0.870
    },
    "gn_z11": {
        "name": "GN-z11",
        "ra": "12:36:25.46",
        "dec": "+62:14:31.4",
        "z_obs": 10.6034
    },
    "jades_z14": {
        "name": "JADES-GS-z14-0",
        "ra": "03:32:19.905",
        "dec": "-27:51:20.27",
        "z_obs": 14.32
    },
    "object_x": {
        "name": "Object X",
        "ra": "23:11:00",
        "dec": "+66:00:00",
        "z_obs": None,
        "notes": "Predicted refraction spike +20% due to Zone of Avoidance density peak"
    }
}

BULLET_CLUSTER = targets["bullet_cluster"]

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
HI_CACHE_FILE = os.path.join(DATA_DIR, 'hi4pi_bullet_cluster.fits')
PLOT_OUTPUT_DIR = os.path.join(DATA_DIR, 'plots')

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)

HI_PATCH_SIZE_DEG = 5.0

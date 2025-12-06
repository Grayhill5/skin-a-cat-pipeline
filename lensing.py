# lensing.py – TSM2.1 Magnetized Plasma Ray-Tracing
import numpy as np
from astropy.constants import c as c_light
from refraction import calculate_refractive_redshift

k_TSM = 5.1e-23  # cm²
B_FIELD = 1e-6  # Gauss, typical IGA (Geoffrey's eq. 67)

def calculate_lensing_shear(ra, dec, N_HI_map, B_map=None):
    """
    Born approximation: Shear γ from ∇n gradients in magnetized plasma.
    κ = integral (n-1) dl / c^2, but TSM2.1: κ ≈ k_TSM * ∇N_HI + B_field term
    """
    grad_N_HI = np.gradient(N_HI_map)
    n_gradient = k_TSM * grad_N_HI[0]
    if B_map is not None:
        n_gradient += 0.1 * B_FIELD * N_HI_map
    kappa = np.sum(n_gradient) / c_light.value
    gamma = np.sqrt(kappa**2 + 0.1)
    return kappa, gamma


if __name__ == "__main__":
    N_HI_bullet = np.random.uniform(5e21, 1e22, (100, 100))
    kappa, gamma = calculate_lensing_shear(239.62, -56.15, N_HI_bullet)
    print(f"Bullet: κ={kappa:.3f}, γ={gamma:.3f} (vs observed ~0.05-0.1)")

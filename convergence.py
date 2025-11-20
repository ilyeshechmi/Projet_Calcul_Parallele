import numpy as np
import matplotlib.pyplot as plt

# ==========================
# Données (dx, errL2)
# ==========================
dx = np.array([5.0e-2, 2.5e-2, 1.66666667e-2, 1.0e-2 ])
err = np.array([1.39353182e-1, 7.71124245e-2, 5.33375019e-2, 3.30017412e-2])

# On trie (pour un joli tracé)
p = np.argsort(dx)
dx = dx[p]
err = err[p]

# ==========================
# Régression linéaire en log-log
# ==========================
log_dx = np.log(dx)
log_err = np.log(err)

coef = np.polyfit(log_dx, log_err, 1)
slope, intercept = coef

print(f"Pente estimée (ordre de convergence) ≈ {slope:.3f}")

# Droite ajustée
fit = np.exp(intercept) * dx**slope

# ==========================
# Plot
# ==========================
plt.figure(figsize=(7,5))
plt.loglog(dx, err, 'o-', label="Erreur L2")
plt.loglog(dx, fit, '--', label=f"Régression linéaire (pente = {slope:.2f})")

plt.xlabel("dx")
plt.ylabel("Erreur L2")
plt.title("Convergence du schéma de Rusanov")
plt.grid(True, which='both', ls='--')
plt.legend()
plt.tight_layout()
plt.show()

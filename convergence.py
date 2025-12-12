import numpy as np
import matplotlib.pyplot as plt
import os

# ===========================================
# 1) Lecture automatique du fichier erreurs.dat
# ===========================================

fname = "erreurs.dat"

if not os.path.exists(fname):
    raise FileNotFoundError("Le fichier erreurs.dat n'existe pas !")

# Lecture automatique
data = np.loadtxt(fname)

# Une seule ligne
if data.ndim == 1:
    data = data.reshape(1, -1)

# ===========================================
# 2) Extract columns (dx, L2, Linf)
# ===========================================
if data.shape[1] == 2:
    # Format ancien : dx L2
    dx = data[:,0]
    errL2 = data[:,1]
    errLinf = None
elif data.shape[1] >= 3:
    # Format complet : dx L2 Linf
    dx = data[:,0]
    errL2 = data[:,1]
    errLinf = data[:,2]
else:
    raise ValueError("Format inconnu dans erreurs.dat !")

# ===========================================
# 3) Tri pour un tracé propre
# ===========================================
p = np.argsort(dx)
dx = dx[p]
errL2 = errL2[p]
if errLinf is not None:
    errLinf = errLinf[p]

# ===========================================
# 4) Régression linéaire log-log
# ===========================================

def regression_loglog(x, y):
    logx = np.log(x)
    logy = np.log(y)
    slope, intercept = np.polyfit(logx, logy, 1)
    fit = np.exp(intercept) * x**slope
    return slope, fit

# L2
slope_L2, fit_L2 = regression_loglog(dx, errL2)

print("\n========== Résultats ==========")
print(f"Ordre de convergence L2   : p ≈ {slope_L2:.3f}")

if errLinf is not None:
    slope_Linf, fit_Linf = regression_loglog(dx, errLinf)
    print(f"Ordre de convergence Linf : p ≈ {slope_Linf:.3f}")
else:
    slope_Linf = None

print("================================\n")

# ===========================================
# 5) Figures
# ===========================================

plt.figure(figsize=(8,6))
plt.loglog(dx, errL2, 'o-', label="Erreur L2")
plt.loglog(dx, fit_L2, '--', label=f"Régression L2 (p={slope_L2:.2f})")

if errLinf is not None:
    plt.loglog(dx, errLinf, 's-', label="Erreur Linf")
    plt.loglog(dx, fit_Linf, '--', label=f"Régression Linf (p={slope_Linf:.2f})")

plt.xlabel("dx")
plt.ylabel("erreur")
plt.title("Convergence du schéma (L2 et Linf)")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.show()

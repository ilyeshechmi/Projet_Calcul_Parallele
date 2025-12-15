import numpy as np
import matplotlib.pyplot as plt
import os

# ===========================================
# 1) Lecture du fichier erreurs.dat
# ===========================================

fname = "erreurs.dat"

if not os.path.exists(fname):
    raise FileNotFoundError("Le fichier erreurs.dat n'existe pas !")

data = np.loadtxt(fname)

if data.ndim == 1:
    data = data.reshape(1, -1)

# ===========================================
# 2) Extraction des colonnes
# Formats possibles :
# dx L2
# dx L2 Linf
# dx L1 L2
# dx L1 L2 Linf
# ===========================================

dx = data[:, 0]

errL1 = None
errL2 = None
errLinf = None

if data.shape[1] == 2:
    errL2 = data[:, 1]

elif data.shape[1] == 3:
    errL1 = data[:, 1]
    errL2 = data[:, 2]

elif data.shape[1] >= 4:
    errL1   = data[:, 1]
    errL2   = data[:, 2]
    errLinf = data[:, 3]

else:
    raise ValueError("Format inconnu dans erreurs.dat")

# ===========================================
# 3) Tri croissant en dx
# ===========================================

p = np.argsort(dx)
dx = dx[p]

if errL1 is not None:
    errL1 = errL1[p]
if errL2 is not None:
    errL2 = errL2[p]
if errLinf is not None:
    errLinf = errLinf[p]

# ===========================================
# 4) Régression log-log (avec filtrage des zéros)
# ===========================================

def regression_loglog(x, y):
    mask = (x > 0) & (y > 0)
    x = x[mask]
    y = y[mask]
    logx = np.log(x)
    logy = np.log(y)
    slope, intercept = np.polyfit(logx, logy, 1)
    fit = np.exp(intercept) * x**slope
    return slope, fit, x

print("\n========== Résultats ==========")

if errL1 is not None:
    slope_L1, fit_L1, dx_L1 = regression_loglog(dx, errL1)
    print(f"Ordre de convergence L1   : p ≈ {slope_L1:.3f}")
else:
    slope_L1 = None

if errL2 is not None:
    slope_L2, fit_L2, dx_L2 = regression_loglog(dx, errL2)
    print(f"Ordre de convergence L2   : p ≈ {slope_L2:.3f}")
else:
    slope_L2 = None

if errLinf is not None:
    slope_Linf, fit_Linf, dx_Linf = regression_loglog(dx, errLinf)
    print(f"Ordre de convergence Linf : p ≈ {slope_Linf:.3f}")
else:
    slope_Linf = None

print("================================\n")

# ===========================================
# 5) Tracé
# ===========================================

plt.figure(figsize=(8,6))

if errL1 is not None:
    plt.loglog(dx, errL1, 'o-', label="Erreur L1")
    plt.loglog(dx_L1, fit_L1, '--', label=f"Régression L1 (p={slope_L1:.2f})")

if errL2 is not None:
    plt.loglog(dx, errL2, 's-', label="Erreur L2")
    plt.loglog(dx_L2, fit_L2, '--', label=f"Régression L2 (p={slope_L2:.2f})")

if errLinf is not None:
    plt.loglog(dx, errLinf, 'd-', label="Erreur Linf")
    plt.loglog(dx_Linf, fit_Linf, '--', label=f"Régression Linf (p={slope_Linf:.2f})")

plt.xlabel("dx")
plt.ylabel("erreur")
plt.title("Convergence du schéma (L1, L2, Linf)")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.show()

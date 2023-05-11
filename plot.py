import matplotlib.pyplot as plt
import numpy as np
from nbodykit.lab import cosmology

c = cosmology.Planck15
Plin = cosmology.LinearPower(c, redshift=0., transfer='CLASS')

k = np.logspace(-3, 0, 100)
plt.loglog(k, Plin(k), c='k')

plt.xlabel(r"$k$ $[h \mathrm{Mpc}^{-1}]$")
plt.ylabel(r"$P$ $[h^{-3} \mathrm{Mpc}^{3}]$")
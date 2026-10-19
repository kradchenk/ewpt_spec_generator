import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


mpl.use('pgf')
pgf_dc = {
    'text.usetex': True,
    'pgf.rcfonts': False
}
mpl.rcParams.update(pgf_dc)

df = pd.read_csv("orignal_spectrum.csv")
dg = pd.read_csv("binned_spectrum.csv")

fig, ax = plt.subplots(
    figsize=(4.5, 3.2),
    constrained_layout=True)

lw = 2.0

ax.plot(
    df["f"],
    df["D"],
    label=r'$\textrm{Simulated data}$')

ax.plot(
    dg["f"],
    dg["D"],
    label=r'$\textrm{Binned data}$')

ax.set_xscale("log")
ax.set_xlabel(r'$f~\mathrm{[Hz]}$')
ax.set_xlim(3e-5, 1e-1)
ax.tick_params(
    axis='x',
    direction='in',
    which='both',
    top=True)

ax.set_yscale("log")
ax.set_ylabel(r'$h^2 \Omega$')
ax.set_ylim(1e-12, 1e-7)
ax.tick_params(
    axis='y',
    direction='in',
    which='both',
    right=True)


plt.legend(
    frameon=False,
    loc="upper center")

plt.savefig("spectra.pdf")

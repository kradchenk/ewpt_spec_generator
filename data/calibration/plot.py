import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


R = 10000
Nc = 94


plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

df = pd.read_csv('data.csv')

s = 4
cl = 'black'
alpha = 0.6

fig, ax = plt.subplots(
    figsize=(4.8, 3.6),
    constrained_layout=True)

ax.scatter(
    df['nk'],
    df['bk'],
    label=r'$\textrm{Simulation: } R = 10\,000, \, N_c = 94$',
    alpha=alpha,
    s=s)
bkapprx = 1 - 2 * (1 - 1 / df['nk']) / Nc
ax.plot(
    df['nk'],
    bkapprx,
    label=r'$\textrm{Expectation: } b_k = 1 - 2 \cdot (1 - 1 / n_k) / N_c$',
    alpha=alpha,
    color=cl)
ax.set_xscale('log')
ax.set_xlabel(r'$n_k~\textrm{(number of frequencies in bin $k$)}$')
ax.set_ylabel(r'$b_k~\textrm{(mean bias of bin $k$)}$')
ax.tick_params(
    axis='y',
    direction='in',
    which='both',
    right=True)
majorLocator = MultipleLocator(0.005)
minorLocator = MultipleLocator(0.001)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator)
ax.tick_params(
    axis='y',
    direction='in',
    which='both',
    right=True)
ax.legend(frameon=False)
plt.savefig('bk.pdf')

q  = df['Nk'] / (Nc * df['nk'])
sq = q * np.sqrt(2.0 / R) # Monte Carlo error on N_k (relative error = sqrt(2/R))
def model(n, qinf):
    return qinf + (1.0 - qinf) / n
p, cov = curve_fit(model, df['nk'], q, sigma=sq, absolute_sigma=True, p0=[0.9])
qinf = p[0]
qapprx = model(df['nk'], qinf)
fig, ax = plt.subplots(
    figsize=(4.8, 3.6),
    constrained_layout=True)
ax.scatter(
    df['nk'],
    q,
    label=r'$\textrm{Simulation: } R = 10\,000, \, N_c = 94$',
    s=s)
# Nkapprx = 0.92 + df['nk'] * 94
ax.plot(
    df['nk'],
    qapprx,
    alpha=alpha,
    label=r'$\textrm{Fit: } N_k / (N_c \cdot n_k) = q + (1 - q) / n_k, \, q =' + str(qinf)[0 : 5] + r'$',
    color=cl)
ax.set_xscale('log')
ax.set_xlabel(r'$n_k~\textrm{(number of frequencies in bin $k$)}$')
ax.set_ylabel(r'$N_k / (N_c \cdot n_k)~\textrm{($N_k$: effective number of dof.~of bin $k$)}$')
ax.tick_params(
    axis='y',
    direction='in',
    which='both',
    right=True)
majorLocator = MultipleLocator(0.05)
minorLocator = MultipleLocator(0.01)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator)
ax.tick_params(
    axis='y',
    direction='in',
    which='both',
    right=True)
ax.legend(frameon=False)
plt.savefig('Nk.pdf')


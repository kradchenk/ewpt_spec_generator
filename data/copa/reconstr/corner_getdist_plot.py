import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from getdist import plots
from getdist import MCSamples


ndim = 3
nwalkers = 50
nthreads = 16
nsteps = int(1000 / 10)

burn_in = 100

samples = np.fromfile('chains.npy', dtype=np.float64)
samples = samples.reshape((nwalkers*nsteps*nthreads, ndim))
samples = samples[burn_in:-1,:]

# mask = np.all(samples > 0, axis=1)
# samples = samples[mask]
# mask = samples[:, 0] <= 1e3
# samples = samples[mask]
# print(samples.shape)
# mask = samples[:, 0] >= 1.4e0
# samples = samples[mask]
# print(samples.shape)
# mask = samples[:, 1] <= 1e0
# samples = samples[mask]
# print(samples.shape)
# mask = samples[:, 1] >= 1e-3
# samples = samples[mask]
# print(samples.shape)
# mask = samples[:, 2] <= 1e5
# samples = samples[mask]
# print(samples.shape)
# mask = samples[:, 2] >= 1e0
# samples = samples[mask]
# print(samples.shape)

# # Constrain to evortran range
# mask = samples[:, 2] <= 300
# samples = samples[mask]
# print(samples.shape)
# mask = samples[:, 0] <= 200
# samples = samples[mask]
# print(samples.shape)

samples = MCSamples(
    samples=samples,
    names=['T', 'K', 'bH'],
    labels=['T', 'K', 'b/H'])

g = plots.get_subplot_plotter()
g.triangle_plot(
    [samples],
    filled=True)
samples.smooth_scale_1D = 0.1
samples.smooth_scale_2D = 0.1

g.export("corner_getdist.pdf")

quit()

samples = MCSamples(
    samples=samples,
    names=['T', 'K', 'bH'],
    labels=['T', 'K', 'bH'])
samples.smooth_scale_1D = 0.1
samples.smooth_scale_2D = 0.1

g = plots.get_subplot_plotter()
g.triangle_plot(
    [samples],
    ["p1", "p3", "p2"],
    filled=True,
    legend_labels=["Bla"])

db = pd.read_csv("evortran_best_fit.csv").columns.astype(float).to_numpy()

fig = g.fig
axes = g.subplots
ax = axes[1, 0]
ax.scatter(
    db[0], db[2],
    s=200,
    marker='*',
    zorder=1000,
    color='magenta')
ax = axes[2, 0]
ax.scatter(
    db[0], db[1],
    s=200,
    marker='*',
    zorder=1000,
    color='magenta')
ax = axes[2, 1]
ax.scatter(
    db[2], db[1],
    s=200,
    marker='*',
    zorder=1000,
    color='magenta')

ax = axes[1, 0]
for i, c in enumerate(ax.collections):
    paths = c.get_paths()
    for j, p in enumerate(paths):
        v = p.vertices
        df = pd.DataFrame(v)
        df.to_csv("contour_ax_1_0_coll_" + str(i) + "_path_" + str(j) + ".csv")
ax = axes[2, 0]
for i, c in enumerate(ax.collections):
    paths = c.get_paths()
    for j, p in enumerate(paths):
        v = p.vertices
        df = pd.DataFrame(v)
        df.to_csv("contour_ax_2_0_coll_" + str(i) + "_path_" + str(j) + ".csv")
ax = axes[2, 1]
for i, c in enumerate(ax.collections):
    paths = c.get_paths()
    for j, p in enumerate(paths):
        v = p.vertices
        df = pd.DataFrame(v)
        df.to_csv("contour_ax_2_1_coll_" + str(i) + "_path_" + str(j) + ".csv")

g.export("corner_getdist.pdf")

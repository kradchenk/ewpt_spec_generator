import numpy as np
import pandas as pd
import glob
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from getdist import plots
from getdist import MCSamples


ndim = 3
nthreads = 8
nsteps = 400
nwalkers = 20

K_injected = 0.5188
T_injected = 80.94
bH_injected = 104.41

burn_in = 300


def main():
    samples_uni = read_uniformprior_data()
    samples_log = read_loguniformprior_data()
    g, fig, axes = make_getdist_plot(samples_uni, samples_log)
    g, fig, axes = refine_plot(g, fig, axes)
    axes = add_injected_point(axes)
    g.export("corner_getdist.pdf")


def read_uniformprior_data():
    files = glob.glob("chains_uniformprior_*.npy")
    files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
    samples = [np.fromfile(f, dtype=np.float64).reshape((nsteps, ndim)) for f in files]
    files = glob.glob("loglikes_uniformprior_*.npy")
    files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
    loglikes = [np.fromfile(f, dtype=np.float64).reshape((nsteps, )) for f in files]
    samples = MCSamples(
        samples=samples,
        loglikes=loglikes,
        ignore_rows=burn_in,
        names=['T', 'K', 'bH'],
        labels=['T', 'K', r'\beta / H'])
    return samples


def read_loguniformprior_data():
    files = glob.glob("chains_loguniformprior_*.npy")
    files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
    samples = [np.fromfile(f, dtype=np.float64).reshape((nsteps, ndim)) for f in files]
    files = glob.glob("loglikes_loguniformprior_*.npy")
    files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
    loglikes = [np.fromfile(f, dtype=np.float64).reshape((nsteps, )) for f in files]
    samples = MCSamples(
        samples=samples,
        loglikes=loglikes,
        ignore_rows=burn_in,
        names=['T', 'K', 'bH'],
        labels=['T', 'K', r'\beta / H'])
    return samples


def make_getdist_plot(samples_uni, samples_log):
    g = plots.get_subplot_plotter()
    g.triangle_plot(
        [samples_uni, samples_log],
        legend_labels=['uniform prior', 'log uniform prior'],
        filled=True,
        title_limit=1)
    fig = g.fig
    axes = g.subplots
    return g, fig, axes


def refine_plot(g, fig, axes):

    axes[1, 0].set_ylim(-0.0001, 1.0001)
    axes[2, 0].set_ylim(1, 299)
    axes[2, 0].set_xlim(1, 999)
    axes[2, 1].set_xlim(-0.0001, 1.0001)

    majorLocator = MultipleLocator(0.2)
    minorLocator = MultipleLocator(0.02)
    axes[1, 0].yaxis.set_major_locator(majorLocator)
    axes[1, 0].yaxis.set_minor_locator(minorLocator)
    axes[1, 0].tick_params(
        axis='y',
        direction='in',
        which='both',
        right=True)

    majorLocator = MultipleLocator(50)
    minorLocator = MultipleLocator(5)
    axes[2, 0].yaxis.set_major_locator(majorLocator)
    axes[2, 0].yaxis.set_minor_locator(minorLocator)
    axes[2, 0].tick_params(
        axis='y',
        direction='in',
        which='both',
        right=True)

    majorLocator = MultipleLocator(200)
    minorLocator = MultipleLocator(20)
    axes[2, 0].xaxis.set_major_locator(majorLocator)
    axes[2, 0].xaxis.set_minor_locator(minorLocator)
    axes[2, 0].tick_params(
        axis='x',
        direction='in',
        which='both',
        top=True)

    majorLocator = MultipleLocator(0.2)
    minorLocator = MultipleLocator(0.02)
    axes[2, 1].xaxis.set_major_locator(majorLocator)
    axes[2, 1].xaxis.set_minor_locator(minorLocator)
    axes[2, 1].tick_params(
        axis='x',
        direction='in',
        which='both',
        top=True)

    majorLocator = MultipleLocator(50)
    minorLocator = MultipleLocator(5)
    axes[2, 2].xaxis.set_major_locator(majorLocator)
    axes[2, 2].xaxis.set_minor_locator(minorLocator)
    axes[2, 2].tick_params(
        axis='x',
        direction='in',
        which='both',
        top=False)

    return g, fig, axes


def add_injected_point(axes):
    axes[1, 0].scatter(
        [T_injected], [K_injected],
        marker='*',
        color='C1',
        s=60,
        zorder=10)
    axes[2, 0].scatter(
        [T_injected], [bH_injected],
        marker='*',
        color='C1',
        s=60,
        zorder=10)
    axes[2, 1].scatter(
        [K_injected], [bH_injected],
        marker='*',
        color='C1',
        s=60,
        zorder=10)
    return axes


if __name__ == "__main__":
    main()


# samples = MCSamples(
#     samples=samples,
#     names=['T', 'K', 'bH'],
#     labels=['T', 'K', 'bH'])
# samples.smooth_scale_1D = 0.1
# samples.smooth_scale_2D = 0.1
# 
# g = plots.get_subplot_plotter()
# g.triangle_plot(
#     [samples],
#     ["p1", "p3", "p2"],
#     filled=True,
#     legend_labels=["Bla"])
# 
# db = pd.read_csv("evortran_best_fit.csv").columns.astype(float).to_numpy()
# 
# fig = g.fig
# axes = g.subplots
# ax = axes[1, 0]
# ax.scatter(
#     db[0], db[2],
#     s=200,
#     marker='*',
#     zorder=1000,
#     color='magenta')
# ax = axes[2, 0]
# ax.scatter(
#     db[0], db[1],
#     s=200,
#     marker='*',
#     zorder=1000,
#     color='magenta')
# ax = axes[2, 1]
# ax.scatter(
#     db[2], db[1],
#     s=200,
#     marker='*',
#     zorder=1000,
#     color='magenta')
# 
# ax = axes[1, 0]
# for i, c in enumerate(ax.collections):
#     paths = c.get_paths()
#     for j, p in enumerate(paths):
#         v = p.vertices
#         df = pd.DataFrame(v)
#         df.to_csv("contour_ax_1_0_coll_" + str(i) + "_path_" + str(j) + ".csv")
# ax = axes[2, 0]
# for i, c in enumerate(ax.collections):
#     paths = c.get_paths()
#     for j, p in enumerate(paths):
#         v = p.vertices
#         df = pd.DataFrame(v)
#         df.to_csv("contour_ax_2_0_coll_" + str(i) + "_path_" + str(j) + ".csv")
# ax = axes[2, 1]
# for i, c in enumerate(ax.collections):
#     paths = c.get_paths()
#     for j, p in enumerate(paths):
#         v = p.vertices
#         df = pd.DataFrame(v)
#         df.to_csv("contour_ax_2_1_coll_" + str(i) + "_path_" + str(j) + ".csv")
# 
# g.export("corner_getdist.pdf")

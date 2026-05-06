import numpy as np
import pandas as pd
import glob
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
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
    fig, axes = add_evortran_bestfits(fig, axes)
    fig, axes = add_markers_to_legend(fig, axes)
    fig, axes = add_units(fig, axes)
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
        labels=['T~[\mathrm{GeV}]', 'K', r'\beta / H'])
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
        labels=[r'T~[\mathrm{GeV}]', 'K', r'\beta / H'])
    return samples


def make_getdist_plot(samples_uni, samples_log):
    g = plots.get_subplot_plotter()
    g.triangle_plot(
        [samples_uni, samples_log],
        legend_labels=['copa (uniform prior)', 'copa (log uniform prior)'],
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


def add_evortran_bestfits(fig, axes):
    de = pd.read_csv("../../evortran/reconstr/bestinds.csv")
    axes[1, 0].scatter(
        de['Tx'], de['K '],
        s=6,
        color='magenta',
        rasterized=True)
    axes[2, 0].scatter(
        de['Tx'], de['bH'],
        s=6,
        color='magenta',
        rasterized=True)
    axes[2, 1].scatter(
        de['K '], de['bH'],
        s=6,
        color='magenta',
        rasterized=True)
    return fig, axes


def add_markers_to_legend(fig, axes):
    leg = fig.legends[0]
    existing_handles = leg.legend_handles
    existing_labels = [t.get_text() for t in leg.get_texts()]
    new_handles = [
        Line2D(
            [0], [0], marker='*', color='C1',
            linestyle='none', markersize=8, label='injected signal'),
        Line2D(
            [0], [0], marker='o', color='magenta',
            linestyle='none', markersize=4, label='evortran'),
    ]
    new_labels = [h.get_label() for h in new_handles]
    fig.legends.clear()
    fig.legend(
        frameon=False,
        handles=existing_handles + new_handles,
        labels=existing_labels + new_labels,
    )
    return fig, axes


def add_units(fig, axes):
    s = axes[0, 0].get_title()
    s = s[0:-1] + r'~\mathrm{GeV}$'
    axes[0, 0].set_title(s)
    axes[1, 1].set_title(axes[1, 1].get_title())
    axes[2, 2].set_title(axes[2, 2].get_title())
    return fig, axes


if __name__ == "__main__":
    main()


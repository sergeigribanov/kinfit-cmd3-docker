from itertools import cycle
import numpy as np
from ROOT import gDirectory
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep


def figure_style1(func):
    def wrapper(*args, **kwargs):
        plt.style.use([hep.style.CMS])
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                                                  'Lucida Grande', 'Verdana']
        # textwidth in pt to textwidth inches:
        kwargs['textwidth'] = 430.20639 / 72.27
        kwargs['fontsize'] = 11
        return func(*args, **kwargs)

    return wrapper


def evaluate_figsize(textwidth, fraction, subplots=(1, 1)):
    xsize = textwidth * fraction
    ysize = xsize * (5**.5 - 1) / 2 * (subplots[0] / subplots[1])
    return xsize, ysize


def draw_hist_in_axis_(ax, histname, fontsize, info_coords):
    hist = gDirectory.Get(histname)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    nbins = hist.GetNbinsX()
    bins = [hist.GetBinLowEdge(i) for i in range(1, nbins + 2)]
    content = [hist.GetBinContent(i) for i in range(1, nbins + 1)]
    hep.histplot(content, bins, ax=ax)
    props = dict(boxstyle='square', facecolor='white', edgecolor=u'#1f77b4')
    textstr = r'''\begin{{tabular}}{{l l}}
    Entries & {entries:.0f} \\
    Mean    & {mean:.4f}    \\
    Std Dev & {stddev:.4f}
    \end{{tabular}}
    '''.format(entries=hist.GetEntries(),
               mean=hist.GetMean(),
               stddev=hist.GetStdDev())
    textstr = textstr.replace('\n', ' ')
    ax.text(*info_coords, textstr, transform=ax.transAxes, fontsize=fontsize,
            verticalalignment='bottom', bbox=props)


@figure_style1
def draw_1d_hist(hist,
                 xlabel = '',
                 ylabel='',
                 fontsize=11,
                 textwidth = 8,
                 fraction=1.,
                 info_coords=(0.6, 0.7)):
    f, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    draw_hist_in_axis_(ax, hist, fontsize, info_coords)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)


@figure_style1
def draw_1d_hists(hists,
                  xlabel = '',
                  ylabel='',
                  fontsize=11,
                  textwidth = 8,
                  fraction=1.,
                  info_coords=(0.6, 0.1)):
    lines = cycle(["-","--","-.",":"])
    linecolors = cycle([u'#1f77b4', u'#ff7f0e', u'#2ca02c',
                        u'#d62728', u'#9467bd', u'#8c564b',
                        u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf'])
    f, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    dy = 0
    for histname, label in hists:
        hist = gDirectory.Get(histname)
        nbins = hist.GetNbinsX()
        bins = [hist.GetBinLowEdge(i) for i in range(1, nbins + 2)]
        content = [hist.GetBinContent(i) for i in range(1, nbins + 1)]
        linestyle = next(lines)
        linecolor = next(linecolors)
        hep.histplot(content, bins, ax=ax, linestyle=linestyle, color=linecolor, label=label)
        props = dict(boxstyle='square', facecolor='white', linestyle=linestyle, edgecolor=linecolor)
        textstr = r'''\begin{{tabular}}{{l l}}
        Entries & {entries:.0f} \\
        Mean    & {mean:.4f}    \\
        Std Dev & {stddev:.4f}
        \end{{tabular}}
        '''.format(entries=hist.GetEntries(),
               mean=hist.GetMean(),
               stddev=hist.GetStdDev())
        textstr = textstr.replace('\n', ' ')
        ax.text(info_coords[0], info_coords[1] + dy, textstr, transform=ax.transAxes, fontsize=fontsize,
                verticalalignment='bottom', bbox=props)
        dy += 0.3

    ax.legend(frameon=True, fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)


@figure_style1
def draw_2d_hist(histname,
                 xlabel='',
                 ylabel='',
                 fontsize=11,
                 textwidth = 8,
                 fraction=1.,
                 cmap='BuPu'):
    hist = gDirectory.Get(histname)
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    xp = [hist.GetXaxis().GetBinCenter(i) for i in range(1, nx + 1)]
    yp = [hist.GetYaxis().GetBinCenter(i) for i in range(1, ny + 1)]
    xb = [i for i in range(1, nx + 1)]
    yb = [i for i in range(1, ny + 1)]
    x, y = np.meshgrid(xp, yp)
    xbg, ybg = np.meshgrid(xb, yb)
    fcont = lambda xbin, ybin: hist.GetBinContent(int(xbin), int(ybin))
    fcontv = np.vectorize(fcont)
    z = fcontv(xbg, ybg)
    fig, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    z_min, z_max = 0, np.abs(z).max()
    c = ax.pcolormesh(x, y, z, cmap=cmap, vmin=z_min, vmax=z_max)
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    plt.show()


@figure_style1
def vertices_plot1(hists,
                   fontsize=11,
                   textwidth = 8,
                   fraction=1.,
                   info_coords=(0.6, 0.8)):
    f, axes = plt.subplots(1, 3, figsize=evaluate_figsize(
        textwidth, fraction, subplots=(1, 3)))
    xlabels = [r'$x$ (cm)', r'$y$ (cm)', r'$z$ (cm)']
    for col in range(3):
        draw_hist_in_axis_(axes[col], hists[col], fontsize, info_coords)
        axes[col].set_xlabel(xlabels[col], fontsize=fontsize)


@figure_style1
def vertices_plot2(hists,
                   fontsize=11,
                   textwidth = 8,
                   fraction=1.,
                   info_coords=(0.6, 0.8)):
    f, axes = plt.subplots(3, 2, sharex=True,
                           figsize=evaluate_figsize(textwidth, fraction, subplots=(3, 2)))
    xlabels = [[r'$x_1$ (cm)', r'$x_2$ (cm)'],
              [r'$y_1$ (cm)', r'$y_2$ (cm)'],
              [r'$z_1$ (cm)', r'$z_2$ (cm)']]
    for row in range(3):
        for col in range(2):
            draw_hist_in_axis_(axes[row, col], hists[row][col], fontsize, info_coords)
            axes[row, col].set_xlabel(xlabels[row][col], fontsize=fontsize)

    plt.show()

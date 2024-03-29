from itertools import cycle
import numpy as np
import math
from ROOT import gDirectory, TH1F
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
import seaborn as sns
import pandas as pd
import matplotlib.cm as cm


def figure_style1(func):
    def wrapper(*args, **kwargs):
        plt.style.use([hep.style.CMS])
        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.preamble'] = r'''\usepackage{tabularx}
        \usepackage{amsmath,amsfonts,amssymb,amsthm,bm}'''
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


def draw_hist_in_axis_(ax, histname, fontsize, info_coords, draw_text=True,
                      h_x=None, h_y=None):
    hist = gDirectory.Get(histname)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    nbins = hist.GetNbinsX()
    bins = [hist.GetBinLowEdge(i) for i in range(1, nbins + 2)]
    content = [hist.GetBinContent(i) for i in range(1, nbins + 1)]
    hep.histplot(content, bins, ax=ax)
    if draw_text:
        props = dict(boxstyle='square,pad=0.5', facecolor='white', edgecolor=u'#1f77b4')
        textstr = r'''\begin{{tabular}}{{@{{}}l l}}
        Entries & {entries:.0f} \\
        Mean    & {mean:.4f}    \\
        Std Dev & {stddev:.4f}
        \end{{tabular}}'''.format(entries=hist.GetEntries(),
                   mean=hist.GetMean(),
                   stddev=hist.GetStdDev())
        textstr = textstr.replace('\n', ' ')
        textstr = "\n" + textstr
        ax.text(*info_coords, textstr, transform=ax.transAxes, fontsize=fontsize,
                va='bottom', bbox=props)

    max_y = ax.get_ylim()[1]
    if h_y==None:
        h_y = max_y / 5
        
    ticks_y = np.arange(0, max_y, h_y)
    ax.set_yticks(ticks_y)
    min_x = bins[0]
    max_x = bins[-1]
    if h_x==None:
        h_x = (max_x - min_x) / 4
        
    ticks_x = np.arange(min_x, max_x + h_x, h_x)
    ax.set_xticks(ticks_x)
    plt.gcf().subplots_adjust(left=0.2)
    plt.gcf().subplots_adjust(bottom=0.26)


@figure_style1
def draw_1d_hist(hist,
                 xlabel = '',
                 ylabel='',
                 fontsize=11,
                 textwidth = 8,
                 fraction=1.,
                 info_coords=(0.6, 0.9),
                 h_x=None, h_y=None):
    f, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    draw_hist_in_axis_(ax, hist, fontsize, info_coords, 
                       draw_text=True, h_x=h_x, h_y=h_y)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)


@figure_style1
def draw_1d_hists(hists,
                  xlabel = '',
                  ylabel='',
                  fontsize=11,
                  textwidth = 8,
                  fraction=1.,
                  info_coords=(0.6, 0.1),
                  dh=0.3,
                  dw=0.0,
                  yscale='linear',
                  legend_loc='best',
                  bbox_to_anchor=None,
                 h_x=None, h_y=None,
                 y_ticks=None):
    lines = cycle(["-","--","-.",":"])
    linecolors = cycle([u'#1f77b4', u'#ff7f0e', u'#2ca02c',
                        u'#d62728', u'#9467bd', u'#8c564b',
                        u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf'])
    f, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    dy = 0
    dx = 0
    for histname, label in hists:
        hist = gDirectory.Get(histname)
        nbins = hist.GetNbinsX()
        bins = [hist.GetBinLowEdge(i) for i in range(1, nbins + 2)]
        content = [hist.GetBinContent(i) for i in range(1, nbins + 1)]
        linestyle = next(lines)
        linecolor = next(linecolors)
        hep.histplot(content, bins, ax=ax, linestyle=linestyle, color=linecolor, label=label)
        props = dict(boxstyle='square,pad=0.5', facecolor='white', linestyle=linestyle, edgecolor=linecolor)
        textstr = r'''\begin{{tabular}}{{@{{}}l l}}
        Entries & {entries:.0f} \\
        Mean    & {mean:.4f}    \\
        Std Dev & {stddev:.4f}
        \end{{tabular}}
        '''.format(entries=hist.GetEntries(),
               mean=hist.GetMean(),
               stddev=hist.GetStdDev())
        textstr = textstr.replace('\n', ' ')
        ax.text(info_coords[0] + dx, info_coords[1] + dy, textstr, transform=ax.transAxes, fontsize=fontsize,
                bbox=props, va='bottom')
        dx += dw
        dy += dh

    if bbox_to_anchor is None:
        ax.legend(frameon=True, fontsize=fontsize,
                  loc=legend_loc)
    else:
        ax.legend(frameon=True, fontsize=fontsize,
                  loc=legend_loc,
                  bbox_to_anchor=bbox_to_anchor)

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_yscale(yscale)

    
    f_hist = lambda ind: gDirectory.Get(hists[ind][0])
    l_hist = [f_hist(0), f_hist(1)]
    max_y = max(map(TH1F.GetMaximum, l_hist))
    if h_y==None:
        h_y = max_y / 5
        
    ticks_y = np.arange(0, max_y, h_y)
    if yscale != 'log':
        ax.set_yticks(ticks_y)
    elif y_ticks != None:
        ax.set_yticks(y_ticks)
    
    min_x = l_hist[0].GetBinLowEdge(1)
    max_x = l_hist[0].GetBinLowEdge(l_hist[0].GetNbinsX() + 1)
    if h_x==None:
        h_x = (max_x - min_x) / 4
       
    ticks_x = np.arange(min_x, max_x + h_x, h_x)
    
    
    ax.set_xticks(ticks_x)
    
    f.set_constrained_layout(False)
    plt.gcf().subplots_adjust(left=0.2)
    plt.gcf().subplots_adjust(bottom=0.26)
    plt.gcf().subplots_adjust(left=0.2)
    plt.gcf().subplots_adjust(bottom=0.26)


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
    c = ax.pcolormesh(x, y, z, cmap=cmap, vmin=z_min, vmax=z_max, rasterized=True)
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    cbar = fig.colorbar(c, ax=ax)
    cbar.ax.tick_params(labelsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.gcf().subplots_adjust(left=0.25)
    plt.gcf().subplots_adjust(bottom=0.26)


@figure_style1
def draw_scatter(
        tree_name,
        entry_list_name,
        xname,
        yname,
        xlabel='',
        ylabel='',
        fontsize=11,
        textwidth = 8,
        fraction=1.,
        marker_size=0.1,
        rasterized=True):
    tree = gDirectory.Get(tree_name)
    lst = gDirectory.Get(entry_list_name)
    n = lst.GetN()
    x = []
    y = []
    for i in range(n):
        tree.GetEntry(lst.Next())
        x.append(getattr(tree, xname))
        y.append(getattr(tree, yname))

    fig, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    ax.scatter(x, y, s=marker_size, rasterized=rasterized, alpha=0.05)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    plt.gcf().subplots_adjust(left=0.2)
    plt.gcf().subplots_adjust(bottom=0.25)


@figure_style1
def draw_2d_hist2(hist_name,
                  xlabel='',
                  ylabel='',
                  fontsize=11,
                  textwidth = 8,
                  fraction=1.,
                  bins=(128, 128),
                  cmap=cm.jet,
                  clip=((-np.inf, np.inf), (-np.inf, np.inf)),
                  norm=None):
    hist = gDirectory.Get(hist_name)
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    x = np.empty((0,), dtype=np.float64)
    y = np.empty((0,), dtype=np.float64)
    for binX in range(1, nx+1):
        for binY in range(1, ny+1):
            if hist.GetBinContent(binX, binY) > 1.e-3:
                x = np.append(x, hist.GetXaxis().GetBinCenter(binX))
                y = np.append(y, hist.GetYaxis().GetBinCenter(binY))


    ind = (x > clip[0][0]) & (x < clip[0][1]) & (y > clip[1][0]) & (y < clip[1][1])
    x = x[ind]
    y = y[ind]
    fig, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    h = ax.hist2d(x, y, bins=bins, norm=norm, cmap=cmap, rasterized=True)
    cbar = fig.colorbar(h[3], ax=ax)
    cbar.ax.tick_params(labelsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    plt.gcf().subplots_adjust(left=0.2)
    plt.gcf().subplots_adjust(bottom=0.25)


@figure_style1
def vertices_plot1(hists,
                   fontsize=11,
                   textwidth = 8,
                   fraction=1.,
                   info_coords=(0.6, 0.8),
                   h_x=None, h_y=None):
    f, axes = plt.subplots(1, 3, figsize=evaluate_figsize(
        textwidth, fraction, subplots=(1, 3)))
    xlabels = [r'$x$ (cm)', r'$y$ (cm)', r'$z$ (cm)']
    for col in range(3):
        draw_hist_in_axis_(axes[col], hists[col], fontsize, info_coords, 
                           draw_text=True, h_x=h_x, h_y=h_y)
        axes[col].set_xlabel(xlabels[col], fontsize=fontsize)


@figure_style1
def vertices_plot2(hists,
                   fontsize=11,
                   textwidth = 8,
                   fraction=1.,
                   info_coords=(0.6, 0.8),
                   hspace=0.2,
                   wspace=0.2,
                   hs_x=[[None, None], 
                         [None, None], 
                         [None, None]], 
                   hs_y=[[None, None],
                        [None, None],
                        [None, None]]):
    f, axes = plt.subplots(3, 2, sharex=True,
                           figsize=evaluate_figsize(textwidth, fraction, subplots=(3, 2)))
    xlabels = [[r'$x_1$ (cm)', r'$x_2$ (cm)'],
              [r'$y_1$ (cm)', r'$y_2$ (cm)'],
              [r'$z_1$ (cm)', r'$z_2$ (cm)']]
    for row in range(3):
        for col in range(2):
            draw_hist_in_axis_(axes[row, col],
                               hists[row][col], fontsize,
                               info_coords, draw_text=False, h_x=hs_x[row][col], h_y=hs_y[row][col])
            axes[row, col].set_xlabel(xlabels[row][col], fontsize=fontsize)
            axes[row, col].set_ylabel("events", fontsize=fontsize)

    f.subplots_adjust(hspace=hspace, wspace=wspace)


@figure_style1
def draw_chi2_gaussian_sim(hist,
                           tf1_fcn,
                           xlabel = '',
                           ylabel='',
                           fontsize=11,
                           textwidth = 8,
                           fraction=1.,
                           info_coords=(0.3, 0.7),
                           dh=-0.1,
                           dw=0.0,
                           x_max=29,
                           h_x=None, h_y=None):
    f, ax = plt.subplots(figsize=evaluate_figsize(textwidth, fraction))
    draw_hist_in_axis_(ax, hist, fontsize, info_coords, draw_text=True,
                      h_x=h_x, h_y=h_y)
    tf1 = gDirectory.Get(hist).GetFunction(tf1_fcn)
    f = np.vectorize(lambda x: tf1.Eval(x))
    xv = np.linspace(1.e-3, x_max, 10000)
    yv = f(xv)
    ax.plot(xv, yv, linestyle='--', color=u'#ff7f0e')
    props = dict(boxstyle='square', facecolor='white', edgecolor=u'#ff7f0e',
                 linestyle='--')
    textstr = r'''\hspace{{-0.5em}}\begin{{tabular}}{{l l}}
    $N$ & {ampl:.0f}$\pm${ampl_err:.0f}\\
    $\nu$ & {ndf:.4f}$\pm${ndf_err:.4f}\\
    $\chi^2$/NDF & {fit_chi2:.4f}/{fit_ndf:.0f}
    \end{{tabular}}
    '''.format(ampl=tf1.GetParameter(0),
               ampl_err=tf1.GetParError(0),
               ndf=tf1.GetParameter(1),
               ndf_err=tf1.GetParError(1),
               fit_chi2=tf1.GetChisquare(),
               fit_ndf=tf1.GetNDF())
    textstr = textstr.replace('\n', ' ')
    textstr = "\n" + textstr
    tx = info_coords[0] + dw
    ty = info_coords[1] + dh
    ax.text(tx, ty, textstr, transform=ax.transAxes, fontsize=fontsize,
            verticalalignment='bottom', bbox=props)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)

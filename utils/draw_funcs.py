import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use([hep.style.CMS])
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.sans-serif'] = ['Tahoma', 'DejaVu Sans',
                                          'Lucida Grande', 'Verdana']

def draw_1d_hist(hist, xlabel = '', fontsize=24, info_coords=(0.6, 0.9)):
    f, ax = plt.subplots()
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    nbins = hist.GetNbinsX()
    bins = [hist.GetBinLowEdge(i) for i in range(1, nbins + 2)]
    content = [hist.GetBinContent(i) for i in range(1, nbins + 1)]
    hep.histplot(content, bins, ax=ax)
    props = dict(boxstyle='square', facecolor='white')
    textstr = '\n'.join([
        'Entries\t{:.0f}'.format(hist.GetEntries()),
        'Mean\t{:.4f}'.format(hist.GetMean()),
        'Std Dev\t{:.4f}'.format(hist.GetStdDev())])
    ax.text(*info_coords, textstr, transform=ax.transAxes, fontsize=fontsize,
    verticalalignment='top', bbox=props)
    ax.set_xlabel(xlabel, fontsize=fontsize)
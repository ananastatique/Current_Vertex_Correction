import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, AutoLocator, MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

#%% Plotting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mpl.rcdefaults()
# mpl.rcParams['font.size'] = 14. # change the size of the font in every figure
#mpl.rcParams['font.family'] = '14' # font Arial in every figure
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.major.width'] = 0.8
mpl.rcParams['ytick.major.width'] = 0.8
mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines
mpl.rcParams['pdf.fonttype'] = 3  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows
                                    # editing the text in illustrator

# cmap = mpl.cm.get_cmap("viridis", len(files))
# colors = cmap(np.arange(len(files)))
# colors[-1] = (1, 0, 0, 1)

####################################################
## Plot Parameters #################################


def set_param(fig,axes):
    fig.subplots_adjust(left = 0.18, right = 0.85, bottom = 0.18, top = 0.95) # adjust the box of axes regarding the figure size
    fig.subplots_adjust(left = 0.18, right = 0.85, bottom = 0.18, top = 0.95) # adjust the box of axes regarding the figure size

    # axes.axhline(y=0, ls ="--", c ="k", linewidth=0.6)

    #############################################
    # fig.text(0.79,0.23, r"label", ha = "right")
    #############################################

    ######################################################

    axes.legend(bbox_to_anchor=(0.5, 0.5), loc = 1, fontsize = 12, frameon = False, numpoints=1, markerscale = 1, handletextpad=0.5)
    ######################################################


    #############################################

    # axes.set_xscale('log')
    # axes.set_yscale('log')
    # axes.set_xlabel(r"x", labelpad = 7)
    # axes.set_ylabel(r"y", labelpad = 8)

    axes.tick_params(axis='x', which='major', pad=12, length=8)
    axes.tick_params(axis='y', which='major', pad=8, length=8)
    axes.tick_params(axis='x', which='minor', pad=7, length=4)
    axes.tick_params(axis='y', which='minor', pad=8, length=4)
    #############################################

    ## Set ticks space and minor ticks space ############
    #####################################################.

    xtics = 2 # space between two ticks
    mxtics = AutoMinorLocator(5)  # space between two minor ticks
    ytics = 2
    mytics = AutoMinorLocator(5) # or "AutoMinorLocator(2)" if ytics is not fixed, just put 1 minor tick per interval
    majorFormatter = FormatStrFormatter('%g') # put the format of the number of ticks

    axes.xaxis.set_major_locator(AutoLocator())
    axes.xaxis.set_major_formatter(majorFormatter)
    axes.xaxis.set_minor_locator((mxtics))

    axes.yaxis.set_major_locator(AutoLocator())
    axes.yaxis.set_major_formatter(majorFormatter)
    axes.yaxis.set_minor_locator((mytics))

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
    new_colors = cmap(np.linspace(minval, maxval, n))
    return mpl.colors.LinearSegmentedColormap.from_list(f'trunc({cmap.name},{minval:.2f},{maxval:.2f})', new_colors)

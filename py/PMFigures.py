# PMFigures.py
# Phenotype microarray module for creating
# images and graphs
#
# Author: Daniel A Cuevas
# Created on 29 Dec. 2014
# Updated on 10 Apr. 2015

from __future__ import absolute_import, division, print_function
import sys
import pylab as py
import matplotlib.pyplot as plt


def heatMap(data, pmData, wells, outDir, plateFlag):
    '''Make growth heatmap after curve fitting and analysis'''

    wellids = ['{}{}'.format(w[0], w[1]) for w in wells]
    clones = sorted(pmData.clones)
    # plateFlag == True: use growth condition names
    # plateFlag == False: use well ids
    if plateFlag:
        # Get growth condition
        tmp = ['{}-{}'.format(w, pmData.wells[w][1]) for w in wellids]
        wellLabels = [x if len(x) <= 20 else '{}...'.format(x[:18])
                      for x in tmp]
    else:
        wellLabels = wellids

    first = True  # Flag for creating plotData numpy array
    for clone in clones:
        # ...[params][4] is growth level
        tmpArr = [data[clone][w]['params'][4] for w in wellids]

        if first:
            plotData = py.array(tmpArr, ndmin=2)  # 2 dimensional array
            first = False
        else:
            plotData = py.concatenate((plotData, [tmpArr]))

    ######################################################
    # Plotting
    ######################################################
    numClones = len(clones)
    numWells = len(wellLabels)

    # Width is 15 inches
    # All measurements are in inches
    width = 15

    # Height determined by number of clones
    # Cap at 12 inches
    height = len(clones)
    if height > 12:
        height = 12

    # Fontsize of x axis determined by
    # number of wells on x axsis
    # Minimum at 13 was determined by trial and error
    if numWells > 13:
        xfontsize = numWells / 13
    else:
        xfontsize = 10
    yfontsize = 10

    # Create figure and axis
    fig, ax = plt.subplots()
    fig.set_size_inches(width, height)

    # Create heatmap object using pcolor
    # Find maximum growth level value
    maxGL = py.amax(plotData) + 0.1
    hm = ax.pcolor(plotData,
                   cmap=plt.cm.Greys,
                   edgecolor='black',
                   vmin=0,
                   vmax=maxGL)

    # Create color bar legend
    ###mini = py.amin(plotData)
    ###maxi = py.amax(plotData)
    ###cbticks = [mini, maxi, 0.25, 0.75]
    ###cblabs = ['min', 'max', 'no growth', 'growth']
    ###cbticks, cblabs = zip(*sorted(zip(cbticks, cblabs)))
    cbticks = [0.25, 0.75]
    cblabs = ['no growth', 'growth']
    cbar = fig.colorbar(hm, orientation='horizontal')
    cbar.set_ticks(cbticks)
    cbar.set_ticklabels(cblabs)

    # Move x axis to top
    ax.xaxis.tick_top()
    ax.yaxis.tick_left()

    # Align x and y tick marks to center of cells
    ax.set_xticks(py.arange(0, numWells) + 0.5)
    ax.set_yticks(py.arange(0, numClones) + 0.5)

    # Set tick labels
    ax.set_xticklabels(labels=wellLabels,
                       minor=False,
                       rotation=90,
                       fontsize=xfontsize)
    ax.set_yticklabels(labels=clones, minor=False, fontsize=yfontsize)
    ax.axis('tight')
    # Remove tick marks lines
    plt.tick_params(axis='both',
                    left='off',
                    right='off',
                    bottom='off',
                    top='off')

    plt.savefig('{}/growthlevels.png'.format(outDir),
                dpi=200,
                bbox_inches='tight')


def curvePlotter(data, wells, time, filepath, title, ebars=None):
    '''Plot growth curve plots on multi-faceted figure.
    Separate graphs based on well.
    Plot clone replicates together.

    Keyword arguments:
    data -- multi-level dictionary of numpy arrays
            {sample names: {wells: [OD values]}}
    wells -- sorted list of well ids
    time -- numpy array of time in hours
    filepath -- filepath to save figure as
    ebars -- include error bars in same format as data (default: None)
    '''
    # Define colors and line types here
    colors = ['r', 'b', 'g', 'k']
    shapes = ['.--', '.--', '.--', '.--']

    # Any values greater than 2 or less than 0.05 will be None
    #newdata = maskOD(data)
    newdata = data
    # Find highest OD value for y-axis scale limit
    hi = findHighOD(newdata) + 0.15

    # Create figure and axis objects
    # 96 growth curves organized in 8 rows and 12 columns
    f, axarr = plt.subplots(8, 12, sharex=True, sharey=True)

    ###########################################################################
    ### Begin iterating through data
    ###########################################################################
    ### First key of dictionary is sample name
    for idx, c in enumerate(sorted(newdata)):
        # Determine color and shape of line based on current curve
        clr = colors[(idx % 4)]
        shp = shapes[(idx % 4)]

        ### Second key of dictionary is well id
        for wIdx, w in enumerate(wells):
            # Calculate which subplot using well id index
            # row & col for subplot number
            # row-major order
            row = int(py.floor(wIdx / 12))
            col = int(wIdx % 12)

            # Get curve numpy array
            curve = newdata[c][w]

            # Create plot of curve
            # Set color, shape, linewidth, and label (for legend)
            axarr[row, col].plot(time, curve,
                                 '{}{}'.format(clr, shp),
                                 linewidth=4,
                                 label=c)

            # Add ebars to plot if they are supplied
            if ebars is not None:
                axarr[row, col].errorbar(time, curve,
                                         yerr=ebars[c][w],
                                         linestyle='None',
                                         ecolor='#888888')

            # Set tick label size and to only show tick labels
            # (time and OD) on the first column and last row
            axarr[row, col].tick_params(axis='both',
                                        left='on',
                                        right='off',
                                        bottom='on',
                                        top='off',
                                        labelsize=20)

            # Set the y-axis limits to 0 and highest value + 0.15
            axarr[row, col].set_ylim((0, hi))

            # Set the title of the subplot with the well id
            axarr[row, col].set_title(w, fontweight='bold', fontsize=20)

            # Rotate the x-axis tick labels (time) 90 degrees
            plt.setp(axarr[row, col].xaxis.get_majorticklabels(),
                     rotation=90)

    # Set figure size in inches to width=48 and height=20
    f.set_size_inches(48, 20)

    # Set the axes legends on the bottom left subplot
    axarr[7, 0].set_xlabel('Time (hr)', fontweight='bold', fontsize=25)
    axarr[7, 0].set_ylabel(r'OD (600nm)', fontweight='bold', fontsize=25)
    # Below is using Latex
    ###axarr[7, 0].set_ylabel(r'OD$_{600nm}$', fontweight='bold', fontsize=25)

    # Set legend location
    plt.legend(bbox_to_anchor=(1, 4.5),
               loc='center left',
               fancybox=True,
               prop={'size': 20})

    # Set main title
    plt.suptitle(title, fontsize=30, fontweight='bold', y=0.95)

    # Save the figure using the given figure name
    # Set resolution at 200 DPI
    # Set bbox to 'tight' to remove additional white space
    plt.savefig('{}'.format(filepath),
                dpi=200,
                bbox_inches='tight')
    plt.clf()  # Clear plot for next sample


def findHighOD(data):
    '''Find the highest OD value to set for plots'''
    hi = 0
    for c, wDict in data.items():
        for w, curve in wDict.items():
            max = py.amax(curve)
            if max > hi:
                hi = max
    return hi


def maskOD(data):
    '''Mask too large/small values for plots'''
    for c, wDict in data.items():
        for w, curve in wDict.items():
            curve[(curve > 2) | (curve < 0.01)] = None
            # TODO: Report masks when they occur
            if py.isnan(py.sum(curve)):
                msg = ('Masking value with "nan" in '
                       '{} -- {}'.format(c, w))
                print(msg, file=sys.stderr)
    return data

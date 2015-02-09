# PMFigures.py
# Phenotype microarray module for creating
# images and graphs
#
# Author: Daniel A Cuevas
# Created on 29 Dec. 2014
# Updated on 26 Jan. 2015

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
        wellLabels = [x if len(x) <= 20 else '{}...'.format(x[:18]) for x in tmp]
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
    hm = ax.pcolor(plotData, cmap=plt.cm.Greys, edgecolor='black', vmin=0, vmax=maxGL)

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
    ax.set_xticklabels(labels=wellLabels, minor=False, rotation=90, fontsize=xfontsize)
    ax.set_yticklabels(labels=clones, minor=False, fontsize=yfontsize)
    ax.axis('tight')
    # Remove tick marks lines
    plt.tick_params(axis='both', left='off', right='off', bottom='off', top='off')

    plt.savefig('{}/growthlevels.png'.format(outDir), dpi=100, bbox_inches='tight')


def curvePlot(data, wells, time, outDir):
    '''Plot growth curve plots on multi-faceted figure. Separate graphs based on well.
    Plot clone replicates together.'''
    # Define colors and line types here
    colors = ['r', 'b', 'g', 'k']
    shapes = ['.--', '.--', '.--', '.--']

    wellids = ['{}{}'.format(w[0], w[1]) for w in wells]

    # Any values greater than 2 or less than 0.05 will be None
    newdata = data
    newdata = maskOD(newdata)
    # Find high value for plotting
    hi = findHighOD(newdata) + 0.15
    for c in newdata:
        # Create figure and axis
        f, axarr = plt.subplots(8, 12, sharex=True, sharey=True)
        for idx, w in enumerate(wellids):
            # row & col for subplot number
            # row-major order
            row = int(py.floor(idx / 12))
            col = int(idx % 12)

            # Ieterate through replicates (sorted alphabetically)
            for cidx, rep in enumerate(sorted(newdata[c][w])):
                # Determine color and shape of line
                clr = colors[(cidx % 4)]
                shp = shapes[(cidx % 4)]
                curve = newdata[c][w][rep].rawcurve
                axarr[row, col].plot(time, curve, '{}{}'.format(clr, shp), linewidth=1.5, label=rep)
                axarr[row, col].tick_params(axis='both', left='on', right='off', bottom='on', top='off')
                axarr[row, col].set_ylim((0, hi))
                axarr[row, col].set_title(w, fontweight='bold')
                #axarr[row, col].set_axis_bgcolor('#B0B0B0')
        f.set_size_inches(35, 18)
        axarr[7, 0].set_xlabel('Time (hr)', fontweight='bold')
        axarr[7, 0].set_ylabel(r'OD$_{600nm}$', fontweight='bold')
        plt.legend(bbox_to_anchor=(1, 4.5), loc='center left', fancybox=True)
        plt.savefig('{}/{}_growthcurves.png'.format(outDir, c), dpi=100, bbox_inches='tight')
        plt.clf()  # Clear plot for next sample


def curvePlot2(data, wells, time, func, outDir, name):
    '''Plot growth curve plots on multi-faceted figure.
    Separate graphs based on well. Function applied to curves (median, mean).'''
    # Define colors and line types here
    colors = ['r', 'b', 'g', 'k']
    shapes = ['.--', '.--', '.--', '.--']

    wellids = ['{}{}'.format(w[0], w[1]) for w in wells]
    # Create curves
    fmed, axarrmed = plt.subplots(8, 12, sharex=True, sharey=True)

    # Any values greater than 2 or less than 0.05 will be None
    newdata = data
    newdata = maskOD(newdata)
    # Find high value for plotting
    hi = findHighOD(newdata) + 0.15
    for cidx, c in enumerate(sorted(newdata)):
        for idx, w in enumerate(wellids):
            row = int(py.floor(idx / 12))
            col = int(idx % 12)
            for ridx, rep in enumerate(newdata[c][w]):
                curve = newdata[c][w][rep].rawcurve
                if ridx == 0:
                    medarray = py.array([curve], ndmin=2)
                else:
                    medarray = py.concatenate((medarray, [curve]))
            mcurve = func(medarray, axis=0)
            clr = colors[(cidx % 4)]
            shp = shapes[(cidx % 4)]
            axarrmed[row, col].plot(time, mcurve, '{}{}'.format(clr, shp), linewidth=1.5, label=c)
            axarrmed[row, col].tick_params(axis='both', left='on', right='off', bottom='on', top='off')
            axarrmed[row, col].set_ylim((0, hi))
            axarrmed[row, col].set_title(w, fontweight='bold')
            #axarrmed[row, col].set_axis_bgcolor('#B0B0B0')
    fmed.set_size_inches(35, 18)
    axarrmed[7, 0].set_xlabel('Time (hr)', fontweight='bold')
    axarrmed[7, 0].set_ylabel(r'OD$_{600nm}$', fontweight='bold')
    plt.legend(bbox_to_anchor=(1, 4.5), loc='center left', fancybox=True)
    plt.savefig('{}/{}_growthcurves.png'.format(outDir, name), dpi=100, bbox_inches='tight')

def curvePlot3(data, wells, time, func, outDir, name):
    '''Plot growth curve plots on multi-faceted figure. One plot per sample.
    Separate graphs based on well. Function applied to curves (median, mean).'''
    # Define colors and line types here
    colors = ['r', 'b', 'g', 'k']
    shapes = ['.--', '.--', '.--', '.--']

    wellids = ['{}{}'.format(w[0], w[1]) for w in wells]
    # Create curves

    # Any values greater than 2 or less than 0.05 will be None
    newdata = data
    newdata = maskOD(newdata)
    # Find high value for plotting
    hi = findHighOD(newdata) + 0.15
    for cidx, c in enumerate(sorted(newdata)):
        fmed, axarrmed = plt.subplots(8, 12, sharex=True, sharey=True)
        for idx, w in enumerate(wellids):
            row = int(py.floor(idx / 12))
            col = int(idx % 12)
            for ridx, rep in enumerate(sorted(newdata[c][w])):
                curve = newdata[c][w][rep].rawcurve
                if ridx == 0:
                    medarray = py.array([curve], ndmin=2)
                else:
                    medarray = py.concatenate((medarray, [curve]))
            mcurve = func(medarray, axis=0)
            clr = colors[3]
            shp = shapes[3]
            axarrmed[row, col].plot(time, mcurve, '{}{}'.format(clr, shp), linewidth=1.5)
            axarrmed[row, col].errorbar(time, mcurve,
                                        yerr=py.std(medarray, axis=0),
                                        linestyle='None',
                                        ecolor='#888888')
            axarrmed[row, col].tick_params(axis='both', left='on', right='off', bottom='on', top='off')
            axarrmed[row, col].set_ylim((0, hi))
#            axarrmed[row, col].set_xlim((0, time))
            axarrmed[row, col].set_title(w, fontweight='bold')
            #axarrmed[row, col].set_axis_bgcolor('#B0B0B0')
        fmed.set_size_inches(35, 18)
        axarrmed[7, 0].set_xlabel('Time (hr)', fontweight='bold')
        axarrmed[7, 0].set_ylabel(r'OD$_{600nm}$', fontweight='bold')
        #plt.legend(bbox_to_anchor=(1, 4.5), loc='center left', fancybox=True)
        plt.savefig('{}/{}_{}_growthcurves.png'.format(outDir, c, name), dpi=100, bbox_inches='tight')


def findHighOD(data):
    '''Find the highest OD value to set for plots'''
    hi = 0
    for c, wDict in data.items():
        for w, rDict in wDict.items():
            for r in rDict:
                max = py.amax(rDict[r].rawcurve)
                if max > hi:
                    hi = max
    return hi


def maskOD(data):
    '''Mask too large/small values for plots'''
    for c, wDict in data.items():
        for w, rDict in wDict.items():
            for r in rDict:
                rDict[r].rawcurve[(rDict[r].rawcurve > 2)
                                  | (rDict[r].rawcurve < 0.01)
                                  ] = None
                # TODO: Report masks when they occur
                if py.isnan(py.sum(rDict[r].rawcurve)):
                    print('Masking value with "nan" in {} -- {} -- {}'.format(c, w, r),
                          file=sys.stderr)
    return data

# PMData.py
# Phenotype microarray module for parsing optical density data
#
# Author: Daniel A Cuevas
# Created on 12 Dec. 2013
# Updated on 16 Jan. 2015

from __future__ import absolute_import, division, print_function
import pylab as py
import sys


class PMData:
    '''Class for parsing phenotype microarray data'''
    def __init__(self, filepath, plateFlag):
        self.filepath = filepath
        self.numClones = 0
        self.numConditions = 0
        self.numFiltered = 0
        self.plateFlag = plateFlag
        self.plateName = {}  # Hash of clone->plateName
        self.replicates = {}  # Hash of clone->[reps]
        self.clones = []  # Set of unique clone names
        if plateFlag:
            self.wells = {}  # Hash of well->(mainsource, condition)
        else:
            self.wells = set([]) # Set of wells
        self.time = []  # Array of time values

        # Primary data structure to access data
        self.dataHash = {}  # clone->rep #->well|->[ODs]
                            #                   |->filter

        self.__beginParse()

    def __beginParse(self):
        '''Initiate parsing on the given PM file'''
        f = open(self.filepath, 'r')
        # Begin iteration through file
        for lnum, l in enumerate(f):
            l = l.rstrip('\n')
            ll = l.split('\t')

            # Line 1: header
            if lnum == 0:
                self.__parseHeader(ll)

            # Line 2+: growth curves
            else:
                self.__parseODCurve(ll)
        f.close()

        # Check each growth curve is the same length
        self.__QACheck()

        # Set number of growth conditions present
        self.numConditions = len(self.wells)

    def __parseHeader(self, ll):
        '''Header line contains data columns and time values'''
        if self.plateFlag:
            self.time = py.array([float(x) for x in ll[5:]])
        else:
            self.time = py.array([float(x) for x in ll[2:]])

    def __parseODCurve(self, ll):
        '''Growth curve parsing method'''
        # Extract curve info
        if self.plateFlag:
            (c, ms, gc, pn, w) = ll[0:5]
            # Add well info
            self.wells[w] = (ms, gc)
            curve = py.array([float(x) for x in ll[5:]])
        else:
            (c, w) = ll[0:2]
            curve = py.array([float(x) for x in ll[2:]])
            self.wells.add(w)

        # Extract clone name and replicate name
        # If no replicate name exists, assign it "1"
        parsedName = c.split('_')
        try:
            (cName, rep) = parsedName[0:2]
        except ValueError as e:
            cName = parsedName[0]
            rep = "1"

        if cName not in self.clones:
            self.clones.append(cName)
            self.numClones += 1
            self.replicates[cName] = []

        if rep not in self.replicates[cName]:
            self.replicates[cName].append(rep)

        # Record plate name
        if self.plateFlag:
            self.plateName[cName] = pn

        # Add curve to primary data hash
        # Initialize data hash when needed
        try:
            self.dataHash[cName]
        except KeyError:
            self.dataHash[cName] = {}
        try:
            self.dataHash[cName][rep]
        except KeyError:
            self.dataHash[cName][rep] = {}

        self.dataHash[cName][rep][w] = {'od': curve, 'filter': False}

    def __QACheck(self):
        '''QA check to ensure stable data set'''
        problems = []
        numTime = len(self.time)
        for clone, repDict in self.dataHash.items():
            for rep, wellDict in repDict.items():
                for w, odDict in wellDict.items():
                    # Find number of values in growth curve
                    numVals = len(odDict['od'])
                    if numVals != numTime:
                        if self.plateFlag:
                            (ms, gc) = self.wells[w]
                            problems.append([clone, rep, ms, gc, w,
                                            'time:{}\tcurve:{}'.format(
                                                numTime, numVals)])

                        else:
                            problems.append([clone, rep, w,
                                            'time:{}\tcurve:{}'.format(
                                                numTime, numVals)])

        # Print out issues
        for p in problems:
            print('\t'.join([str(x) for x in p]), file=sys.stderr)

    def getCloneReplicates(self, clone, w, applyFilter=False):
        '''Retrieve all growth curves for a clone+well'''
        # Check if any other replicates should be returned
        # retArray is a 2xN multidimensional numpy array
        retArray = py.array([])
        first = True
        for rep in self.replicates[clone]:
            # Get replicate
            filterMe = self.dataHash[clone][rep][w]['filter']
            currCurve = self.dataHash[clone][rep][w]['od']

            # Check if filter is enabled and curve should be filtered
            if applyFilter and filterMe:
                continue

            # Create multidimensional array if first
            elif first:
                retArray = py.array([currCurve])
                first = False

            # Append to multidimensional array if not first
            else:
                retArray = py.concatenate((retArray,
                                           py.array([currCurve])))

        return retArray

    def getODCurve(self, clone, w, rep):
        '''Retrieve a single OD curve'''
        return self.dataHash[clone][rep][w]['od']

    def getFiltered(self):
        '''Retrieve array of all growth curves labeled as filtered'''
        # Create array of tuples for each replicate curve labeled as filtered
        # Format: [(clone, main source, growth condition,
        #           replicate #, [OD values]), (next...), ...]
        ret = []

        # Iterate through clones
        for clone, repDict in self.dataHash.items():
            # Iterate through replicates
            for rep, wellDict in repDict.items():
                # Iterate through wells
                for w, odDict in wellDict.items():
                    # Check if filter is set to True
                    if odDict['filter']:
                        ret.append((clone, rep, w, odDict['od']))

        return ret

    def setFilter(self, clone, rep, w, filter):
        '''Set filter for specific curve'''
        oldFilter = self.dataHash[clone][rep][w]['filter']
        self.dataHash[clone][rep][w]['filter'] = filter

        # If filter changed from False to True, increment number of filtered
        if not oldFilter and filter:
            self.numFiltered += 1

        # If filter changed from True to False, decrement number of filtered
        elif oldFilter and not filter:
            self.numFiltered -= 1

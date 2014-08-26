# PMData.py
# Phenotype microarray module for parsing optical density data
#
# Author: Daniel A Cuevas
# Created on 12 Dec. 2013
# Updated on 26 Aug. 2014

import pylab as py
import sys


class PMData:
    '''Class for parsing phenotype microarray data'''
    def __init__(self, filepath):
        self.filepath = filepath
        self.numClones = 0
        self.numConditions = 0
        self.numFiltered = 0
        self.replicates = {}  # Hash of clone->[reps]
        self.clones = []  # Set of unique clone names
        self.wells = {}  # Hash of well->(mainsource, condition)
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
        self.time = py.array([float(x) for x in ll[4:]])

    def __parseODCurve(self, ll):
        '''Growth curve parsing method'''
        # Extract curve info
        (c, ms, gc, w) = ll[0:4]

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

        # Add well info
        self.wells[w] = (ms, gc)

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

        curve = py.array([float(x) for x in ll[4:]])
        self.dataHash[cName][rep][w] = {'od': curve, 'filter': False}

    def __QACheck(self):
        '''QA check to ensure stable data set'''
        problems = []
        numTime = len(self.time)
        for clone, repDict in self.dataHash.items():
            for rep, wellDict in repDict.items():
                for w, odDict in wellDict.items():
                    (ms, gc) = self.wells[w]
                    # Find number of values in growth curve
                    numVals = len(odDict['od'])
                    if numVals != numTime:
                        problems.append([clone, rep, ms, gc, w,
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

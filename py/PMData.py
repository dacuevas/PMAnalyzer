# PMData.py
# Phenotype microarray module for parsing optical density data
#
# Author: Daniel A Cuevas
# Created on 12 Dec. 2013
# Updated on 05 Mar. 2014

import pylab as py
import sys


class PMData:
    '''Class for parsing phenotype microarray data'''
    def __init__(self, filepath):
        self.filepath = filepath
        self.numClones = 0
        self.numConditions = 0
        self.numFiltered = 0
        self.numReplicates = {}  # Hash of clone->{rep. count}
        self.clones = []  # Set of unique clone names
        self.conditions = {}  # Hash of source->[conditions]
        self.wells = {}  # Hash of source->{condition]->well
        self.time = []  # Array of time values

        self.clonesNU = []  # Array of clones (non-unique)
        self.sourcesNU = []  # Array of sources (non-unique)
        self.conditionsNU = []  # Array of conditions (non-unique)

        # Primary data structure to access data
        self.dataHash = {}  # clone->{rep #}|->{source}->{condition}->[ODs]
                            #                                      |->{filter}
        self.__beginParse()

    def __beginParse(self):
        '''Initiate parsing on the given PM file'''
        f = open(self.filepath, 'r')
        # Begin iteration through file
        for lnum, l in enumerate(f):
            l = l.rstrip('\n')
            ll = l.split('\t')

            # Line 1: clone names
            if lnum == 0:
                self.__parseClones(ll)

            # Line 2: source names
            elif lnum == 1:
                self.__parseSources(ll)

            # Line 3: condition substrates
            elif lnum == 2:
                self.__parseConditions(ll)

            # Line 4: well indicies [A-H][1-12]
            elif lnum == 3:
                self.__parseWells(ll)

            # Line 5+: OD values
            else:
                self.__parseOD(ll)
        f.close()

        # Check each growth curve is the same length
        self.__QACheck()

    def __parseClones(self, ll):
        '''Clone line parsing method'''
        # All non-unique clones (order preserved)
        self.clonesNU = ll[1:]

        # Unique set of clones
        self.clones = set(ll[1:])
        self.numClones = len(self.clones)

    def __parseSources(self, ll):
        '''Main sources line parsing method'''
        # All non-unique main sources (order preserved)
        self.sourcesNU = ll[1:]

        # Initialize conditions hash
        self.conditions = {s: [] for s in set(self.sourcesNU)}

    def __parseConditions(self, ll):
        '''Growth conditions line parsing method'''
        # All non-unique growth conditions (order preserved)
        self.conditionsNU = ll[1:]

        # Add to conditions hash
        [self.conditions[self.sourcesNU[idx]].append(c)
         for idx, c in enumerate(self.conditionsNU)]

        # Duplicate conditions created for each source
        # in above method - must remove for unique set
        for source in self.conditions:
            self.conditions[source] = set(self.conditions[source])
            self.numConditions += len(self.conditions[source])

        # Initialize main data hash
        prevClone = ""
        prevCond = ""
        numRep = 1
        if self.numClones == 1:
            totalReps = len(self.conditionsNU) / self.numConditions

        for idx, clone in enumerate(self.clonesNU):
            currCond = self.conditionsNU[idx]
            # Keep track of replicate numbers
            # Reset rep number if we reach a new clone name
            if self.numClones == 1:
                numRep = (idx % totalReps) + 1

            else:
                numRep = numRep + 1 if (clone == prevClone and
                                        currCond == prevCond) else 1

            # Update replicate count for clone
            self.numReplicates[clone] = numRep
            prevClone = clone
            prevCond = currCond

            if clone not in self.dataHash:
                self.dataHash[clone] = {}

            if numRep not in self.dataHash[clone]:
                self.dataHash[clone][numRep] = {}

            # Add condition to main data hash
            # Pre-set filter to False
            for source, sourceList in self.conditions.items():
                self.dataHash[clone][numRep][source] =\
                    {cond: {'filter': False, 'od': py.array([])}
                     for cond in sourceList}

    def __parseWells(self, ll):
        '''Well line parsing method'''
        # Store as source->{condition}->well
        for idx, well in enumerate(ll[1:]):
            source = self.sourcesNU[idx]
            cond = self.conditionsNU[idx]
            try:
                # Intialize hash if the source was not yet added
                self.wells[source]
            except KeyError:
                self.wells[source] = {}
            self.wells[source][cond] = well

    def __parseOD(self, ll):
        '''OD data lines parsing method'''
        ll = [float(x) for x in ll]

        # Add the current time
        self.time.append(ll[0])
        numRep = 1
        prevClone = ""
        prevCond = ""
        if self.numClones == 1:
            totalReps = len(self.conditionsNU) / self.numConditions

        for idx, od in enumerate(ll[1:]):
            clone = self.clonesNU[idx]
            source = self.sourcesNU[idx]
            condition = self.conditionsNU[idx]

            # Check which clone + replicate we are observing
            if self.numClones == 1:
                numRep = (idx % totalReps) + 1

            else:
                numRep = numRep + 1 if (clone == prevClone and
                                        condition == prevCond) else 1
            prevClone = clone
            prevCond = condition

            # Append OD reading to array
            self.dataHash[clone][numRep][source][condition]['od'] =\
                py.append(self.dataHash[clone][numRep][source]
                          [condition]['od'], od)

    def __QACheck(self):
        '''QA check to ensure stable data set'''
        problems = []
        numTime = len(self.time)
        for clone, repDict in self.dataHash.items():
            for rep, sourceDict in repDict.items():
                for source, condDict in sourceDict.items():
                    for cond, odDict in condDict.items():

                        # Find number of values in growth curve
                        numVals = len(odDict['od'])
                        if numVals != numTime:
                            problems.append([clone, rep, source, cond,
                                             'time:{}\tcurve:{}'.format(
                                                 numTime, numVals)])

        # Print out issues
        for p in problems:
            print >> sys.stderr, '\t'.join([str(x) for x in p])

    def getCloneReplicates(self, clone, source, condition, applyFilter=False):
        '''Retrieve all growth curves for a clone+source+condition'''
        # Check if any other replicates should be returned
        # retArray is a 2xN multidimensional numpy array
        retArray = py.array([])
        first = True
        for i in xrange(1, self.numReplicates[clone] + 1):
            # Get replicate
            filterMe = self.dataHash[clone][i][source][condition]['filter']
            currCurve = self.dataHash[clone][i][source][condition]['od']

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

    def getFiltered(self):
        '''Retrieve array of all growth curves labeled as filtered'''
        # Create array of tuples for each replicate curve labeled as filtered
        # Format: [(clone, main source, growth condition,
        #           replicate #, [OD values]), (next...), ...]
        ret = []

        # Iterate through clones
        for clone, repDict in self.dataHash.items():

            # Iterate through replicates
            for rep, sourceDict in repDict.items():

                # Iterate through main sources
                for source, condDict in sourceDict.items():

                    # Iterate through growth conditions
                    for cond, odDict in condDict.items():

                        # Check if filter is set to True
                        if odDict['filter']:
                            ret.append((clone, source, cond, rep,
                                        odDict['od']))

        return ret

    def setFilter(self, clone, rep, source, condition, filter):
        '''Set filter for specific curve'''
        oldFilter = self.dataHash[clone][rep][source][condition]['filter']
        self.dataHash[clone][rep][source][condition]['filter'] = filter

        # If filter changed from False to True, increment number of filtered
        if not oldFilter and filter:
            self.numFiltered += 1

        # If filter changed from True to False, decrement number of filtered
        elif oldFilter and not filter:
            self.numFiltered -= 1

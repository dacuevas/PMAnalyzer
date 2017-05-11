# PMData.py
# Phenotype microarray module for parsing optical density data
#
# Author: Daniel A Cuevas
# Created on 12 Dec 2013
# Updated on 11 May 2017

from __future__ import absolute_import, division, print_function
import pandas as pd


class PMData:
    """Class for parsing phenotype microarray data"""
    def __init__(self, filepath, plateFlag):
        self.plateFlag = plateFlag
        self.replicates = {}  # Hash of clone->[reps]
        self.clones = []  # Set of unique clone names
        self.wells = []

        # Primary data structure to access data
        self.DF = pd.DataFrame()

        self.__loadData(filepath)
        self.__init()

    def __loadData(self, filepath):
        """Load data into Pandas DataFrame"""
        indices = ["sample", "rep", "well", "time"]
        self.DF = pd.read_csv(filepath, delimiter="\t", index_col=indices,
                              dtype={"sample": str, "rep": str})
        self.DF = self.DF.sort_index(level=[0, 1, 2, 3])

    def __init(self):
        """Initialize all class variables"""
        self.__sortWells()

    def __sortWells(self):
        """Sort wells numerically rather than alphanumerically"""
        #self.wells = [(x[0], int(x[1:]))
        #              for x in self.DF.index.levels[2]]
        #self.wells = sorted(self.wells, key=itemgetter(0, 1))
        #self.wells = ["{}{}".format(x[0], x[1]) for x in self.wells]
        self.wells = self.DF.index.get_level_values("well").unique()

    def __QACheck(self):
        """QA check to ensure stable data set"""
        pass

    def getSampleNames(self):
        return self.DF.index.get_level_values("sample").unique()

    def getNumSamples(self):
        return len(self.getSampleNames())

    def getReplicates(self, sample):
        return self.DF.loc[sample].index.get_level_values("rep")

    def getWells(self):
        """
        Return well information
        With plate info: return a DataFrame
        Without plate info: return an Index array
        """
        if self.plateFlag:
            # Grab only the mainsource and compound columns
            wells = self.DF[["mainsource", "compound"]]

            # Remove duplicate items by grouping by the well id and
            wells = wells.groupby(level="well", sort=False).last()
        else:
            wells = self.wells
        return wells

    def getNumWells(self):
        return len(self.getWells())

    def getODCurve(self, sample, well, rep):
        """Retrieve a single OD curve"""
        return self.DF.loc[(sample, rep, well, slice(None))]["od"]

    def getMedianCurves(self):
        """Return DataFrame median curves for each sample"""
        df = self.DF.median(level=["sample", "well", "time"])
        if self.plateFlag:
            leftMerge = df.reset_index()
            rightMerge = self.DF[["mainsource", "compound"]].reset_index(
                level=[1, 3], drop=True).reset_index().drop_duplicates()
            df = pd.merge(
                leftMerge,
                rightMerge,
                on=["sample", "well"],
                how="left"
            ).set_index(["sample", "well"])
        return df

    def getMeanCurves(self):
        """Return DataFrame of mean curves for each sample"""
        df = self.DF.mean(level=["sample", "well", "time"])
        if self.plateFlag:
            leftMerge = df.reset_index()
            rightMerge = self.DF[["mainsource", "compound"]].reset_index(
                level=[1, 3], drop=True).reset_index().drop_duplicates()
            df = pd.merge(
                leftMerge,
                rightMerge,
                on=["sample", "well"],
                how="left"
            ).set_index(["sample", "well"])
        return df

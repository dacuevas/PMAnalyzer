# PMUtil.py
# Phenotype microarray utility functions
#
# Author: Daniel A Cuevas
# Created on 27 Jan. 2015
# Updated on 27 Jan. 2015

from __future__ import absolute_import, division, print_function
import sys
import time
import datetime


def timeStamp():
    '''Return time stamp'''
    t = time.time()
    fmt = '[%Y-%m-%d %H:%M:%S]'
    return datetime.datetime.fromtimestamp(t).strftime(fmt)


def printStatus(msg):
    '''Print status message'''
    print('{} {}'.format(timeStamp(), msg), file=sys.stderr)
    sys.stderr.flush()

#
# For the solar flare finder paper.  Calculates probability
# mass distributions as a function of the instrument(s)
# selected based on their individual flare detection
# probabilities.
#
import os
import copy
import itertools
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scipy.stats import binom
from instruments import xrt, megsa, megsb, eis, sot, rhessi, iris

# Image root
imgroot = os.path.expanduser('~/solar_flare_finder/img')

# CSV root
csvroot = os.path.expanduser('~/solar_flare_finder/csv')

# Number of flares
nf = 934


# Simple storage
class SolarFlareFinder(object):
    def __init__(self, instrument, fdp):
        self.instrument = instrument
        self.fdp = fdp


class FlareDetectionProbability(object):
    def __init__(self, timerange, p):
        self.timerange = timerange
        self.p = p

timerange = [None, None]


summary_names = ['mean', 'median', 'mode', 'lower68', 'higher68', 'n']

# Numbers taken from the paper
xrt_sff = SolarFlareFinder(xrt, FlareDetectionProbability(timerange, 0.57))
megsa_sff = SolarFlareFinder(megsa, FlareDetectionProbability(timerange, 1.0))
megsb_sff = SolarFlareFinder(megsb, FlareDetectionProbability(timerange, 0.12))
eis_sff = SolarFlareFinder(eis, FlareDetectionProbability(timerange, 0.06))
sot_sff = SolarFlareFinder(sot, FlareDetectionProbability(timerange, 0.13))
rhessi_sff = SolarFlareFinder(rhessi, FlareDetectionProbability(timerange, 0.58))
iris_sff = SolarFlareFinder(iris, FlareDetectionProbability(timerange, 0.11))

# Instruments we want to consider together
all_sff = (xrt_sff, megsa_sff, megsb_sff, eis_sff, sot_sff, rhessi_sff, iris_sff)

# Generate tabular data header
all_sff_names = [sff.instrument.name for sff in all_sff]
header = copy.deepcopy(all_sff_names)
for summary_name in summary_names:
    header.append(summary_name)

# Create the table
t = Table(names=header)

# Calculate some properties of the PMF for multiple instruments.

# Go through all the number of instruments
for i in range(1, len(all_sff)+1):

    # Calculate all the permutations of 'i' instrumens
    joints = itertools.combinations(all_sff, i)
    for joint in joints:

        # Calculate the joint name and the joint probability
        p_joint = 1.0
        instr_joint = []
        joint_names = []
        for instr in joint:
            p_joint = p_joint * instr.fdp.p
            instr_joint.append(instr.instrument.name)

        # Calculate the binomial distribution
        x = np.arange(0, nf+1)
        b = binom(nf, p_joint)
        y = b.pmf(x)

        # Calculate summaries of the distribution
        mode = x[np.argmax(y)]
        mean = np.float64(b.stats(moments='m'))
        interval = b.interval(alpha=0.6827)
        big_interval = b.interval(alpha=0.9999)
        median = np.int(b.median())

        summary_stats = {'mean': mean, 'median': median, 'mode': mode,
                         'lower68': interval[0], 'higher68': interval[1],
                         'n': i}

        # Make the plot and write it out
        instr_joint = ", ".join(instr_joint)
        xlim = [0.99*big_interval[0], np.max([10, big_interval[1]])]
        plt.close('all')
        plt.ion()
        plt.plot(x, y, label='B({}, {:.4f})'.format(nf, p_joint))
        plt.axvline(mode, color='k', label='most probable ({})'.format(mode))
        plt.axvline(mean, color='k', linestyle='dotted', label='mean ({:.4f})'.format(mean))
        plt.axvline(median, color='k', linestyle='dashed', label='median ({})'.format(median))
        interval68 = '(${}'.format(interval[0]) + r'\rightarrow' + '{}$)'.format(interval[1])
        plt.axvline(interval[0], color='r', linestyle='dashed', label='68% interval {}'.format(interval68))
        plt.axvline(interval[1], color='r', linestyle='dashed')
        plt.xlim(xlim)
        plt.xlabel('number of flares')
        plt.ylabel('probability')
        plt.title('{:s} probability mass function'.format(instr_joint))
        plt.legend()
        filename = 'binomial_i{:n}_{:s}'.format(i, instr_joint)
        filepath = '{}/{}'.format(imgroot, filename)
        plt.savefig(filepath)

        # Update the table
        results_list = []
        for sff_name in all_sff_names:
            if sff_name in instr_joint:
                results_list.append(1)
            else:
                results_list.append(0)

        for summary_name in summary_names:
            results_list.append(summary_stats[summary_name])
        t.add_row(results_list)

# Write the table out
filename = 'summary_stats.csv'
filepath = '{}/{}'.format(csvroot, filename)
t.write(filepath)

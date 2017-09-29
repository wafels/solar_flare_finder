#
# For the solar flare finder paper.  Calculates probability
# mass distributions as a function of the instrument(s)
# selected based on their individual flare detection
# probabilities.
#
import os
import itertools
import matplotlib.pyplot as plt
import numpy as np
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
joint = (rhessi_sff, xrt_sff, eis_sff)

# Calculate some properties of the PMF for multiple
# instruments.

for i in range(1, len(all_sff)+1):
    joints = itertools.permutations(all_sff, i)
    for joint in joints:

        p_joint = 1.0
        instr_joint = []
        for instr in joint:
            p_joint = p_joint * instr.fdp.p
            instr_joint.append(instr.instrument.name)

        instr_joint = ", ".join(instr_joint)
        x = np.arange(0, nf+1)
        b = binom(nf, p_joint)
        y = b.pmf(x)
        most_probable = x[np.argmax(y)]

        mean = b.stats(moments='m')
        interval = b.interval(alpha=0.6827)
        median = b.median()

        # Make the plot
        plt.close('all')
        plt.ion()
        plt.plot(x, y, label='B({}, {})'.format(nf, p_joint))
        plt.axvline(most_probable, color='k', label='most probable ({})'.format(most_probable))
        plt.axvline(mean, color='k', linestyle='dotted', label='mean ({})'.format(mean))
        plt.axvline(median, color='k', linestyle='dashed', label='median ({})'.format(median))
        interval68 = '(${}'.format(interval[0]) + r'\rightarrow' + '{}$)'.format(interval[1])
        plt.axvline(interval[0], color='r', linestyle='dashed', label='68% interval {}'.format(interval68))
        plt.axvline(interval[1], color='r', linestyle='dashed')
        plt.xlabel('number of flares')
        plt.ylabel('probability')
        plt.title('{:s} probability mass function'.format(instr_joint))
        plt.legend()
        filename = 'binomial_i{:n}_{:s}'.format(i, instr_joint)
        filepath = '{}{}'.format(imgroot, filename)
        plt.savefig(filepath)


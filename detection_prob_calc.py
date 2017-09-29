#
# For the solar flare finder paper.  Calculates probability
# mass distributions as a function of the instrument(s)
# selected based on their individual flare detection
# probabilities.
#
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom
from instruments import xrt, megsa, megsb, eis, sot, rhessi, iris

# Number of flares
nf = 934


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
joint = (iris_sff, rhessi_sff, xrt_sff)

# Calculate some properties of the PMF for multiple
# instruments.
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

mean, var, skew, kurt = binom.stats(nf, p_joint, moments='mvsk')


# Make the plot
plt.close('all')
plt.ion()
plt.plot(x, y, label='binomial({}, {})'.format(nf, p_joint))
plt.axvline(most_probable, color='k', linestyle=':', label='most probable ({})'.format(most_probable))
plt.xlabel('number of flares')
plt.ylabel('probability')
plt.title('{:s} probability mass function'.format(instr_joint))
plt.legend()
plt.show()


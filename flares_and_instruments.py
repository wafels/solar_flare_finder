from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from sunpy.time import parse_time

# Friday 20 January 2017 flare numbers
numbers_20170120 = {"num": np.asarray([338, 1394, 2557, 1903, 779, 541, 28, 2]),
                    "date": "2017/01/20"}

# Analysis of the flare numbers
num = numbers_20170120["num"]
frac = num / np.sum(num)
n = np.arange(0, len(frac))
mu = np.sum(n*frac)
poi = poisson(mu)

# Plot the measured distribution and the estimate distribuion
linewidth = 3
plt.ion()
plt.ylabel('fraction')
plt.xlabel('number of instruments')
plt.plot(n, frac, label='fraction observed', linewidth=linewidth)
mu_string = '$\hat\mu$'
mu_label = '{:s}={:3.1f}'.format(mu_string, mu)
plt.plot(n, poi.pmf(n), color='k', label='Poisson distribution ({:s})'.format(mu_label), linewidth=linewidth)
plt.axvline(mu, linestyle=':', color='k', label='mean ({:s}) number of instruments'.format(mu_string), linewidth=linewidth)
plt.legend(framealpha=0.5)

# Define the instrument class
class Instrument:
    def __init__(self, duty_cycle, disk_coverage, name,
                 observation_start=None, observation_end=None):
        self.duty_cycle = duty_cycle
        self.disk_coverage = disk_coverage
        self.name = name
        self.observation_start = observation_start
        self.observation_end = observation_end

# Define each of the instruments and the full instrument list
rhessi = Instrument([0.50, 0.50], [1.0, 1.0], 'RHESSI')
eve_a = Instrument([1.0, 1.0], [1.0, 1.0], 'MEGS-A')
eve_b = Instrument([1.0, 1.0], [1.0, 1.0], 'MEGS-B')
eis = Instrument([0.1, 1.0], [0.1, 1.0], 'EIS')
sot = Instrument([0.1, 1.0], [0.1, 1.0], 'SOT')
xrt = Instrument([0.1, 1.0], [0.1, 1.0], 'XRT')
iris = Instrument([0.1, 1.0], [0.1, 1.0], 'IRIS')

instruments = [rhessi, eve_a, eve_b, eis, sot, xrt, iris]

# Calculate estimates of the co-observation probability
lower_limit = upper_limit = midpoint = 1.0
for instrument in instruments:
    lower_limit = lower_limit * instrument.disk_coverage[0] * instrument.duty_cycle[0]
    upper_limit = upper_limit * instrument.disk_coverage[1] * instrument.duty_cycle[1]
    midpoint = midpoint * 0.5 * (instrument.disk_coverage[0] + instrument.disk_coverage[1]) * \
        0.5 * (instrument.duty_cycle[0] + instrument.duty_cycle[1])

# Count how many instruments were observing each day.
study_start_date = parse_time('2010/05/01')
study_end_date = parse_time('2016/10/31')
number_observing = []
this_date = study_start_date
while this_date <= study_end_date:
    n = 0
    for instrument in instruments:
        if (this_date >= instrument.observation_start) and (this_date <= instrument.observation_end):
            n += 1
    number_observing.append(n)
    this_date = this_date + timedelta(days=1)
number_observing = np.asarray(number_observing)

# Save the number of instruments observing per day as a function of time

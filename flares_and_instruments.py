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
plt.figure(1)
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
    def __init__(self, duty_cycle, disk_coverage, name, operation=None):
        self.duty_cycle = duty_cycle
        self.disk_coverage = disk_coverage
        self.name = name
        default_operation_end = parse_time('2099/12/31')
        if operation is not None:
            self.operation = [parse_time(operation[0])]
            if operation[1] is None:
                self.operation.append(default_operation_end)
            else:
                self.operation.append(parse_time(operation[1]))

# Define each of the instruments and the full instrument list
rhessi = Instrument([0.50, 0.50], [1.0, 1.0], 'RHESSI', operation=['2002/02/05', None])
eve_a = Instrument([1.0, 1.0], [1.0, 1.0], 'MEGS-A', operation=['2010/05/01', '2014/05/26'])
eve_b = Instrument([1.0, 1.0], [1.0, 1.0], 'MEGS-B', operation=['2010/05/01', None])
eis = Instrument([0.1, 1.0], [0.1, 1.0], 'EIS', operation=['2006/09/22', None])
sot = Instrument([0.1, 1.0], [0.1, 1.0], 'SOT', operation=['2006/09/22', None])
xrt = Instrument([0.1, 1.0], [0.1, 1.0], 'XRT', operation=['2006/09/22', None])
iris = Instrument([0.1, 1.0], [0.1, 1.0], 'IRIS', operation=['2013/06/27', None])

instruments = [rhessi, eve_a, eve_b, eis, sot, xrt, iris]

# Calculate estimates of the co-observation probability
lower_limit = upper_limit = midpoint = 1.0
for instrument in instruments:
    lower_limit = lower_limit * instrument.disk_coverage[0] * instrument.duty_cycle[0]
    upper_limit = upper_limit * instrument.disk_coverage[1] * instrument.duty_cycle[1]
    midpoint = midpoint * 0.5 * (instrument.disk_coverage[0] + instrument.disk_coverage[1]) * \
        0.5 * (instrument.duty_cycle[0] + instrument.duty_cycle[1])

# Study start and end date
study_start_date = parse_time('2010/05/01')
study_end_date = parse_time('2016/10/31')

# Count how many instruments were observing each day.
instrument_start_date = parse_time('2010/05/01')
instrument_start_date = parse_time('2001/05/01')
instrument_end_date = parse_time('2016/10/31')

number_observing = []
date_no = []
this_date = instrument_start_date
while this_date <= instrument_end_date:
    n = 0
    for instrument in instruments:
        operation = instrument.operation
        if (this_date >= operation[0]) and (this_date <= operation[1]):
            n += 1
    number_observing.append(n)
    date_no.append(this_date)
    this_date = this_date + timedelta(days=1)
number_observing = np.asarray(number_observing)

# Plot the number of instruments that were observing as a function of time
plt.figure(2)
plt.plot(date_no, number_observing, linewidth=linewidth)
plt.xlabel('observation time (year)')
plt.ylabel('number of instruments available to observe flares')
plt.ylim(np.max([0, np.min(number_observing)-1]), np.max(number_observing)+1)
plt.xlim(instrument_start_date, parse_time('now'))
plt.axvline(study_start_date, linestyle=":", color='k', label='study start/end time')
plt.axvline(study_end_date, linestyle=":", color='k')
plt.legend(loc='upper left')



# Save the number of instruments observing per day as a function of time

from sunpy.time import parse_time


# Define the instrument class
class Instrument(object):
    def __init__(self, observatory=None, instrument=None, detector=None, name=None,
                 duty_cycle=None, disk_coverage=None, operation=None):
        self.duty_cycle = duty_cycle
        self.disk_coverage = disk_coverage
        self.observatory = observatory
        self.instrument = instrument
        self.detector = detector
        self.name = name

        self._default_operation_end = parse_time('2099/12/31')
        if operation is not None:
            self.operation = [parse_time(operation[0])]
            if operation[1] is None:
                self.operation.append(self._default_operation_end)
            else:
                self.operation.append(parse_time(operation[1]))

# Define each of the instruments and the full instrument list
rhessi = Instrument(name='RHESSI',
                    duty_cycle=[0.50, 0.50],
                    disk_coverage=[1.0, 1.0],
                    operation=['2002/02/05', None])

megsa = Instrument(name='MEGS-A',
                   duty_cycle=[1.0, 1.0],
                   disk_coverage=[1.0, 1.0],
                   operation=['2010/05/01', '2014/05/26'])

megsb = Instrument(name='MEGS-B',
                   duty_cycle=[1.0, 1.0],
                   disk_coverage=[1.0, 1.0],
                   operation=['2010/05/01', None])

eis = Instrument(name='EIS',
                 duty_cycle=[0.1, 1.0],
                 disk_coverage=[0.1, 1.0],
                 operation=['2006/09/22', None])

sot = Instrument(name='SOT',
                 duty_cycle=[0.1, 1.0],
                 disk_coverage=[0.1, 1.0],
                 operation=['2006/09/22', None])

xrt = Instrument(name='XRT',
                 duty_cycle=[0.1, 1.0],
                 disk_coverage=[0.1, 1.0],
                 operation=['2006/09/22', None])

iris = Instrument(name='IRIS',
                  duty_cycle=[0.1, 1.0],
                  disk_coverage=[0.1, 1.0],
                  operation=['2013/06/27', None])

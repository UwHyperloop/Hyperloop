import numpy as np

from openmdao.core.component import Component

class Mission(Component):
    '''Place holder for real mission analysis. Could consider a pseudospectral optimal control approach'''
    def __init__(self):
        self.add_param('max_velocity', 308.0, desc='Maximum travel speed for pod', units='m/s')
        self.add_param('tube_len', 563270.0, desc='length of one trip', units='m')
        self.add_param('pwr_marg', 0.3, desc='fractional extra energy requirement')
        self.add_param('pwr_req', 420.0, desc='average power requirement for a mission', units='kW')

        self.add_output('time_mission', 0.0, desc='travel time to make one trip', units='s')
        self.add_output('energy', 0.0, desc='total energy storage requirement', units='kW*h')

    def solve_nonlinear(self, params, unknowns, resids):
        '''a VERY coarse approximation that takes the speed profile given in the original proposal as given. It just serves as a place holder for now. It's better than nothing, but a real analysis is needed here'''
        t1 = (300 - 0) * 1609 / (9.81 * 0.5) / 3600 # time needed to accelerate
        t2 = (555 - 300) * 1609 / (9.81 * 0.5) / 3600 # time needed to accelerate
        t3 = (params['max_velocity'] - 555) * 1609 / (9.81 * 0.5) / 3600 # time needed to accelerate
        #speed profile data from hyperloop alpha proposal, pg43
        dataStart = np.array([
            [0, 0],
            [t1, 300 * 1.609], # (300[mi/h] * 1609[m/mi]) / (9.81[m/s] * 0.5) / 3600[s/h]
            [167, 300 * 1.609],
            [167 + t2, 555 * 1.609], # (555 - 300[mi/h] * 1609[m/mi]) / (9.81[m/s] * 0.5) / 3600[s/h]
            [435, 555 * 1.609],
            [435 + t3, params['max_velocity']]])
        startUp = np.trapz(dataStart[:, 1] / 3600, dataStart[:, 0]) * 1000 # km covered during start up Los Angeles Grapevine
        dataEnd = np.array([
            [0, params['max_velocity']],
            [t3, 555 * 1.609],
            [t3 + 100, 555 * 1.609],
            [t3 + 100 + t2, 300 * 1.609],
            [t3 + 100 + t2 + 400, 300 * 1.609],
            [t3 + 100 + t2 + 400 + t1, 0]])
        windDown = np.trapz(dataEnd[:, 1] / 3600, dataEnd[:, 0]) * 1000 # km covered during wind down along I-580 to SF      
        len_middle = params['tube_length'] - (startUp + windDown)
        time_middle = middleLength / params['max_velocity']
        unknowns['time_mission'] = time_middle + 435 + t3 + t3 + 100 + t2 + 400 + t1

        unknowns['energy'] = (params['pwr_req'] * unknowns['time_mission'] / 3600.0) * (1 + params['pwr_marg']) # convert to hours

class SubscaleMission(Component):
    '''Place holder for real mission analysis. Could consider a pseudospectral optimal control approach'''
    def __init__(self):
        self.add_param('launch_v', 90.0, desc='maximum travel speed for pod', units='m/s')
        self.add_param('launch_time', 5.0, desc='time spent accelerating', units='s')
        self.add_param('tube_len', 1600.0, desc='length of one trip', units='m')
        self.add_param('pwr_marg', 0.3, desc='fractional extra energy requirement')
        self.add_param('pwr_req', 420.0, desc='average power requirement for a mission', units='kW')

        self.add_output('time_mission', 0.0, desc='travel time to make one trip', units='s')
        self.add_output('energy', 0.0, desc='total energy storage requirement', units='kW*h')


if __name__ == "__main__":
    from openmdao.core.problem import Problem
    from openmdao.core.group import Group

    p = Problem(root=Group())
    p.root.add('comp', Mission())
    p.setup()
    p.run()

from math import pi

from openmdao.core.component import Component

THICKNESS_RATIO = 0.23 / 111.5 # ratio given on pg27 of original proposal

class TubeStructural(Component):
    '''Placeholder for real structural calculations to size the tube wall Thickness''' 
    def __init__(self):
        super(TubeStructural, self).__init__()
        self.add_param('tube_P', 99.0, desc='static pressure in tube', units='Pa')
        self.add_param('tube_T', 292.1, desc='static temperature in tube', units='degK')
        self.add_param('tube_r', 0.9, desc='inner radius of tube', units='m')
        self.add_param('fill_area', 0.214, desc='cross sectional area filled with solid e.g. concrete floor', units='m**2')

        self.add_output('tube_r_outer', 1.0, desc='outer radius of tube', units='m')
        self.add_output('tube_area', 2.33, desc='cross sectional area inside of tube', units='m**2')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['tube_r_outer'] = params['tube_r'] * (1.0 + THICKNESS_RATIO)
        unknowns['tube_area'] = pi * params['tube_r'] ** 2 - params['fill_area']

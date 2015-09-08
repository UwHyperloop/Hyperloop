from openmdao.core.component import Component

THICKNESS_RATIO = 0.23 / 111.5 # ratio given on pg27 of original proposal

class TubeStructural(Component):
    '''Placeholder for real structural calculations to size the tube wall Thickness''' 
    def __init__(self):
        super(TubeStructural, self).__init__()
        self.add_param('Ps', 99.0, desc='static pressure in the tube', units='Pa')
        self.add_param('r_inner', 3.0, desc='inner radius of tube', units='m')

        self.add_output('r_outer', 3.006, desc='outer radius of tube', units='m')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['r_outer'] = params['r_inner'] * (1.0 + THICKNESS_RATIO)

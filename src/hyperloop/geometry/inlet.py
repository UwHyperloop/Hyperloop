from math import pi, sqrt

from openmdao.core.component import Component

class InletGeom(Component):
    '''Calculates the dimensions for the inlet and compressor entrance'''
    def __init__(self):
        super(InletGeom, self).__init__()
        self.add_param('wall_thickness', 0.05, desc='thickness of inlet wall', units='m')
#        self.add_param('area_in', 0.0, desc='flow area required at front of inlet', units='m**2')
        self.add_param('area_out', 0.0, desc='flow area required at back of inlet', units='m**2')
        self.add_param('hub_to_tip', 0.4, desc='hub to tip ratio for compressor')
        self.add_param('cross_section', 1.4, desc='cross sectional area of passenger capsule', units='m**2')
        self.add_param('tube_area', 2.33, desc='cross sectional area inside of tube', units='m**2')

        self.add_output('r_back_inner', 0.0, desc='inner radius of back of inlet', units='m')
        self.add_output('r_back_outer', 0.0, desc='outer radius of back of inlet', units='m')
        self.add_output('bypass_area', 0.0, desc='available flow area round capsule', units='m**2')
        self.add_output('area_frontal', 0.0, desc='total capsule frontal area', units='m**2')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['r_back_inner'] = sqrt(params['area_out'] / pi / (1.0 - params['hub_to_tip'] ** 2))
        unknowns['r_back_outer'] = unknowns['r_back_inner'] + params['wall_thickness']
        unknowns['bypass_area'] = params['tube_area'] - params['cross_section']
        unknowns['area_frontal'] = pi * (unknowns['r_back_outer']) ** 2

if __name__ == '__main__':
    from openmdao.core.problem import Problem
    from openmdao.core.group import Group

    p = Problem(root=Group())
    p.root.add('comp', InletGeom())
    p.setup()
    p.run()

    for var_name, units in (('r_back_inner', 'm'), ('r_back_outer', 'm'), ('bypass_area', 'm**2'), ('area_frontal', 'm**2')):
        print '%s (%s): %f' % (var_name, units, p.root.comp.unknowns[var_name])

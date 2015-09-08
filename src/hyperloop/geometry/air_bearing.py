from math import pi, sin

from openmdao.core.component import Component

class AirBearing(Component):
    def __init__(self):
        super(AirBearing, self).__init__()
        self.add_param('tube_radius', 4.0, desc='radius of tube', units='m')
        self.add_param('capsule_mass', 15000.0, desc='mass of capsule', units='kg')
        self.add_param('P_in', 9.4, desc='air injection pressure for bearings', units='kPa')
        self.add_param('n_bearings', 7, desc='number of rows of bearing pads')
        self.add_param('sweep_angle', 4.0, desc='sweep angle of a single pad on tube wall', units='deg')
        self.add_param('eff', 0.5, desc='efficiency of the bearings (load capacity / grounding force)')

        self.add_output('total_area', 0.0, desc='total required bearing area', units='m**2')
        self.add_output('bearing_area', 0.0, desc='reqiured area per bearing', units='m**2')
        self.add_output('bearing_len', 0.0, desc='required length per bearing', units='m')
        self.add_output('bearing_width', 0.0, desc='linear width of bearing', units='m')

    def solve_nonlinear(self, params, unknowns, resids):
        arc_len = params['tube_radius'] * params['sweep_angle'] * pi / 180.0
        unknowns['total_area'] = params['capsule_mass'] * 9.81 / (params['P_in'] * 1000.0) / params['eff'] #convert to Pa from kPa
        unknowns['bearing_width'] = 2.0 * params['tube_radius'] * sin(params['sweep_angle'] * pi / 180.0 / 2.0)
        total_len = unknowns['total_area'] / arc_length / 2.0 # divide by 2 because there are 2 parallel skis
        unknowns['bearing_len'] = total_len / params['n_bearings'] / 2.0
        unknowns['bearing_area'] = unknowns['total_area'] / params['n_bearings'] / 2.0

if __name__ == '__main__':
    from openmdao.core.problem import Problem
    from openmdao.core.group import Group

    p = Problem(root=Group())
    p.root.add('comp', AirBearing())
    p.setup()
    p.run()

    print 'total_area (m**2): %f' % p.root.comp.unknowns['total_area']
    print 'area_per_bearing (m**2): %f' % p.root.comp.unknowns['area_per_bearing']
    print 'length_per_bearing (m): %f' % p.root.comp.unknowns['length_per_bearing']
    print 'bearing_width (m): %f' % p.root.comp.unknowns['bearing_width']

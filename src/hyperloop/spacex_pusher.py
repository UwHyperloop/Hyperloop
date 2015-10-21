from openmdao.core.component import Component

import math

class Pusher(Component):
    """Calculates rough final launch velocity of a sub-scale pod launched with the SpaceX wheeled pusher, based on numbers provided in October 2015 Tube Spec draft and assuming a logarithmic relationship between pod mass and acceleration"""

    def __init__(self):
        super(Pusher, self).__init__()

        self.add_param('g', 9.807, desc='gravity constant', units='m/s')
        self.add_param('displacement', 243.8, desc='maximum acceleration distance', units='m')
        self.add_param('pod_mass', 500.0, desc='pod mass including payload', units='kg')
        self.add_param('A', -0.651*9.807, desc='A in acceleration=A*ln(pod mass)+B')
        self.add_param('B', 6.0025*9.807, desc='B in acceleration=A*ln(pod mass)+B')

        self.add_output('pod_a', 0.0, desc='acceleration of pod', units='m/s**2')
        self.add_output('t', 0.0, desc='time required to complete acceleration', units='s')
        self.add_output('pod_V', 0.0, desc='final velocity after acceleration', units='m/s')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['pod_a'] = params['A'] * math.log(params['pod_mass']) + params['B']
        unknowns['t'] = math.sqrt(2 * params['displacement'] / unknowns['pod_a'])
        unknowns['pod_V'] = unknowns['pod_a'] * unknowns['t']


if __name__ == "__main__":
    from openmdao.core.problem import Problem
    from openmdao.core.group import Group
    from openmdao.units.units import convert_units as cu

    p = Problem(root=Group())
    p.root.add('pusher', Pusher())

    p.setup(check=False)

    print ''
    p['pusher.pod_mass'] = float(raw_input('Pod mass with payload (kg): '))

    p.run()

    print ''
    print 'Acceleration:', p['pusher.pod_a'], 'm/s**2'
    print 'Time:', p['pusher.t'], 's'
    print 'Final velocity:', p['pusher.pod_V'], 'm/s'
    print '               ', cu(p['pusher.pod_V'], 'm/s', 'ft/s'), 'ft/s'
    print '               ', cu(p['pusher.pod_V'], 'm/s', 'mi/h'), 'mi/h'

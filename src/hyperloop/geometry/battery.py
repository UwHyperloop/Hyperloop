from openmdao.core.component import Component

class Battery(Component):
    def __init__(self):
        super(Battery, self).__init__()
        self.add_param('time_mission', 2100.0, desc='travel time', units='s')
        self.add_param('cross_section', 1.3, desc='available cross section area for battery pack', units='m**2')
        self.add_param('energy', 0.0, desc='total energy storage requirements', units='kW*h')
        # from http://en.wikipedia.org/wiki/Lithium-ion_battery
        # e ranges from 0.100 to 0.265 kW*h/kg; U ranges from 250 to 739 kW*h/m**3
        self.add_param('e', 0.182, desc='specific energy of Li-ion battery', units='kW*h/kg')
        self.add_param('U', 494.0, desc='energy density of Li-ion battery', units='kW*h/m**3')

        self.add_output('mass', 0.0, desc='total mass of batteries', units='kg')
        self.add_output('volume', 0.0, desc='total volume of batteries', units='m**3')
        self.add_output('len', 0.0, desc='required length of battery pack', units='m')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['mass'] = params['energy'] / params['e']
        unknowns['volume'] = params['energy'] / params['U']
        unknowns['len'] = unknowns['volume'] / params['cross_section']

if __name__ == '__main__':
    from openmdao.core.problem import Problem
    from openmdao.core.group import Group

    p = Problem(root=Group())
    p.root.add('comp', Battery())
    p.setup()
    p.run()

    print 'mass (Kg): %f' % p.root.comp.unknowns['mass']
    print 'energy (kW*hr): %f' % p.root.comp.unknowns['energy']
    print 'volume (m**3): %f' % p.root.comp.unknowns['volume']
    print 'length (m): %f' % p.root.comp.unknowns['len']

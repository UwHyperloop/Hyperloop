from math import pi
import numpy as np
from scipy.optimize import newton

import pylab

from openmdao.core.component import Component
from openmdao.core.group import Group
from openmdao.units.units import convert_units as cu
from openmdao.components.param_comp import ParamComp
from openmdao.components.constraint import ConstraintComp

from pycycle.set_total import SetTotal
from pycycle.thermo_static import SetStaticMN, SetStaticPs
from pycycle.constants import AIR_MIX

from geometry.tube_structure import TubeStructural
from geometry.inlet import InletGeom

class AreaRatio(Component):
    def __init__(self):
        super(AreaRatio, self).__init__()
        self.add_param('tube_r', 0.0, desc='inner radius of tube', units='m')
        self.add_param('inlet_area', 0.0, desc='area of inlet at largest point', units='m**2')
        self.add_param('gamma', 1.41, desc='ratio of specific heats')
        self.add_param('Mach', 1.0, desc='travel Mach')
        self.add_param('Mach_bypass', 0.95, desc='Mach of air passing around pod') # used in TubeLimitFlow

        self.add_output('bypass_area', 0.0, desc='area between tube wall and pod', units='m**2')
        self.add_output('AR', 0.0, desc='ratio between tube area and bypass area')

        self.add_state('AR_resid', 0.0, desc='AR - target AR')

    def solve_nonlinear(self, params, unknowns, resids):
        self.apply_nonlinear(params, unknowns, resids)

    def apply_nonlinear(self, params, unknowns, resids):
        tube_area = pi * (params['tube_r'] ** 2)
        unknowns['bypass_area'] = tube_area - params['inlet_area']
        AR_target = tube_area / unknowns['bypass_area']
        g = params['gamma']
        g_exp = (g + 1.0) / (2.0 * (g - 1.0))
        unknowns['AR'] = ((g + 1.0) / 2.0) ** (-1.0 * g_exp) * ((1.0 + (g - 1.0) / 2.0 * params['Mach'] ** 2) ** g_exp) / params['Mach']
        resids['AR_resid'] = unknowns['AR'] - AR_target

class TubeThermo(Component):
    def __init__(self):
        super(TubeThermo, self).__init__()
        self.add_param('Ps', 99.0, desc='static pressure in tube', units='Pa')
        self.add_param('Ts', 292.1, desc='static temperature in tube', units='degK')
        self.add_param('Mach', 1.0, desc='travel Mach')
        self.add_param('gamma', 1.41, desc='ratio of specific heats')

        self.add_output('Pt', 0.0, desc='total pressure in tube', units='Pa')
        self.add_output('Tt', 0.0, desc='total temperature in tube', units='degK')

    def solve_nonlinear(self, params, unknowns, resids):
        multiplier = (1.0 + (params['gamma'] - 1.0) / 2.0 * params['Mach'] ** 2)
        unknowns['Pt'] = params['Ps'] * multiplier ** (params['gamma'] / (params['gamma'] - 1.0))
        unknowns['Tt'] = params['Ts'] * multiplier

class TubeAero(Component):
    def __init__(self):
        super(TubeAero, self).__init__()
        self.add_param('velocity_tube', 0.0, desc='travel speed where choking occurs', units='m/s')
        self.add_param('velocity_bypass', 0.0, desc='bypass speed where choking occurs', units='m/s')
        self.add_param('bypass_area', 0.0, desc='bypass area', units='m**2')
        self.add_param('tube_r', 0.0, desc='inner radius of tube', units='m')
        self.add_param('rho_tube', 0.0, desc='density in tube', units='g/cm**3')
        self.add_param('rho_bypass', 0.0, desc='density in bypass', units='g/cm**3')

        self.add_output('tube_area', 0.0, desc='cross sectional area of tube', units='m**2')
        self.add_output('W_tube', 0.0, desc='tube demand flow', units='kg/s')
        self.add_output('W_kant', 0.0, desc='Kantrowitz limit flow', units='kg/s')
        self.add_output('W_excess', 0.0, desc='Excess mass flow above the Kantrowitz limit', units='kg/s')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['tube_area'] = pi * params['tube_r'] ** 2
        unknowns['W_tube'] = params['rho_tube'] * params['velocity_tube'] * pi * params['tube_r'] ** 2
        unknowns['W_kant'] = params['rho_bypass'] * params['velocity_bypass'] * params['bypass_area']
        unknowns['W_excess'] = unknowns['W_tube'] - unknowns['W_kant']

class TubeLimitFlow(Group):
    '''Finds the limit velocity for a body traveling through a tube'''
    def __init__(self):
        super(TubeLimitFlow, self).__init__()
        self.add('Mach_param', ParamComp('Mach', 1.0), promotes=['Mach'])
        self.add('Mach_con', ConstraintComp('Mach > 0.0'))
        self.add('tube_struct', TubeStructural())
        self.add('inlet', InletGeom())
        self.add('AR_comp', AreaRatio(), promotes=['Mach_bypass', 'bypass_area'])
        self.add('tube_thermo', TubeThermo())
        self.add('tube_total', SetTotal(init_reacts=AIR_MIX, mode='T'))
        self.add('tube_static', SetStaticMN(init_reacts=AIR_MIX))
        self.add('bypass_total', SetTotal(init_reacts=AIR_MIX, mode='T'))
        self.add('bypass_static', SetStaticMN(init_reacts=AIR_MIX))
        self.add('tube_aero', TubeAero(), promotes=['tube_area', 'W_tube', 'W_kant', 'W_excess'])

        self.connect('Mach', 'AR_comp.Mach')
        self.connect('Mach', 'Mach_con.Mach')
        self.connect('tube_struct.r_inner', 'AR_comp.tube_r')
        self.connect('inlet.area_in', 'AR_comp.inlet_area')
        self.connect('Mach', 'tube_thermo.Mach')
        self.connect('tube_thermo.Pt', 'tube_total.P')
        self.connect('tube_thermo.Tt', 'tube_total.T')
        self.connect('tube_thermo.Pt', 'bypass_total.P')
        self.connect('tube_thermo.Tt', 'bypass_total.T')
        self.connect('tube_total.S', 'tube_static.S')
        self.connect('tube_total.h', 'tube_static.ht')
        self.connect('tube_total.P', 'tube_static.Pt')
        self.connect('tube_total.gamma', 'tube_static.gamt')
        self.connect('bypass_total.S', 'bypass_static.S')
        self.connect('bypass_total.h', 'bypass_static.ht')
        self.connect('bypass_total.P', 'bypass_static.Pt')
        self.connect('bypass_total.gamma', 'bypass_static.gamt')
        self.connect('Mach', 'tube_static.MN_target')
        self.connect('Mach_bypass', 'bypass_static.MN_target')
        self.connect('tube_static.V', 'tube_aero.velocity_tube')
        self.connect('bypass_static.V', 'tube_aero.velocity_bypass')
        self.connect('tube_total.rho', 'tube_aero.rho_tube')
        self.connect('bypass_total.rho', 'tube_aero.rho_bypass')
        self.connect('bypass_area', 'tube_aero.bypass_area')
        self.connect('tube_struct.r_inner', 'tube_aero.tube_r')

def plot_data(p, c='b'):
    '''utility function to make the Kantrowitz Limit Plot''' 
    Machs = []
    W_tube = []
    W_kant = []
    for Mach in np.arange(.2, 1.1, .1):
        p['comp.Mach'] = Mach
        p.run()
        Machs.append(Mach)
        W_kant.append(p['comp.W_kant'])
        W_tube.append(p['comp.W_tube'])
    print 'Area in:', p['comp.inlet.area_in']
    fig = pylab.plot(Machs, W_tube, '-', label="%3.1f Req." % (p['comp.tube_area'] / p['comp.inlet.area_in']), lw=3, c=c)
    pylab.plot(Machs, W_kant, '--', label="%3.1f Limit" % (p['comp.tube_area'] / p['comp.inlet.area_in']), lw=3, c=c)
    pylab.tick_params(axis='both', which='major', labelsize=15)
    pylab.xlabel('Pod Mach Number', fontsize=18)
    pylab.ylabel('Flow Rate (kg/sec)', fontsize=18)
    pylab.title('Tube Flow Limits for Three Area Ratios', fontsize=20)
    return fig

if __name__ == '__main__':
    from openmdao.core.problem import Problem
    from openmdao.core.group import Group
    p = Problem(root=Group())
    comp = p.root.add('comp', TubeLimitFlow())
    p.setup()

    p['comp.tube_struct.r_inner'] = 100.0
    plot_data(p, c='b')

    p['comp.tube_struct.r_inner'] = 150.0
    plot_data(p, c='g')

    p['comp.tube_struct.r_inner'] = 200.0
    plot_data(p, c='r')

    pylab.legend(loc='best')
    
    pylab.gcf().set_size_inches(11, 5.5)
    pylab.gcf().savefig('test2png.png', dpi=130)
    pylab.show()

from math import pi
import numpy as np
from scipy.optimize import newton

import pylab

from openmdao.core.component import Component
from openmdao.units.units import convert_units as cu

from pycycle import flowstation

class TubeLimitFlow(Component):
    '''Finds the limit velocity for a body traveling through a tube'''
    def __init__(self):
        super(TubeLimitFlow, self).__init__()
        self.add_param('r_tube', 1.115, desc='required radius for tube', units='m')
        self.add_param('r_inlet', 0.737, desc='radius of inlet at its largest point', units='m')
        self.add_param('Ps', 99.0, desc='static pressure in tube', units='Pa')
        self.add_param('Ts', 292.1, desc='static temperature in tube', units='degK')
        self.add_param('Mach_pod', 1.0, desc='travel Mach of pod')
        self.add_param('Mach_bypass', 0.95, desc='Mach of air passing around pod')
        
        self.add_output('limit_velocity', 0.0, desc='pod travel speed where flow choking occurs', units='m/s')
        self.add_output('limit_Mach', 0.0, desc='pod travel mach where flow choking occurs')
        self.add_output('W_excess', 0.0, desc='excess tube mass flow above Kantrowitz limit', units='kg/s')
        self.add_output('W_tube', 0.0, desc='tube demand flow', units='kg/s')
        self.add_output('W_kant', 0.0, desc='Kantrowitz limit flow', units='kg/s')
        self.add_output('tube_area', 0.0, desc='area of tube', units='ft**2')
        self.add_output('inlet_area', 0.0, desc='area of inlet', units='ft**2')

    def solve_nonlinear(self, params, unknowns, resids):
        r_tube = cu(params['r_tube'], 'm', 'ft')
        r_inlet = cu(params['r_inlet'], 'm', 'ft')
        unknowns['tube_area'] = pi * (r_tube ** 2)
        unknowns['inlet_area'] = pi * (r_inlet ** 2)
        bypass_area = unknowns['tube_area'] - unknowns['inlet_area']
        Ts = cu(params['Ts'], 'degK', 'degR')
        Ps = cu(params['Ps'], 'Pa', 'psi')
        area_ratio_target = unknowns['tube_area'] / bypass_area
        #iterate over pod speed until the area ratio = A_tube / A_bypass
        out = [None] # makes out[0] a pointer
        def f(Mach):
            out[0] = flow_tube = flowstation.solve(Ts=Ts, Ps=Ps, Mach=Mach)
            g = flow_tube.gamt
            g_exp = (g + 1.0) / (2.0 * (g - 1.0))
            AR = ((g + 1.0) / 2.0) ** (-1.0 * g_exp) * ((1.0 + (g - 1.0) / 2.0 * Mach ** 2) ** g_exp) / Mach
            return AR - area_ratio_target
        #Solve for Mach where AR = AR_target
        unknowns['limit_Mach'] = newton(f, 0.3)
        flow_tube = out[0]
        unknowns['limit_velocity'] = cu(flow_tube.Vflow, 'ft/s', 'm/s')
        #excess mass flow calculation
        flow_tube = flowstation.solve(Ts=Ts, Ps=Ps, Mach=params['Mach_pod'])
        unknowns['W_tube'] = cu(flow_tube.rhos * flow_tube.Vflow * unknowns['tube_area'], 'lbm/s', 'kg/s')
        flow_tube = flowstation.solve(Tt=flow_tube.Tt, Pt=flow_tube.Pt, Mach=params['Mach_bypass'])
        unknowns['W_kant'] = cu(flow_tube.rhos * flow_tube.Vflow * bypass_area, 'lbm/s', 'kg/s')
        unknowns['W_excess'] = unknowns['W_tube'] - unknowns['W_kant']

def plot_data(p, comp, c='b'):
    '''utility function to make the Kantrowitz Limit Plot''' 
    Machs = []
    W_tube = []
    W_kant = []
    for Mach in np.arange(.2, 1.1, .1):
        comp.params['Mach_pod'] = Mach
        p.run()
        Machs.append(Mach)
        W_kant.append(comp.unknowns['W_kant'])
        W_tube.append(comp.unknowns['W_tube'])
    fig = pylab.plot(Machs, W_tube, '-', label="%3.1f Req." % (comp.unknowns['tube_area'] / comp.unknowns['inlet_area']), lw=3, c=c)
    pylab.plot(Machs, W_kant, '--', label="%3.1f Limit" % (comp.unknowns['tube_area'] / comp.unknowns['inlet_area']), lw=3, c=c)
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

    comp.params['r_tube'] = 100.0
    plot_data(p, comp, c='b')

    comp.params['r_tube'] = 150.0
    plot_data(p, comp, c='g')

    comp.params['r_tube'] = 200.0
    plot_data(p, comp, c='r')

    pylab.legend(loc='best')
    
    pylab.gcf().set_size_inches(11, 5.5)
    pylab.gcf().savefig('test2png.png', dpi=130)
    pylab.show()

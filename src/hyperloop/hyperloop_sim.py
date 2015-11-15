from openmdao.core.group import Group
from openmdao.core.component import Component
from openmdao.components.indep_var_comp import IndepVarComp
from openmdao.components.exec_comp import ExecComp
from openmdao.drivers.scipy_optimizer import ScipyOptimizer
from openmdao.units.units import convert_units as cu

from pycycle.components.flow_start import FlowStart

from cycle.compression_system import CompressionSystem
from cycle.splitter import SplitterW
from geometry.pod import Pod
from aero import Aero

from math import pi, sqrt


class TubeFlow(Component):
    '''Calculates initial flow values approaching the pod.'''

    def __init__(self):
        super(TubeFlow, self).__init__()

        self.add_param('pod_MN', 0.5, desc='travel Mach of the pod')
        self.add_param('gamma', 1.41, desc='ratio of specific heats of air')
        self.add_param('tube_T', 292.6, desc='static temperature of tube', units='degK')
        self.add_param('tube_P', 99.0, desc='static pressure of tube', units='Pa')
        self.add_param('tube_area', 2.0, desc='cross sectional area of tube', units='m**2')
        self.add_param('R', 286.0, desc='specific gas constant for flow', units='m**2/s**2/degK')

        self.add_output('W', 0.0, desc='weight flow entering compression system', units='kg/s')
        self.add_output('Pt', 0.0, desc='total pressure of flow entering compression system',
                units='Pa')
        self.add_output('Tt', 0.0, desc='total temperature of flow entering compression system',
                units='degK')

        self.add_param('converted_bypass_area', 0.0, units='m**2')
        self.add_param('converted_inlet_area', 0.0, units='m**2')

    def solve_nonlinear(self, params, unknowns, resids):
        gam = params['gamma']
        MN = params['pod_MN']
        Ts = params['tube_T']
        Ps = params['tube_P']
        R = params['R']
        multiplier = (1.0 + (gam - 1.0) / 2.0 * MN ** 2)
        unknowns['Pt'] = Ps * multiplier ** (gam / (gam - 1.0))
        unknowns['Tt'] = Ts * multiplier
        unknowns['W'] = Ps / R / Ts * params['tube_area'] * MN * sqrt(gam * R * Ts)


class BypassFlow(Component):
    '''Calculates flow values around the pod (in the external bypass).'''

    def __init__(self):
        super(BypassFlow, self).__init__()

        self.add_param('rhot', 0.0, desc='total density of air in bypass', units='kg/m**3')
        self.add_param('Tt', 99.0, desc='total temperature of air in bypass', units='degK')
        self.add_param('bypass_MN', 0.9, desc='Mach number of air in bypass')
        self.add_param('gamma', 1.41, desc='specific heat ratio for air')
        self.add_param('R', 286.0, desc='specific gas constant for flow', units='m**2/s**2/degK')
        self.add_param('bypass_area', 0.5, desc='available bypass area', units='m**2')
        self.add_param('total_W', 0.0, desc='total mass flow through bypass and compression system',
                units='kg/s')
        self.add_param('percent_into_bypass', 1.0 - 1e-4,
                desc='proportion of tube flow to force through the bypass until choked')

        self.add_output('bypass_W', 0.0, desc='mass flow through bypass', units='kg/s')

    def solve_nonlinear(self, params, unknowns, resids):
        gam = params['gamma']
        MN = params['bypass_MN']
        multiplier = (1.0 + (gam - 1.0) / 2.0 * MN ** 2)
        rhos = params['rhot'] * multiplier ** (1.0 / (1.0 - gam))
        Ts = params['Tt'] / multiplier
        Vflow = MN * sqrt(gam * params['R'] * Ts)
        unknowns['bypass_W'] = min(rhos * Vflow * params['bypass_area'],
                params['total_W'] * params['percent_into_bypass'])


class HyperloopSim(Group):
    def __init__(self):
        super(HyperloopSim, self).__init__()

        pod_promotes = ('cross_section', 'bypass_area', 'tube_P', 'tube_T', 'tube_r', 'tube_area',
                'fill_area')
        tube_fl_promotes = ('pod_MN', 'tube_T', 'tube_P', 'tube_area', 'converted_bypass_area',
                'converted_inlet_area')
        bypass_fl_promotes = ('bypass_W', 'bypass_area', 'bypass_MN', 'percent_into_bypass')

        self.add('pod', Pod(), promotes=pod_promotes)
        self.add('tube_flow', TubeFlow(), promotes=tube_fl_promotes)
        self.add('start', FlowStart())
        self.add('bypass_flow', BypassFlow(), promotes=bypass_fl_promotes)
        self.add('split', SplitterW(mode='area'))
        self.add('compression_system', CompressionSystem())

        # non-essential boundary params here to provide definite default values
        self.add('tube_P_param', IndepVarComp('tube_P', 99.0, units='Pa'), promotes=['*'])
        self.add('pod_MN_param', IndepVarComp('pod_MN', 0.3), promotes=['*'])
        self.add('internal_bypass_MN_param', IndepVarComp('internal_bypass_MN', 0.9),
                promotes=['*'])
        self.add('comp2_mouth_MN_param', IndepVarComp('comp2_mouth_MN', 0.5), promotes=['*'])
        self.add('air_bearing_W_param', IndepVarComp('air_bearing_W', 0.2, units='kg/s'),
                promotes=['*'])
        self.add('comp1_exit_MN_param', IndepVarComp('comp1_exit_MN', 0.8), promotes=['*'])
        self.add('comp2_exit_MN_param', IndepVarComp('comp2_exit_MN', 0.8), promotes=['*'])
        self.add('inlet_area_param', IndepVarComp('inlet_area', 0.785, units='m**2'),
                promotes=['*'])

        # connect known flow values to FlowStart component
        self.connect('tube_flow.W', 'start.W')
        self.connect('tube_flow.Pt', 'start.P')
        self.connect('tube_flow.Tt', 'start.T')
        self.connect('pod_MN', 'start.MN_target')

        # connect approaching flow values to external bypass flow values
        self.connect('start.Fl_O:tot:rho', 'bypass_flow.rhot')
        self.connect('start.Fl_O:tot:T', 'bypass_flow.Tt')
        self.connect('start.Fl_O:stat:W', 'bypass_flow.total_W')

        self.connect('bypass_W', 'split.W1')
        self.connect('bypass_area', 'split.area_out1_target')
        self.connect('inlet_area', 'split.area_out2_target')
        self.connect('comp2_mouth_MN', 'compression_system.split.MN_out1_target')
        self.connect('comp1_exit_MN', 'compression_system.comp1_funnel.MN_out_target')
        self.connect('comp2_exit_MN', 'compression_system.comp2_funnel.MN_out_target')
        self.connect('internal_bypass_MN', 'compression_system.split.MN_out2_target')
        self.connect('air_bearing_W', 'compression_system.split.W1')

        self.connect('tube_P', 'compression_system.nozzle.Ps_exhaust')

        self.connect('split.Fl_O2:stat:MN', 'compression_system.diffuser.MN_out_target')
                # no diffuser

        self.add('comp1_cfm_calc', ExecComp('comp1_cfm = W * 60.0 / rho'), promotes=['comp1_cfm'])
        self.connect('split.Fl_O2:stat:W', 'comp1_cfm_calc.W')
        self.connect('split.Fl_O2:stat:rho', 'comp1_cfm_calc.rho')
        
        CompressionSystem.connect_flow(self, 'start.Fl_O', 'split.Fl_I',
                connect_FAR=False)
        CompressionSystem.connect_flow(self, 'split.Fl_O2',
                'compression_system.inlet.Fl_I', connect_stat=False,
                connect_FAR=False)

    @staticmethod
    def p_factory(tube_P=99.0, tube_T=292.6, pod_MN=0.2, inlet_area=0.33,
        cross_section=0.82, tube_r=0.9, fill_area=0.214, bypass_MN=0.9):
        '''
        Sets up an OpenMDAO system for a basic scenario and returns the top-
        level problem.

        Parameters
        ----------
        tube_P : float
            Tube pressure in Pa.
        tube_T : float
            Tube temperature in degK.
        pod_MN : float
            Travel Mach of the pod.
        inlet_area : float
            Area of inlet to compression system in m**2.
        cross_section : float
            Cross-sectional area of the pod at its widest point in m**2.
        tube_r : float
            Inner radius of tube in m.
        fill_area : float
            Cross-sectional area of tube filled by concrete floor in m**2.
        bypass_MN : float
            Desired maximum Mach number of air passing around pod.

        Returns
        -------
        openmdao.core.problem.Problem
        '''

        from openmdao.core.problem import Problem

        g = HyperloopSim()
        p = Problem(root=g)
        
        p.setup(check=False)

        # tube flow
        p['tube_P'] = tube_P
        p['tube_T'] = tube_T
        p['pod_MN'] = pod_MN
        p['tube_r'] = tube_r
        p['fill_area'] = fill_area
        p['bypass_MN'] = bypass_MN

        # compression system
        p['inlet_area'] = inlet_area
        p['cross_section'] = cross_section
        p['compression_system.comp1.PR_design'] = 1.5
        p['compression_system.comp1.eff_design'] = 0.8
        p['comp1_exit_MN'] = 0.35 # keep internal MN greater than or equal to MN of bypass to avoid
                # trailing vacuum
        p['air_bearing_W'] = 1e-4 # negligible
        p['internal_bypass_MN'] = 0.9
        p['comp2_mouth_MN'] = 0.8
        p['compression_system.nozzle.dPqP'] = 0.0
        p['compression_system.comp2.PR_design'] = 1.0
        p['compression_system.comp2.eff_design'] = 1.0
        p['comp2_exit_MN'] = 0.8

        return p


if __name__ == "__main__":
    print 'Setting up...'

    # Set up OpenMDAO problem `p` with pod traveling at Mach 0.2
    p = HyperloopSim.p_factory(pod_MN=0.35, inlet_area=0.4735, cross_section=0.8538)
    # Change some problem parameters

    p['percent_into_bypass'] = 1.0 - p['inlet_area'] / p['tube_area']

    print 'Setup complete.'
    print ''
    print 'Ambient tube pressure:', ' ' * 17, p['tube_P'], 'Pa'
    print 'Ambient tube temperature:', ' ' * 14, cu(p['tube_T'], 'degK', 'degC'), 'degC'
    print 'Maximum pod cross-section:', ' ' * 13, p['cross_section'], 'm**2'
    print 'Inlet area:', ' ' * 28, p['inlet_area'], 'm**2'
    print 'Comp 1 PR:', ' ' * 29, p['compression_system.comp1.PR_design']
    print 'Comp 2 PR:', ' ' * 29, p['compression_system.comp2.PR_design']
    print ''
    print 'WARNING: Temp values for compressors are wrong. With temp, proceed with caution.'
    print 'Optimizing...'

    # Run problem
    p.run()

    # Display some results
    print ''
    print 'Velocity:', cu(p['start.Fl_O:stat:V'], 'ft/s', 'm/s'), 'm/s'
    print 'Tube area:', ' ' * 29, p['tube_area'], 'm**2'
    print 'Ambient tube pressure:', ' ' * 17, p['tube_P'], 'Pa'
    print 'Ambient tube temperature:', ' ' * 14, cu(p['tube_T'], 'degK', 'degC'), 'degC'
    print 'Maximum pod cross-section:', ' ' * 13, p['cross_section'], 'm**2'
    print ''
    print 'Total mass flow:', ' ' * 23, p['split.split_calc.W_in'], 'kg/s'
    print 'Mass flow through bypass:', ' ' * 14, p['bypass_W'], 'kg/s'
    print 'Mass flow through compression system:', ' ' * 2, p['split.split_calc.W2'], 'kg/s'
    print ''
    print 'Comp 1 PR:', ' ' * 23, p['compression_system.comp1.PR_design']
    print 'Comp 1 pwr req:', ' ' * 24, -cu(p['compression_system.comp1.power'], 'hp', 'W'), 'W'
    print 'Total pressure comp 1 exit:', ' ' * 3, cu(p['compression_system.comp1.Fl_O:tot:P'],
            'psi', 'Pa'), 'Pa'
    print 'Area comp 1 exit:', ' ' * 22, p['compression_system.comp1_funnel.Fl_O:stat:area'], 'in**2'
    print 'CFM into comp 1:', ' ' * 23, p['comp1_cfm'], 'ft**3/min'
    print ''
    print 'Area internal bypass:', ' ' * 18, cu(p['compression_system.split.Fl_O2:stat:area'],
            'inch**2', 'm**2'), 'm**2'
    print ''
    print 'Real flow h', p['compression_system.comp1.real_flow.h']
    print 'Real flow P', p['compression_system.comp1.real_flow.P']
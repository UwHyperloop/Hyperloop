from openmdao.core.component import Component
from openmdao.core.group import Group
from openmdao.units.units import convert_units as cu
from openmdao.components.exec_comp import ExecComp
from openmdao.components.param_comp import ParamComp

from pycycle import species_data
from pycycle.constants import AIR_MIX
from pycycle.components.flow_start import FlowStart
from pycycle.components.inlet import Inlet
from pycycle.components.compressor import Compressor
from pycycle.components.nozzle import Nozzle
from pycycle.flowstation import FlowIn

from splitter import SplitterW
from transmogrifier import Transmogrifier

class Performance(Component):
    def __init__(self):
        super(Performance, self).__init__()
        self.add_param('C1_pwr', 0.0, units='hp')
        self.add_param('C2_pwr', 0.0, units='hp')
        self.add_param('Fg', 0.0, desc='gross thrust', units='lbf')
        self.add_param('F_ram', 0.0, units='lbf')
        self.add_param('Ps_bearing_target', 0.0, units='psi')
        self.add_param('Ps_bearing', 0.0, units='psi')
        
        self.add_output('pwr', 0.0, desc='total power required', units='hp')
        self.add_output('Fnet', 0.0, desc='net force', units='lbf')
        self.add_state('Ps_bearing_resid', 0.0, units='psi')

    def solve_nonlinear(self, params, unknowns, resids):
        self.apply_nonlinear(params, unknowns, resids)

    def apply_nonlinear(self, params, unknowns, resids):
        unknowns['pwr'] = params['C1_pwr'] + params['C2_pwr']
        unknowns['Fnet'] = params['Fg'] + params['Fram']
        resids['Ps_bearing_resid'] = params['Ps_bearing'] - params['Ps_bearing_target']

class CompressionSystem(Group):

    @staticmethod
    def connect_flow(group, Fl_O_name, Fl_I_name, connect_stat=True, connect_FAR=True):
        for v_name in ('h', 'T', 'P', 'rho', 'gamma', 'Cp', 'Cv', 'S'):
            for prefix in (('tot', 'stat') if connect_stat else ('tot',)):
                group.connect('%s:%s:%s' % (Fl_O_name, prefix, v_name), '%s:%s:%s' % (Fl_I_name, prefix, v_name))
        if connect_stat:
            for stat in ('V', 'MN', 'area', 'W'):
                group.connect('%s:stat:%s' % (Fl_O_name, stat), '%s:stat:%s' % (Fl_I_name, stat))
        else:
            group.connect('%s:stat:W' % Fl_O_name, '%s:stat:W' % Fl_I_name)
        if connect_FAR:
            group.connect('%s:FAR' % Fl_O_name, '%s:FAR' % Fl_I_name)

    def __init__(self):
        super(CompressionSystem, self).__init__()

        self.thermo_data = species_data.janaf
        self.elements = AIR_MIX

        gas_thermo = species_data.Thermo(self.thermo_data, init_reacts=self.elements)
        self.gas_prods = gas_thermo.products
        self.num_prod = len(self.gas_prods)

        self.add('MN_param', ParamComp('pod_MN', 0.5), promotes=['*'])
        self.add('gam_param', ParamComp('gamma', 1.41), promotes=['*'])
        self.add('T_param', ParamComp('tube_T', 292.6, units='degK'), promotes=['*'])
        self.add('P_param', ParamComp('tube_P', 99.0, units='Pa'), promotes=('*',))
        self.add('area_param', ParamComp('inlet_area', 2.0, units='m**2'), promotes=['*'])

        self.add('W_start', ExecComp('W = P / R / T * area * MN * math.sqrt(gamma * R * T)'), promotes=['R'])
        self.add('Pt_start', ExecComp('Pt = Ps * (1.0 + (gamma - 1.0) / 2.0 * MN ** 2) ** (gamma / (gamma - 1.0))'))
        self.add('Tt_start', ExecComp('Tt = Ts * (1.0 + (gamma - 1.0) / 2.0 * MN ** 2)'))

        self.add('start', FlowStart())
        self.add('inlet', Inlet())
        self.add('diffuser', Transmogrifier())
        self.add('comp1', Compressor())
        self.add('comp1_funnel', Transmogrifier(mode='MN')) # calculates statics based on exit Mach
        self.add('split', SplitterW(mode='MN'))
        self.add('nozzle', Nozzle(elements=AIR_MIX))
        self.add('comp2', Compressor())
        self.add('comp2_funnel', Transmogrifier()) # calculates statics based on exit area
        self.add('perf', Performance())

        self.connect('pod_MN', 'W_start.MN')
        self.connect('gamma', 'W_start.gamma')
        self.connect('tube_T', 'W_start.T')
        self.connect('tube_P', 'W_start.P')
        self.connect('inlet_area', 'W_start.area')
        self.connect('W_start.W', 'start.W')

        self.connect('pod_MN', 'Pt_start.MN')
        self.connect('gamma', 'Pt_start.gamma')
        self.connect('tube_P', 'Pt_start.Ps')
        self.connect('Pt_start.Pt', 'start.P')

        self.connect('pod_MN', 'Tt_start.MN')
        self.connect('gamma', 'Tt_start.gamma')
        self.connect('tube_T', 'Tt_start.Ts')
        self.connect('Tt_start.Tt', 'start.T')

        self.connect('pod_MN', 'start.MN_target')

        conn_fl = CompressionSystem.connect_flow
        conn_fl(self, 'start.Fl_O', 'inlet.Fl_I', connect_FAR=False)
        conn_fl(self, 'inlet.Fl_O', 'diffuser.Fl_I', connect_stat=False)
        conn_fl(self, 'diffuser.Fl_O', 'comp1.Fl_I')
        conn_fl(self, 'comp1.Fl_O', 'comp1_funnel.Fl_I', connect_stat=False)
        conn_fl(self, 'comp1_funnel.Fl_O', 'split.Fl_I')
        conn_fl(self, 'split.Fl_O1', 'comp2.Fl_I')
        conn_fl(self, 'comp2.Fl_O', 'comp2_funnel.Fl_I', connect_stat=False)
        conn_fl(self, 'split.Fl_O2', 'nozzle.Fl_I', connect_FAR=False)

        self.connect('comp1.power', 'perf.C1_pwr')
        self.connect('comp2.power', 'perf.C2_pwr')
        self.connect('comp2_funnel.Fl_O:stat:P', 'perf.Ps_bearing')
        self.connect('nozzle.Fg', 'perf.Fg')
        self.connect('inlet.F_ram', 'perf.F_ram')


if __name__ == "__main__":
    from openmdao.core.problem import Problem

    g = CompressionSystem()
    p = Problem(root=g)

    p.setup()

    p['tube_P'] = 99.0
    p['tube_T'] = 292.6
    p['gamma'] = 1.41 # constant for air
    p['pod_MN'] = 0.5
    p['R'] = 286.0 # constant for air
    p['inlet_area'] = 2.0

    p['inlet.ram_recovery'] = 1.0

    p['diffuser.area_out_target'] = p['inlet_area'] # no diffuser if equal to inlet_area

    p['comp1.PR_design'] = 12.47
    p['comp1.eff_design'] = 0.8

    p['comp1_funnel.MN_out_target'] = 0.6

    p['split.W1'] = 0.44 # weight flow requirement of the bearing system goes here
    p['split.MN_out1_target'] = 0.6
    p['split.MN_out2_target'] = 1.0

    p['nozzle.dPqP'] = 0.0
    p['nozzle.Ps_ideal'] = 99.0 # TODO is this right??

    p['comp2.PR_design'] = 5.0
    p['comp2.eff_design'] = 0.8

    p['comp2_funnel.MN_out_target'] = 0.6

    p.run()

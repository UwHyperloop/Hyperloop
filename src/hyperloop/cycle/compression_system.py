from openmdao.core.component import Component
from opnemdao.core.group import Group
from openmdao.units.units import convert_units as cu

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
        self.add_param('inlet_area', 0.0, desc='cross section of flow entering inlet', units='m**2')
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
    # Required data: Fl_I:stat:V, inlet_area, Fl_I:stat:rho, 

    @staticmethod
    def connect_flow(group, Fl_O_name, Fl_I_name):
        for v_name in ('h', 'T', 'P', 'rho', 'gamma', 'Cp', 'Cv', 'S', 'n', 'n_moles'):
            for prefix in ('tot', 'stat'):
                group.connect('%s:%s:%s' % (Fl_O_name, prefix, v_name), '%s:%s:%s' % (Fl_I_name, prefix v_name))
        for stat in ('V', 'Vsonic', 'MN', 'area', 'W', 'Wc'):
            group.connect('%s:stat:%s' % (Fl_O_name, stat), '%s:stat:%s' % (Fl_I_name, stat))
        group.connect('%s:FAR' % Fl_O_name, '%s:FAR' % Fl_I_name)

    def __init__(self):
        super(CompressionSystem, self).__init__()

        flow_in = FlowIn('Fl_I', self.num_prod)
        self.add('flow_in', flow_in, promotes=flow_in.flow_in_vars)
        self.add('W_start', ExecComp('W = V * rho * area'))
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

        self.connect('Fl_I:stat:V', 'W_start.V')
        self.connect('Fl_I:stat:rho', 'W_start.rho')
        self.connect('inlet_area', 'W_start.area')
        self.connect('W_start.W', 'start.W')

        conn_fl = CompressionSystem.connect_flow
        conn_fl(self, 'start', 'inlet.Fl_I')
        conn_fl(self, 'inlet.Fl_O', 'diffuser.Fl_I')
        conn_fl(self, 'diffuser.Fl_O', 'comp1.Fl_I')
        conn_fl(self, 'comp1.Fl_O', 'comp1_funnel.Fl_I')
        conn_fl(self, 'comp1_funnel.Fl_O', 'split.Fl_I')
        conn_fl(self, 'split.Fl_O1', 'comp2.Fl_I')
        conn_fl(self, 'comp2.Fl_O', 'comp2_funnel.Fl_I')
        conn_fl(self, 'split.Fl_O2', 'nozzle.Fl_I')

        self.connect('comp1.power', 'perf.C1_pwr')
        self.connect('comp2.power', 'perf.C2_pwr')
        self.connect('comp2.Fl_O:stat:P', 'perf.Ps_bearing')
        self.connect('nozzle.Fg', 'perf.Fg')
        self.connect('inlet.F_ram', 'perf.F_ram')


if __name__ == "__main__":
    from openmdao.core.problem import Problem

    p = Problem()
    g = p.root = CompressionSystem()

    p.setup()

    g.params['Fl_I:stat:P'] = cu(0.01436, 'lbf', 'Pa')
    g.params['Fl_I:stat:T'] = cu(525.6, 'degR', 'degK')

    g.params['inlet.ram_recovery'] = 1.0

    g.params['diffuser.area_out_target'] = g.params['inlet_area']

    g.params['comp1.PR_design'] = 12.47
    g.params['comp1.eff_design'] = 0.8

    g.params['comp1_funnel.MN_out_target'] = 0.6

    g.params['split.W1'] = 0.44 # weight flow requirement of the bearing system goes here
    g.params['split.MN_out1_target'] = 0.6
    g.params['split.MN_out2_target'] = 1.0

    g.params['nozzle.dPqP'] = 0.0
#    g.params['Ps_exhaust'] = 

    g.params['comp2.PR_design'] = 5.0
    g.params['comp2.eff_design'] = 0.8

    g.params['comp2_funnel.MN_out_target'] = 0.6

    p.run()




# TODO still working..............


        self.connect('duct2.Fl_O.Ps', 'perf.Ps_bearing')

        #Input variable pass_throughs to the assembly boundary
        #Compress -> Tube
        self.connect('W_in', 'tube.W')
        self.connect('Ts_tube','tube.Ts')
        self.connect('Ps_tube', 'tube.Ps')
        self.connect('Mach_pod_max', 'tube.Mach')
        #Compress -> Inlet
        self.connect('Mach_c1_in', 'inlet.MNexit_des')
        #Compress -> C1
        self.connect('c1_PR_des','comp1.PR_des')
        #Compress -> C2
        self.connect('c2_PR_des','comp2.PR_des')
        #Compress -> Splitter
        self.connect('W_bearing_in', 'split.W1_des')
        #Compress -> Perf
        self.connect('Ps_bearing', 'perf.Ps_bearing_target')
        
        #Output variable pass_throughs to the assembly boundary
        self.connect('tube.Fl_O.rhot', 'rho_air') #promoted for aero calc
        self.connect('tube.Fl_O.area', 'area_inlet_in')
        self.connect('tube.Fl_O.Vflow', 'speed_max')
        self.connect('inlet.Fl_O.area', 'area_c1_in')
        self.connect('comp1.Fl_O.area', 'area_c1_out')
        self.connect('nozzle.Fl_O.area', 'nozzle_flow_area')
        self.connect('nozzle.Fl_O', 'nozzle_Fl_O')
        self.connect('duct2.Fl_O', 'bearing_Fl_O')
        self.connect('perf.F_net','F_net')
        self.connect('perf.pwr', 'pwr_req') 
        self.connect('perf.Ps_bearing_residual', 'Ps_bearing_residual')



if __name__ == "__main__": 
    from math import pi
    from openmdao.main.api import set_as_top

    hlc = set_as_top(CompressionSystem())
    hlc.Mach_pod_max = 1
    hlc.run()

    print "pwr: ", hlc.comp1.pwr+hlc.comp2.pwr,hlc.comp1.pwr,hlc.comp2.pwr 
    print "tube area:", hlc.tube.Fl_O.area 
    print "tube Ps", hlc.tube.Fl_O.Ps, hlc.tube.Fl_O.Pt
    print "tube Rhos", hlc.tube.Fl_O.rhos
    print "tube W", hlc.tube.W
    print "inlet W", hlc.inlet.Fl_I.W
    print "tube rad: ", (hlc.tube.Fl_O.area/pi)**.5
    print "tube V: ", hlc.tube.Fl_O.Vflow, hlc.tube.Fl_O.Mach

    fs = hlc.tube.Fl_O





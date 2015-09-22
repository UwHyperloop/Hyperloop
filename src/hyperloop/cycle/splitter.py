from openmdao.core.group import Group
from openmdao.core.component import Component

from scipy.optimize.newton import newton

from pycycle.constants import g_c, AIR_FUEL_MIX, AIR_MIX
from pycycle.set_total import SetTotal
from pycycle.thermo_static import SetStaticMN
from pycycle import species_data
from pycycle.flowstation import FlowIn, PassThrough


class SplitterWCalc(Component):
    """Calculates statics based on weight flow"""

    def __init__(self, mode="MN"):
        super(SplitterWCalc, self).__init__()

        if mode == "MN" or mode == "area":
            pass
        else:
            raise ValueError('Only "MN" and "area" are allowed values for mode')

        self.mode = mode

        self.add_param('W1', 0.44, desc='Weight flow for Fl_O1', units='kg/s')
        if mode == 'MN':
            self.add_param('MN_out1_target', 1.0, desc='Exit Mach for Fl_O1', units='m**2')
            self.add_param('MN_out2_target', 1.0, desc='Exit Mach for Fl_O2', units='m**2')

        elif mode == 'area':
            self.add_param('area_out1_target', 1.0, desc='Exit area for Fl_O1', units='m**2')
            self.add_param('area_out2_target', 1.0, desc='Exit area for Fl_O2', units='m**2')

            self.add_param('exit_area1', 1.0, desc='Proposed exit area for Fl_O1')
            self.add_param('exit_area2', 1.0, desc='Proposed exit area for Fl_O2')

            self.add_state('area_resid1', 0.0, desc='exit_area_1 - area by Mach', units='m**2')
            self.add_state('area_resid2', 0.0, desc='exit_area_2 - area by Mach', units='m**2')

    def solve_nonlinear(self, params, unknowns, resids):
        self.apply_nonlinear(self, params, unknowns, resids)

    def apply_nonlinear(self, params, unknowns, resids):
        if mode == 'area':
            resids['area_resid1'] = params['exit_area1'] - params['area_out1_target']
            resids['area_resid2'] = params['exit_area2'] - params['area_out2_target']

        elif mode == 'MN':
            pass


class SplitterW(Group):
    """A Group that models a splitter with known weight flows"""

    def __init__(self, thermo_data=species_data.janaf, elements=AIR_MIX, mode="MN"):
        super(Splitter, self).__init__()

        if mode == "MN" or mode == "area":
            pass
        else:
            raise ValueError('Only "MN" and "area" are allowed values for mode')

        self.mode = mode

        self.thermo_data = thermo_data
        self.elements = elements

        gas_thermo = species_data.Thermo(thermo_data, init_reacts=elements)
        self.gas_prods = gas_thermo.products
        self.num_prod = len(self.gas_prods)

        # Create inlet flowstation
        flow_in = FlowIn('Fl_I', self.num_prod)
        self.add('flow_in', flow_in, promotes=flow_in.flow_in_vars)

        out1_stat = SetStaticMN(thermo_data, elements, 'Fl_O1:stat')
        self.add('out1_stat', out1_stat, promotes=out1_stat.flow_out_vars)

        self.connect('Fl_I:tot:P', 'out1_stat.Pt')
        self.connect('Fl_I:tot:S', 'out1_stat.S')

        out2_stat = SetStaticMN(thermo_data, elements, 'Fl_O2:stat')
        self.add('out2_stat', out2_stat, promotes=out2_stat.flow_out_vars)

        self.connect('Fl_I:tot:P', 'out2_stat.Pt')
        self.connect('Fl_I:tot:S', 'out2_stat.S')

        # total vars
        for v_name in ('h','T','P','S','rho','gamma','Cp','Cv','n','n_moles'):
            for n in range(2):
                self.add('%s_passthru_%d' % (v_name, n), PassThrough('Fl_I:tot:%s' % v_name, 'Fl_O%d:tot:%s' % v_name, 0.0), promotes=['*'])

        self.add('split_calc', SplitterWCalc(mode), promotes=['*'])

        if mode == 'area':
            self.connect('Fl_O1:stat:area', 'exit_area1')
            self.connect('Fl_O2:stat:area', 'exit_area2')

        elif mode == 'MN':
            self.connect('MN_out1_target', 'out1_stat.MN_target')
            self.connect('MN_out2_target', 'out2_stat.MN_target')

        self.add('FAR_passthru', PassThrough('Fl_I:FAR', 'Fl_O:FAR', 0.0), promotes=['*'])

    def solve_nonlinear(self, params, unknowns, resids):
        self.out1_static.params['W'] = params['W1']
        self.out2_static.params['W'] = params['Fl_I:stat:W'] - params['W1']

        if self.mode == 'area':
            def f(Mach_guess, exit_n):
                flow = self.out1_stat if exit_n == 1 else self.out2_stat
                flow.params['MN_target'] = Mach_guess
                flow.solve_nonlinear(flow.params, flow.unknowns, flow.resids)
                self.split_calc.apply_nonlinear(self.split_calc.params, self.split_calc.unknowns, self.split_calc.resids)
                return resids['area_resid%d' % exit_n]
            
            flow = self.out2_stat
            flow.params['MN_target'] = params['Fl_I:stat:Mach']
            flow.solve_nonlinear(flow.params, flow.unknowns, flow.resids)

            for n in (1, 2):
                newton(f, params['Fl_I:stat:Mach'], args=(n,))

        elif self.mode == 'MN':
            super(SplitterW, self).solve_nonlinear(params, unknowns, resids)


if __name__ == "__main__":
#    from openmdao.core.problem import Problem
#
#    p = Problem()
#    p.root = Splitter()
#
#    p.setup()
#    p.run()
    print '__main__ not implemented'

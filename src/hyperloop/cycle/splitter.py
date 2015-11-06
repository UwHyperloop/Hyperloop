import numpy as np

from openmdao.core.group import Group
from openmdao.core.component import Component
from openmdao.components.indep_var_comp import IndepVarComp

from scipy.optimize import newton

from pycycle.constants import g_c, AIR_FUEL_MIX, AIR_MIX
from pycycle.set_total import SetTotal
from pycycle.thermo_static import SetStaticMN, SetStaticArea
from pycycle import species_data
from pycycle.flowstation import FlowIn, PassThrough


class SplitterWCalc(Component):
    """Calculates statics based on weight flow"""

    def __init__(self, mode):
        super(SplitterWCalc, self).__init__()

        self.mode = mode

        self.add_param('W_in', 1.0, desc='total weight flow in', units='kg/s')
        self.add_param('W1', 0.44, desc='weight flow for Fl_O1', units='kg/s')
        if mode == 'MN':
            self.add_param('MN_out1_target', 0.5, desc='Mach number of Fl_O1')
            self.add_param('MN_out2_target', 0.5, desc='Mach number of Fl_O2')
        elif mode == 'area':
            self.add_param('area_out1_target', 0.5, desc='Area of Fl_O1')
            self.add_param('area_out2_target', 0.5, desc='Area of Fl_O2')

        self.add_output('W2', 0.56, desc='Weight flow for Fl_O1', units='kg/s')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['W2'] = params['W_in'] - params['W1']


class SplitterW(Group):
    """A Group that models a splitter with known weight flows"""

    def __init__(self, thermo_data=species_data.janaf, elements=AIR_MIX, mode='MN'):
        super(SplitterW, self).__init__()

        self.mode = mode

        self.thermo_data = thermo_data
        self.elements = elements

        gas_thermo = species_data.Thermo(thermo_data, init_reacts=elements)
        self.gas_prods = gas_thermo.products
        self.num_prod = len(self.gas_prods)

        # Create inlet flowstation
        flow_in = FlowIn('Fl_I', self.num_prod)
        self.add('flow_in', flow_in, promotes=flow_in.flow_in_vars)

        if mode == 'MN':
            promotes = ('W1', 'MN_out1_target', 'MN_out2_target')
        elif mode == 'area':
            promotes = ('W1', 'area_out1_target', 'area_out2_target')
        self.add('split_calc', SplitterWCalc(mode), promotes=promotes)

        if mode == 'MN':
            out1_stat = SetStaticMN(thermo_data, elements, 'Fl_O1:stat')
            out2_stat = SetStaticMN(thermo_data, elements, 'Fl_O2:stat')
        elif mode == 'area':
            out1_stat = SetStaticArea(thermo_data, elements, 'Fl_O1:stat')
            out2_stat = SetStaticArea(thermo_data, elements, 'Fl_O2:stat')
        self.add('out1_stat', out1_stat, promotes=out1_stat.flow_out_vars)
        self.add('out2_stat', out2_stat, promotes=out2_stat.flow_out_vars)

        if mode == 'MN':
            self.connect('MN_out1_target', 'out1_stat.MN_target')
            self.connect('MN_out2_target', 'out2_stat.MN_target')
        elif mode == 'area':
            self.connect('area_out1_target', 'out1_stat.area_target')
            self.connect('area_out2_target', 'out2_stat.area_target')

        self.connect('Fl_I:stat:W', 'split_calc.W_in')

        self.connect('W1', 'out1_stat.W')
        self.connect('Fl_I:tot:h', 'out1_stat.ht')
        self.connect('Fl_I:tot:S', 'out1_stat.S')
        self.connect('Fl_I:tot:n', 'out1_stat.n_guess')

        self.connect('split_calc.W2', 'out2_stat.W')
        self.connect('Fl_I:tot:h', 'out2_stat.ht')
        self.connect('Fl_I:tot:S', 'out2_stat.S')
        self.connect('Fl_I:tot:n', 'out2_stat.n_guess')

        # total vars
        for v_name in ('h','T','P','S','rho','gamma','Cp','Cv','n_moles'):
            for n in (1, 2):
                self.add('%s_passthru%d' % (v_name, n), PassThrough('Fl_I:tot:%s' % v_name, 'Fl_O%d:tot:%s' % (n, v_name), 0.0), promotes=['*'])
        for n in (1, 2):
            self.add('n_passthru%d' % n, PassThrough('Fl_I:tot:n', 'Fl_O%d:tot:n' % n, np.zeros(self.num_prod)), promotes=['*'])

        for n in (1, 2):
            self.add('FAR_passthru%d' % n, PassThrough('Fl_I:FAR', 'Fl_O%d:FAR' % n, 0.0), promotes=['*'])


if __name__ == "__main__":
#    from openmdao.core.problem import Problem
#
#    p = Problem()
#    p.root = Splitter()
#
#    p.setup()
#    p.run()
    print '__main__ not implemented'

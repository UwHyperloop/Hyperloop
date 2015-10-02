import numpy as np

from openmdao.core.group import Group
from openmdao.core.component import Component

from scipy.optimize import newton

from pycycle.constants import g_c, AIR_FUEL_MIX, AIR_MIX
from pycycle.set_total import SetTotal
from pycycle.thermo_static import SetStaticMN
from pycycle import species_data
from pycycle.flowstation import FlowIn, PassThrough


class SplitterWCalc(Component):
    """Calculates statics based on weight flow"""

    def __init__(self, mode="MN"):
        super(SplitterWCalc, self).__init__()

        if mode == "MN":
            pass
        else:
            raise ValueError('Only "MN" is an allowed value for mode')

        self.mode = mode

        self.add_param('W1', 0.44, desc='Weight flow for Fl_O1', units='kg/s')
        if mode == 'MN':
            self.add_param('MN_out1_target', 1.0, desc='Exit Mach for Fl_O1', units='m**2')
            self.add_param('MN_out2_target', 1.0, desc='Exit Mach for Fl_O2', units='m**2')

    def solve_nonlinear(self, params, unknowns, resids):
        self.apply_nonlinear(self, params, unknowns, resids)

    def apply_nonlinear(self, params, unknowns, resids, metadata):
        if self.mode == 'MN':
            pass


class SplitterW(Group):
    """A Group that models a splitter with known weight flows"""

    def __init__(self, thermo_data=species_data.janaf, elements=AIR_MIX, mode="MN"):
        super(SplitterW, self).__init__()

        if mode == "MN":
            pass
        else:
            raise ValueError('Only "MN" is an allowed value for mode')

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
        for v_name in ('h','T','P','S','rho','gamma','Cp','Cv','n_moles'):
            for n in (1, 2):
                self.add('%s_passthru%d' % (v_name, n), PassThrough('Fl_I:tot:%s' % v_name, 'Fl_O%d:tot:%s' % (n, v_name), 0.0), promotes=['*'])
        for n in (1, 2):
            self.add('n_passthru%d' % n, PassThrough('Fl_I:tot:n', 'Fl_O%d:tot:n' % n, np.zeros(self.num_prod)), promotes=['*'])

        self.add('split_calc', SplitterWCalc(mode), promotes=['*'])

        if mode == 'MN':
            self.connect('MN_out1_target', 'out1_stat.MN_target')
            self.connect('MN_out2_target', 'out2_stat.MN_target')

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

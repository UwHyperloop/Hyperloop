import numpy as np

from openmdao.core.group import Group
from openmdao.core.component import Component
from openmdao.solvers.newton import Newton

# from scipy.optimize import brentq

from pycycle.constants import AIR_MIX
from pycycle.thermo_static import SetStaticMN
from pycycle import species_data
from pycycle.flowstation import FlowIn, PassThrough


class TransmogrifierCalc(Component):

    def __init__(self, mode="area"):
        super(TransmogrifierCalc, self).__init__()

        assert mode == "MN"
        
        self.mode = mode

        if mode == 'MN':
            self.add_param('MN_out_target', 1.0, desc='Mach at exit')

    def solve_nonlinear(self, params, unknowns, resids):
        self.apply_nonlinear(params, unknowns, resids)

    def apply_nonlinear(self, params, unknowns, resids):
        if self.mode == 'MN':
            pass


class Transmogrifier(Group):
    """Calculates statics for a component with changing area (e.g. a diffuser or converging nozzle); can be used alone or in sequence with a compressor, splitter, etc. to model area change across the component"""

    def __init__(self, thermo_data=species_data.janaf, elements=AIR_MIX, mode="MN"):
        super(Transmogrifier, self).__init__()

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

        if mode == 'MN':
            promotes = ['MN_out_target']
        self.add('calc', TransmogrifierCalc(mode), promotes=promotes)

        # Calculate statics based on MN
        set_stat = SetStaticMN(thermo_data, elements, 'Fl_O:stat')
        self.add('set_stat', set_stat, promotes=set_stat.flow_out_vars)

        self.connect('Fl_I:tot:gamma', 'set_stat.gamt')
        self.connect('Fl_I:tot:h', 'set_stat.ht')
        self.connect('Fl_I:tot:P', 'set_stat.Pt')
        self.connect('Fl_I:tot:S', 'set_stat.S')
        self.connect('Fl_I:stat:W', 'set_stat.W')
        self.connect('Fl_I:tot:n', 'set_stat.n_guess')

        if mode == 'MN':
            self.connect('MN_out_target', 'set_stat.MN_target')

        for v_name in ('h', 'T', 'P', 'rho', 'gamma', 'Cp', 'Cv', 'S', 'n_moles'):
            self.add('passthru_tot:%s' % v_name, PassThrough('Fl_I:tot:%s' % v_name, 'Fl_O:tot:%s' % v_name, 0.0), promotes=['*'])
        self.add('passthru_tot:n', PassThrough('Fl_I:tot:n', 'Fl_O:tot:n', np.zeros(self.num_prod)), promotes=['*'])
        self.add('passthru_FAR', PassThrough('Fl_I:FAR', 'Fl_O:FAR', 0.0), promotes=['*'])

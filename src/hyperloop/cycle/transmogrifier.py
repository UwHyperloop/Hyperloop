from openmdao.core.group import Group
from openmdao.core.component import Component

from scipy.optimize.newton import brentq

from pycycle.constants import AIR_MIX
from pycycle.thermo_static import SetStaticMN
from pycycle import species_data
from pycycle.flowstation import FlowIn, PassThrough


class TransmogrifierCalc(Component):

    def __init__(self):
        super(TransmogrifierCalc, self).__init__()

        self.add_param('area_out_target', 0.0, desc='area at exit', units='m**2')
        self.add_param('area_out', 0.0, desc='calculated area at exit', units='m**2')
        self.add_param('shock', False, desc='True if component makes supersonic flow subsonic')

        self.add_state('area_resid', 0.0, desc='area_out - area_out_target', units='m**2')

    def solve_nonlinear(self, params, unknowns, unknowns):
        self.apply_nonlinear(params, unknowns, resids)

    def apply_nonlinear(self, params, unknowns, resids):
        resids['area_resid'] = params['area_out'] - params['area_out_target']


class Transmogrifier(Group):
    """Boringly calculates statics for a component with changing area (e.g. a diffuser or converging nozzle)"""

    def __init__(self, thermo_data=species_data.janaf, elements=AIR_MIX):
        super(Transmogrifier, self).__init__()

        self.thermo_data = thermo_data
        self.elements = elements

        gas_thermo = species_data.Thermo(thermo_data, init_reacts=elements)
        self.gas_prods = gas_thermo.products
        self.num_prod = len(self.gas_prods)

        # Create inlet flowstation
        flow_in = FlowIn('Fl_I', self.num_prod)
        self.add('flow_in', flow_in, promotes=flow_in.flow_in_vars)

        # Calculate statics based on MN
        set_stat = SetStaticMN(thermo_data, elements, 'Fl_O:stat')
        self.add('set_stat', set_stat, promotes=set_stat.flow_out_vars)

        self.connect('Fl_I:tot:P', 'set_stat.Pt')
        self.connect('Fl_I:tot:S', 'set_stat.S')
        self.connect('Fl_I:stat:W', 'set_stat.W')

        self.add('calc', TransmogrifierCalc(), promotes=['area_out_target', 'shock'])

        self.connect('Fl_O:stat:area', 'calc.area_out')

        for v_name in ('h', 'T', 'P', 'rho', 'gamma', 'Cp', 'Cv', 'S', 'n', 'n_moles'):
            self.add('passthru_tot:%s' % v_name, PassThrough('Fl_I:tot:%s', 'Fl_O:tot:%s', 0.0), promotes=['*'])

        self.add('passthru_FAR', PassThrough('Fl_I:FAR', 'Fl_O:FAR', 0.0), promotes=['*'])

    def solve_nonlinear(self, params, unknowns, resids):
        def f(MN_exit):
            params['set_stat.MN_target'] = MN_exit
            super(Transmogrifier, self).solve_nonlinear(params, unknowns, resids)
            return resids['calc.area_resid']

        if params['Fl_I:stat:MN'] < 1.0 or params['shock']:
            brentq(f, 1e-4, 1.0)
        else:
            brentq(f, 1.0, 50.0)

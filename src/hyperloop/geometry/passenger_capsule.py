from openmdao.core.component import Component

class PassengerCapsule(Component):
    '''Place holder component for passenger capsule sizing and structural analysis. Currently, just assume the baseline shape from the original proposal'''
    def __init__(self):
        super(PassengerCapsule, self).__init__()
        self.add_param('n_rows', 14, desc='number of rows of seats')
        self.add_param('row_len', 1.5, desc='length of each row of seats', units='m')

        self.add_output('capsule_len', 0.0, desc='overall length of passenger capsule', units='m')
        self.add_output('cross_section', 0.0, desc='cross sectional area of passenger capsule', units='m**2')

    def solve_nonlinear(self, params, unknowns, resids):
        # TODO replace 10% "fudge factor" with more specific accomodations
        unknowns['capsule_len'] = 1.1 * params['n_rows'] * params['row_len'] # 10% fudge factor
        unknowns['cross_section'] = 1.4 # page 15 of the original proposal

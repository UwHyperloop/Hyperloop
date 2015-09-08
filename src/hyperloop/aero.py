from openmdao.core.component import Component

class Aero(Component):
    '''Placeholder for real aerodynamic calculations of the capsule''' 
    def __init__(self):
        self.add_param('coef_drag', 1.0, desc='capsule drag coefficient')
        self.add_param('area_frontal', 18.0, desc='frontal area of capsule', units='m**2')
        self.add_param('velocity_capsule', 600.0, desc='capsule velocity', units='m/s')
        self.add_param('rho', 0.0, desc='tube air density', units='kg/m**3')
        self.add_param('gross_thrust', 0.0, desc='nozzle gross thrust', units='N')
        
        self.add_output('net_force', 0.0, desc='net force with drag considerations', units='N')
        self.add_output('drag', 0.0, desc='drag force', units='N')

    def solve_nonlinear(self, params, unknowns, resids):
        # drag = Cd * rho * Velocity ** 2 * Area / 2.0
        unknowns['drag'] = params['coef_drag'] * params['rho'] * params['velocity_capsule'] ** 2 * params['area_frontal'] / 2.0
        unknowns['net_force'] = params['gross_thrust'] - unknowns['drag']

'''
tubeModel.py - 
    Determines the steady state temperature of the hyperloop tube.
    Calculates Q released/absorbed by hyperloop tube due to:
    Internal Convection, Tube Conduction, Ambient Natural Convection, Solar Flux In, Radiation Out
    
-original calculations from Jeff Berton, ported and extended by Jeff Chin

Compatible with OpenMDAO v1.0.5
'''

from math import log, pi, sqrt, e

from openmdao.core.group import Group
from openmdao.units.units import convert_units as cu

from pycycle import flowstation
from pycycle.cycle_component import CycleComponent

class TubeWallTemp(CycleComponent):
    '''Calculates Q released/absorbed by the hyperloop tube'''
    def __init__(self):
        super(TubeWallTemp, self).__init__()
        self._add_flowstation('flow_nozzle')
        self._add_flowstation('flow_bearings')

        self.add_param('r_tube_outer', 3.006, desc='outer radius of tube', units='m')
        self.add_param('tube_len', 482803.0, desc='length of one trip', units='m')
        self.add_param('n_pods', 34, desc='number of Pods in the tube at a given time')
        self.add_param('temp_boundary', 322.0, desc='average temperature of tube wall', units='degK')
        self.add_param('temp_ambient', 305.6, desc='average temperature of outside air', units='degK')

        self.add_param('insolation', 1000.0, desc='solar irradiation at sea level on a clear day', units='W/m**2')
        self.add_param('nn_incidence_factor', 0.7, desc='non-normal incidence factor')
        self.add_param('reflectance', 0.5, desc='solar reflectance index')
        self.add_param('emissivity', 0.5, desc='emissivity of the tube', units='W')
        self.add_param('sb_const', 0.00000005670373, desc='Stefan-Boltzmann constant', units='W/m**2/degK**4')
        
        self.add_output('radiating_area', 337486.1, desc='tube radiating area', units='m**2')
        self.add_output('GrDelTL3', 1946216.7, desc='heat radiated to the outside', units='1/ft**3/F')
        self.add_output('Pr', 0.707, desc='Prandtl')
        self.add_output('Gr', 12730351223.0, desc='Grashof #')
        self.add_output('Ra', 8996312085.0, desc='Rayleigh #')
        self.add_output('Nu', 232.4543713, desc='Nusselt #')
        self.add_output('k', 0.02655, desc='thermal conductivity', units='W/m/degK')
        self.add_output('h', 0.845464094, desc='heat radiated to the outside', units='W/m**2/degK')
        self.add_output('convection_area', 3374876.115, desc='convection area', units='W')
        self.add_output('Qradiated_nat_convection_per_area', 286900419.0, desc='heat radiated per area to the outside via natural convection', units='W/m**2/degK')
        self.add_output('Qradiated_nat_convection_tot', 286900419.0, desc='total heat radiated to the outside via natural convection', units='W')
        self.add_output('Qsolar_per_area', 350.0, desc='solar heat rate absorbed per area', units='W/m**2')
        self.add_output('Qsolar_tot', 2902611.636, desc='total solar heat rate absorbed', units='W/m**2')
        self.add_output('heat_rate_per_pod', 519763.0, desc='heating due to a single pod', units='W')
        self.add_output('heat_rate_tot', 17671942.0, desc='heating rate due to all pods', units='W')
        self.add_output('Qradiated_per_area', 31.6, desc='heat radiated to the outside per area', units='W/m**2')
        self.add_output('Qradiated_tot', 106761066.5, desc='total heat radiated to the outside', units='W/m**2')
        self.add_output('area_viewing', 1074256.0, desc='effective area hit by sun', units='m**2')
        self.add_output('Qout_tot', 286900419.0, desc='total heat released via radiation', units='W')
        self.add_output('Qin_tot', 286900419.0, desc='total heat absorbed/added via pods and solar absorption', units='W')
        self.add_output('Q_resid', 0.0, desc='residual of Qin_tot and Qout_tot', units='W')

    def solve_nonlinear(self, params, unknowns, resids):
        self._clear_unknowns('flow_nozzle', unknowns)
        self._clear_unknowns('flow_bearings', unknowns)
        self._solve_flow_vars('flow_nozzle', params, unknowns)
        self._solve_flow_vars('flow_bearings', params, unknowns)
        # Q = mdot * cp * deltaT
        Qbearing = cu(unknowns['flow_bearings:out:W'],'lbm/s', 'kg/s') * cu(unknowns['flow_bearings:out:Cp'], 'Btu/lbm/degR', 'J/kg/K') * (cu(unknowns['flow_bearings:out:Tt'], 'degR', 'degK') - params['temp_boundary'])
        Qnozzle = cu(unknowns['flow_nozzle:out:W'], 'lbm/s', 'kg/s') * cu(unknowns['flow_nozzle:out:Cp'], 'Btu/lbm/degR', 'J/kg/K') * (cu(unknowns['flow_nozzle:out:Tt'], 'degR', 'degK') - params['temp_boundary'])
        unknowns['heat_rate_per_pod'] = Qnozzle + Qbearing
        unknowns['heat_rate_tot'] = unknowns['heat_rate_per_pod'] * params['n_pods']
        # Determine thermal resistance of outside via natural or forced convection
        # Prandtl # (Pr) = viscous diffusion rate / thermal diffusion rate = Cp * dyanamic viscosity / thermal conductivity
        # Pr << 1: thermal diffusivity dominates; Pr >> 1: momentum diffusivity dominates
        # SI units (https://mdao.grc.nasa.gov/publications/Berton-Thesis.pdf pg51)
        if params['temp_ambient'] < 400.0:
            unknowns['GrDelTL3'] = 4.178e19 * params['temp_ambient'] ** -4.639
            unknowns['Pr'] = 1.23 * params['temp_ambient'] ** -0.09685
            unknowns['k'] = 0.0001423 * params['temp_ambient'] ** 0.9138
        else:
            unknowns['GrDelTL3'] = 4.985e18 * params['temp_ambient'] ** -4.284
            unknowns['Pr'] = 0.59 * params['temp_ambient'] ** 0.0239
            unknowns['k'] = 0.0002494 * params['temp_ambient'] ** 0.8152
        # Grashof # (Gr) < 10^8: laminar; Gr > 10^9: turbulent
        unknowns['Gr'] = unknowns['GrDelTL3'] * (params['temp_boundary'] - params['temp_ambient']) * (2.0 * params['r_tube_outer']) ** 3
        # Rayleigh #: buoyancy driven flow (natural convection)
        unknowns['Ra'] = unknowns['Pr'] * unknowns['Gr']
        if unknowns['Ra'] <= 1e12: # valid in specific flow regime
            # Nusselt # (Nu) = convective heat transfer / conductive heat transfer
            unknowns['Nu'] = (0.6 + 0.387 * unknowns['Ra'] ** (1.0 / 6.0) / (1.0 + (0.559 / unknowns['Pr']) ** (9.0 / 16.0)) ** (8.0 / 27.0)) ** 2 # 3rd Ed. of Introduction to Heat Transfer by Incropera and DeWitt, equations (9.33) and (9.34) on page 465
        else:
            raise Exception('Rayleigh number outside of acceptable range.')
        unknowns['h'] = unknowns['k'] * unknowns['Nu'] / (2.0 * params['r_tube_outer']) # h = k * Nu / characteristic length
        unknowns['convection_area'] = pi * params['tube_len'] * 2.0 * params['r_tube_outer']
        unknowns['Qradiated_nat_convection_per_area'] = unknowns['h'] * (params['temp_boundary'] - params['temp_ambient'])
        unknowns['Qradiated_nat_convection_tot'] = unknowns['Qradiated_nat_convection_per_area'] * unknowns['convection_area']
        unknowns['area_viewing'] = params['tube_len'] * 2.0 * params['r_tube_outer'] # sun hits an effective rectangular cross section
        unknowns['Qsolar_per_area'] = (1.0 - params['reflectance']) * params['nn_incidence_factor'] * params['insolation']
        unknowns['Qsolar_tot'] = unknowns['Qsolar_per_area'] * unknowns['area_viewing']
        unknowns['radiating_area'] = unknowns['convection_area']
        unknowns['Qradiated_per_area'] = params['sb_const'] * params['emissivity'] * (params['temp_boundary'] ** 4 - params['temp_ambient'] ** 4) # P / A = SB * emmisitivity * (T ** 4 - To ** 4)
        unknowns['Qradiated_tot'] = unknowns['radiating_area'] * unknowns['Qradiated_per_area']
        unknowns['Qout_tot'] = unknowns['Qradiated_tot'] + unknowns['Qradiated_nat_convection_per_area']
        unknowns['Qin_tot'] = unknowns['Qsolar_tot'] + unknowns['heat_rate_tot']
        unknowns['Q_resid'] = abs(unknowns['Qout_tot'] - unknowns['Qin_tot'])

if __name__ == '__main__':
    from openmdao.core.group import Group
    from openmdao.core.problem import Problem
    from openmdao.drivers.scipy_optimizer import ScipyOptimizer
    from openmdao.components.param_comp import ParamComp
    from openmdao.components.constraint import ConstraintComp

    g = Group()
    p = Problem(root=g, driver=ScipyOptimizer())
    p.driver.options['optimizer'] = 'COBYLA'

    g.add('temp_boundary', ParamComp('T', 340.0))
    g.add('tube_wall', TubeWallTemp())
    g.connect('temp_boundary.T', 'tube_wall.temp_boundary')
    g.add('con', ConstraintComp('temp_boundary > 305.7', out='out'))
    g.connect('temp_boundary.T', 'con.temp_boundary')
    
    p.driver.add_param('temp_boundary.T', low=0.0, high=10000.0)
    p.driver.add_objective('tube_wall.Q_resid')
    p.driver.add_constraint('con.out')

    p.setup()

    g.tube_wall.params['flow_nozzle:in:Tt'] = 1710.0
    g.tube_wall.params['flow_nozzle:in:Pt'] = 0.304434211
    g.tube_wall.params['flow_nozzle:in:W'] = 1.08
    g.tube_wall.params['flow_bearings:in:W'] = 0.0
    g.tube_wall.params['r_tube_outer'] = 2.22504 / 2.0
    g.tube_wall.params['tube_len'] = 482803.0
    g.tube_wall.params['n_pods'] = 34
    g.tube_wall.params['temp_ambient'] = 305.6

    p.run()

    print '\nCompleted tube heat flux model calculations...\n'
    print 'Compress Q:             %g\nSolar Q:                %g\nRadiation Q:            %g\nConvection Q:           %g' % (g.tube_wall.unknowns['heat_rate_tot'], g.tube_wall.unknowns['Qsolar_tot'], g.tube_wall.unknowns['Qradiated_tot'], g.tube_wall.unknowns['Qradiated_nat_convection_tot'])
    print 'Equilibrium wall temp.: %g K or %g F' % (g.unknowns['temp_boundary.T'], cu(g.unknowns['temp_boundary.T'], 'degK', 'degF'))
    print 'Ambient temp.:          %g K or %g F' % (g.tube_wall.params['temp_ambient'], cu(g.tube_wall.params['temp_ambient'], 'degK', 'degF'))
    print 'Q out:                  %g W\nQ in:                   %g W\nError:                  %3.9f%%\n' % (g.tube_wall.unknowns['Qout_tot'], g.tube_wall.unknowns['Qin_tot'], (g.tube_wall.unknowns['Qout_tot'] - g.tube_wall.unknowns['Qin_tot']) / g.tube_wall.unknowns['Qout_tot'] * 100.0)

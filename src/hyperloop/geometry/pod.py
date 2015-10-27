from openmdao.core.group import Group

from inlet import InletGeom
from battery import Battery
from passenger_capsule import PassengerCapsule
from tube_structure import TubeStructural
from aero import Aero

class Pod(Group):
    def __init__(self):
        super(Pod, self).__init__()
        capsule = self.add('capsule', PassengerCapsule(), promotes=['n_rows', 'row_len'])
        tube = self.add('tube', TubeStructural(), promotes=['tube_P', 'tube_T', 'tube_r', 'tube_area', 'fill_area'])
        inlet = self.add('inlet', InletGeom(), promotes=['hub_to_tip', 'tube_area', 'bypass_area', 'cross_section', 'area_frontal'])
        battery = self.add('battery', Battery(), promotes=['time_mission', 'energy', 'cross_section'])
        aero = self.add('aero', Aero(), promotes=['rho', 'gross_thrust', 'net_force', 'coef_drag', 'velocity_capsule', 'area_frontal'])

if __name__ == '__main__':
    from openmdao.core.problem import Problem

    p = Problem(root=Pod())
    p.setup()
    p.run()

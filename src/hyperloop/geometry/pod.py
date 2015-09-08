from openmdao.core.group import Group

from inlet import InletGeom
from battery import Battery
from passenger_capsule import PassengerCapsule
from tube_structure import TubeStructural
from aero import Aero

class Pod(Group):
    def __init__(self):
        super(Pod, self).__init__()
        capsule = self.add('capsule', PassengerCapsule(), promotes=['n_rows', 'row_len', 'cross_section'])
        tube = self.add('tube', TubeStructural())
        inlet = self.add('inlet', InletGeom(), promotes=['hub_to_tip', 'area_bypass'])
        battery = self.add('battery', Battery(), promotes=['time_mission', 'energy'])
        aero = self.add('aero', Aero(), promotes=['rho', 'gross_thrust', 'net_force', 'coef_drag', 'velocity_capsule'])
        self.connect('cross_section', 'inlet.cross_section')
        self.connect('cross_section', 'battery.cross_section')
        self.connect('inlet.area_frontal', 'aero.area_frontal')

if __name__ == '__main__':
    from openmdao.core.problem import Problem

    p = Problem(root=Pod())
    p.setup()
    p.run()

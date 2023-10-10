import numpy
import trimesh
from stl_util import STL_Util

materials = {
    1: {'name': 'ABS', 'mass': 1.04},
    2: {'name': 'PLA', 'mass': 1.25},
    3: {'name': '3k CFRP', 'mass': 1.79},
    4: {'name': 'Plexiglass', 'mass': 1.18},
    5: {'name': 'Alumide', 'mass': 1.36},
    6: {'name': 'Aluminum', 'mass': 2.68},
    7: {'name': 'Brass', 'mass': 8.6},
    8: {'name': 'Bronze', 'mass': 9.0},
    9: {'name': 'Copper', 'mass': 9.0},
    10: {'name': 'Gold_14K', 'mass': 13.6},
    11: {'name': 'Gold_18K', 'mass': 15.6},
    12: {'name': 'Polyamide_MJF', 'mass': 1.01},
    13: {'name': 'Polyamide_SLS', 'mass': 0.95},
    14: {'name': 'Rubber', 'mass': 1.2},
    15: {'name': 'Silver', 'mass': 10.26},
    16: {'name': 'Steel', 'mass': 7.86},
    17: {'name': 'Titanium', 'mass': 4.41},
    18: {'name': 'Resin', 'mass': 1.2}
}

def main():
    lower_brace = trimesh.load('lower_2.stl')
    scene = trimesh.Scene(lower_brace)
    scene.add_geometry(geometry=lower_brace)
    print(scene.area/1000, "cm^2")
    print(scene.volume/1000, "cm^3")
    scene.moment_inertia
    # print(scene.density)
    # scene.show()

if __name__ == "__main__":
    main()
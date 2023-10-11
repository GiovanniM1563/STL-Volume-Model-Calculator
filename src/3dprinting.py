import math
import struct
import sys
import numpy as np
import trimesh
from stl_util import STL_Util
from model import Model

"""
# Load your 3D mesh (replace 'your_mesh.obj' with your mesh file)
mesh = trimesh.load_mesh('../res/lower_1.stl')

# Find the bottommost vertex
bottommost_vertex = find_bottommost_vertex(mesh)

print("Bottommost Vertex:", bottommost_vertex)
"""

def main():
    denture = Model({'Resin': 1.2})
    stl_helper = STL_Util()

    try:
        if len(sys.argv) == 1:
            print("Define model to calculate volume eg: python measure_volume.py torus.stl")
        else:
            denture.reference = sys.argv[1]
            stl_helper.read_3d_model(path=denture.reference)
            denture.triangles_sum = stl_helper.get_sum_triangles(ref=denture.reference) 
            denture.mesh = trimesh.load_mesh(denture.reference)

            stl_helper.calculate_geometric_attributes(model=denture)
            print('--------- Denture 3D Model Details ---------')
            print('Surface Area:', denture.surface_area, 'cm^2')
            print('Volume:', denture.volume)
            print('Density:', denture.density)
            print('Mass:', denture.mass)
    except Exception as e:
        print(e)


    # lower_brace = trimesh.load('lower_2.stl')
    # scene = trimesh.Scene(lower_brace)
    # scene.add_geometry(geometry=lower_brace)
    # print(scene.area/1000, "cm^2")
    # print(scene.volume/1000, "cm^3")
    # scene.moment_inertia
    # # print(scene.density)
    # # scene.show()

if __name__ == "__main__":
    main()
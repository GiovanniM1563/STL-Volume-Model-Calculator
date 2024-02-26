import math
import struct
import sys
import traceback
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
        # if len(sys.argv) = 1:
        #     print("Define model to calculate volume eg: python measure_volume.py torus.stl")
        # else:

        # Change 
        
        denture.reference = '../res/Upper_1.stl'
        stl_helper.read_3d_model(path=denture.reference)
        denture.triangles_sum = stl_helper.get_sum_triangles(ref=denture.reference) 
        
        # Load the mesh before beginning complex calculations
        denture.mesh = trimesh.load_mesh(denture.reference)
        # denture.mesh.show()

        stl_helper.calculate_geometric_attributes(model=denture)
        stl_helper.find_top_most_vertex(model=denture)
        print('--------- Denture 3D Model Details ---------')
        print('Surface Area using Linear Algebra:', denture.surface_area, 'cm^2')
        # print('Surface Area using Trimesh:', denture.mesh.area/1000, 'cm^2')
        print('Volume using Tetrahedron Formula:', denture.volume, 'cm^3')
        print('Volume using Trimesh:', denture.mesh.volume/1000, 'cm^3')
        print('Density:', denture.density, 'g/cm^3')
        print('Mass:', denture.mass, 'g')
    except:
        print(traceback.format_exc())

if __name__ == "__main__":
    main()
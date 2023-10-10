import trimesh
import numpy as np

# Load your 3D mesh (replace 'your_mesh.obj' with your mesh file)
mesh = trimesh.load_mesh('LowerBase.stl')

# Define the normal vector
normal = np.array([0, 0, 1])

# Define the desired thickness of the slice in the Z-axis direction
thickness = 5 # Adjust this value as needed, Higher = Finer/thinner 

# Set the Z-coordinate of the slicing plane origin to control the thickness
point_on_plane = np.array([0, 0, thickness])

# Perform the slice using the defined plane
sliced_mesh = trimesh.intersections.slice_mesh_plane(mesh, plane_normal=normal, plane_origin=point_on_plane)

# Visualize the sliced mesh
sliced_mesh.show()
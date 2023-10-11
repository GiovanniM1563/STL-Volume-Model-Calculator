import trimesh
import numpy as np

# Load your 3D mesh (replace 'your_mesh.obj' with your mesh file)
mesh = trimesh.load_mesh('LowerBase.stl')

slices = []

lowest_z = find_bottommost_vertex(mesh)

# Define the normal vector
normal = np.array([0, 0, 1])

# Define the desired number of slices and thickness of each slice
num_slices = 100  # Adjust the number of slices as needed
slice_thickness = 0.1  # Adjust the slice thickness as needed

# Iterate through different heights and create slices
for i in range(num_slices):
    # Calculate the Z-coordinate of the slicing plane origin for the current slice
    z_coord = i * -slice_thickness
    
    if not np.allclose(z_coord, lowest_z):
        # Set the Z-coordinate of the slicing plane origin
        point_on_plane = np.array([0, 0, z_coord])
            
        # Perform the slice using the defined plane
        sliced_mesh = trimesh.intersections.slice_mesh_plane(mesh, plane_normal=normal, plane_origin=point_on_plane, face_index=None, cap=True, cached_dots=None,  engine=None )
        
        # Visualize the current slice
        # sliced_mesh.show()

        slices.append(sliced_mesh)
    else:
        break

surface_areas = []


for slice in slices:
    surface_area = slice.area
    surface_areas.append(surface_area)

print(max(surface_areas))
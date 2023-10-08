import trimesh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Code finds the top most plane of the mesh. useful for slicing as you can go all the way to the bottom.  

# Load your 3D mesh (replace 'your_mesh.obj' with your mesh file)
mesh = trimesh.load_mesh('LowerBase.stl')

# Compute vertex normals of the mesh
normals = mesh.vertex_normals

# Determine the direction that represents the "top" of the object
desired_normal_direction = [0.0, 0.0, 1.0]  # Specify the desired normal direction

# Define a tolerance to consider points on the "top" surface
tolerance = 0.01  # Adjust the tolerance as needed

# Collect the indices of vertices with normals pointing towards the desired normal direction
dot_products = np.dot(normals, desired_normal_direction)
top_surface_indices = np.where(dot_products > 1 - tolerance)[0]

# Extract the vertices corresponding to the top surface
top_surface_vertices = mesh.vertices[top_surface_indices]

# Create a plane that is perpendicular to the Z-axis
plane_normal = [0.0, 0.0, 1.0]  # Normal vector along the Z-axis

# Compute the centroid of the top surface vertices to use as the point on the plane
plane_origin = np.mean(top_surface_vertices, axis=0)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the top surface vertices
ax.scatter(top_surface_vertices[:, 0], top_surface_vertices[:, 1], top_surface_vertices[:, 2], c='b', marker='o', label='Top Surface Vertices')

# Plot the plane
xx, yy = np.meshgrid(np.linspace(top_surface_vertices[:, 0].min(), top_surface_vertices[:, 0].max(), 10),
                     np.linspace(top_surface_vertices[:, 1].min(), top_surface_vertices[:, 1].max(), 10))
zz = (-plane_normal[0] * xx - plane_normal[1] * yy - plane_origin.dot(plane_normal)) / plane_normal[2]
ax.plot_surface(xx, yy, zz, color='r', alpha=0.5, label='Fitted Plane')

# Set labels and legend
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

# Show the plot
plt.show()
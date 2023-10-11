import trimesh
import numpy as np

def find_bottommost_vertex(mesh):
    """
    Find the bottommost vertex of a trimesh object.

    Parameters:
        mesh (trimesh.base.Trimesh): The input trimesh object.

    Returns:
        numpy.ndarray: The coordinates of the bottommost vertex.
    """
    # Get the vertices of the mesh
    vertices = mesh.vertices

    # Find the vertex with the lowest Z-coordinate (assuming Z is the vertical axis)
    bottommost_vertex = vertices[np.argmin(vertices[:, 2])]

    return bottommost_vertex

# Load your 3D mesh (replace 'your_mesh.obj' with your mesh file)
mesh = trimesh.load_mesh('../res/lower_1.stl')

# Find the bottommost vertex
bottommost_vertex = find_bottommost_vertex(mesh)

print("Bottommost Vertex:", bottommost_vertex)
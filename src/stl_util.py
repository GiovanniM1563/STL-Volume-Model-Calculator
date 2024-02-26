import math
import struct
import sys
import traceback
from io import BufferedReader

import trimesh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from trimesh import creation, transformations
from scipy.spatial.transform import Rotation


class STL_Util:
    """
    STL_Util

    Utility class for 3D Model STL calculations

    Calculates
        1. Surface Area using numpy
        2. Volume using Tetrahedron volume formula
        3. Density 
        4. Mass of 3D Model based off material weight
    """
    def __init__(self) -> None:
        self.fb = []  # debug list
        self.file = None

    def volume_of_triangle(self, p1, p2, p3) -> float:
        """
        Calculate volume for the 3D mesh using Tetrahedron volume
        based on: http://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
        """
        # 0, 1, 2 = X, Y, Z
        v321 = p3[0] * p2[1] * p1[2]
        v231 = p2[0] * p3[1] * p1[2]
        v312 = p3[0] * p1[1] * p2[2]
        v132 = p1[0] * p3[1] * p2[2]
        v213 = p2[0] * p1[1] * p3[2]
        v123 = p1[0] * p2[1] * p3[2]
        return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)
    
    def area_of_triangle(self, p1, p2, p3) -> float:
        """
        Calculates volume for 3D mesh by summing total area of sum triangles
        based on: https://stackoverflow.com/questions/26312570/calculate-surface-area-of-a-3d-mesh

        Parameters:
            p1: Triangle edge 
            p2: 
            p3: 
        
        Returns:
            float: Area of triangle face
        """
        # 0, 1, 2 = X, Y, Z
        ax = p2[0] - p1[0]
        ay = p2[1] - p1[1]
        az = p2[2] - p1[2]
        bx = p3[0] - p1[0]
        by = p3[1] - p1[1]  
        bz = p3[2] - p1[2]
        cx = ay*bz - az*by
        cy = az*bx - ax*bz
        cz = ax*by - ay*bx
        v = np.matrix([[ax - cx, ay - cy, az - cz], [bx - cx, by - cy, bz - cz]])
        return 0.5 * np.sqrt(np.abs(np.linalg.det(v*v.T)))
    
    def find_bottom_most_vertex(self, model) -> np.ndarray:
        """
        Find the bottommost vertex of a trimesh object.

        Parameters:
            mesh (trimesh.base.Trimesh): The input trimesh object.

        Returns:
            numpy.ndarray: The coordinates of the bottommost vertex.
        """
        # Get the vertices of the mesh
        vertices = model.mesh.vertices

        # Find the vertex with the lowest Z-coordinate (assuming Z is the vertical axis)
        bottom_most_vertex = vertices[np.argmin(vertices[:, 2])]
        return bottom_most_vertex

    def find_top_most_vertex(self, model) -> np.ndarray:
        # Load your 3D mesh (replace 'your_mesh.obj' with your mesh file)
        # mesh = trimesh.load_mesh('LowerBase.stl')
        slices = []
        lowest_z = self.find_bottom_most_vertex(model)[2]

        # Define the desired number of slices and thickness of each slice
        slice_thickness = 1  # Adjust the slice thickness as needed
        num_slices = 90  # Adjust the number of slices as needed

        # User-defined angle for plane orientation in degrees
        slice_angle_degrees = 45 # Adjust the angle as needed

        # Convert the user-defined angle to radians
        slice_angle_radians = np.radians(slice_angle_degrees)

        # Define the normal vector for the initial slicing plane
        normal = np.array([0, 0, 1])

        # Iterate through different heights and create slices
        for i in range(num_slices):
            # Calculate the Z-coordinate of the slicing plane origin for the current slice
            z_coord = i * -slice_thickness

            if not np.allclose(z_coord, lowest_z):
                # Set the Z-coordinate of the slicing plane origin
                point_on_plane = np.array([0, 0, z_coord])

                # Create a rotation matrix for the user-defined angle
                rotation_matrix = Rotation.from_euler('x', slice_angle_radians, degrees=False)

                # Rotate the normal vector to orient the plane
                rotated_normal = rotation_matrix.apply(normal)

                # Perform the slice using the defined plane
                sliced_mesh = trimesh.intersections.slice_mesh_plane(
                    model.mesh,
                    plane_normal=rotated_normal,
                    plane_origin=point_on_plane,
                    face_index=None,
                    cap=True,
                    cached_dots=None,
                    engine=None
                )

                # Create a slice on the opposite side by rotating the plane 180 degrees
                reversed_rotation_matrix = Rotation.from_euler('x', np.pi, degrees=False)
                reversed_normal = reversed_rotation_matrix.apply(rotated_normal)

                # Adjust the Z-coordinate to add the slice thickness
                reversed_point_on_plane = point_on_plane + np.array([0, 0, slice_thickness])

                # Perform the opposite slice
                reversed_sliced_mesh = trimesh.intersections.slice_mesh_plane(
                    sliced_mesh,
                    plane_normal=reversed_normal,
                    plane_origin=reversed_point_on_plane,
                    face_index=None,
                    cap=True,
                    cached_dots=None,
                    engine=None
                )

                slices.append(reversed_sliced_mesh)
            else:
                break

        total_area = 0
        slice_area = {}
        garbage = set()
        name = trimesh.Scene()
        for slice in slices:
            for i in np.arange(len(slice.triangles)):
                x1, y1, z1 = slice.triangles[i][0][0], slice.triangles[i][0][1], slice.triangles[i][0][2]
                x2, y2, z2 = slice.triangles[i][0][0], slice.triangles[i][0][1], slice.triangles[i][0][2]
                x3, y3, z3 = slice.triangles[i][0][0], slice.triangles[i][0][1], slice.triangles[i][0][2]
                x = x1 + x2 + x3
                y = y1 + y2 + y3
                z = z1 + z2 + z3
                key = (x,y,z)
                # print(key)
                if key not in slice_area:
                    p1 = slice.triangles[i][0]
                    p2 = slice.triangles[i][1]
                    p3 = slice.triangles[i][2]
                    vertex_area = self.area_of_triangle(p1, p2, p3)
                    slice_area[key] = [p1, p2, p3, vertex_area]
                    total_area += vertex_area
                else:
                    garbage.add(key)

            # print('t area:', total_area / 1000, 'cm^2')
        name.add_geometry(slices)
        name.show()
        print('Final model Surface Area:', total_area / 1000, 'cm^2')
        print('seen this many duplicate points: ', len(garbage))

    # Documentation required here to understand what this is doing
    def unpack(self, sig, byte):
        s = self.file.read(byte)
        return struct.unpack(sig, s)
    
    def read_3d_model(self, path) -> None:
        self.file = open(file=path, mode='rb')
        self.file.seek(self.file.tell() + 80)

    # Might need to refactor this to make it more readable
    def read_triangle(self, model) -> list:
        n = self.unpack(sig="<3f", byte=12)
        p1 = self.unpack(sig="<3f", byte=12)
        p2 = self.unpack(sig="<3f", byte=12)
        p3 = self.unpack(sig="<3f", byte=12)
        self.unpack("<h", 2)
        model.normals.append(n)
        model.points.append(p1)
        model.points.append(p2)
        model.points.append(p3)
        # model.byte_count.append(b[0])
        
        # TODO: Refactor this
        l = len(model.points)
        model.triangles.append((l, l + 1, l + 2))
        return [p1, p2, p3, n]

    def get_sum_triangles(self, ref) -> bytes:
        length = struct.unpack("@i", self.file.read(4))
        return length[0]
    
    def convert_cm3_to_in3(self, v) -> float:
        return v * 0.0610237441

    def calculate_geometric_attributes(self, model) -> None:
        try:
            # We want to slice model, so this for loop will change once we implement slicing
            for _ in range(model.triangles_sum):
                edge = self.read_triangle(model=model)
                # f.write('{}{}{}{}\n'.format(edge[0], edge[1], edge[2], edge[3]))
                model.surface_area += self.area_of_triangle(edge[0], edge[1], edge[2])
                model.volume += self.volume_of_triangle(edge[0], edge[1], edge[2])
                # print(model.volume, edge[0], edge[1], edge[2])
            model.volume /= 1000
            model.surface_area /= 1000
            model.mass = model.material['Resin'] * model.volume
            model.density = model.mass / model.volume

        except:
            print(traceback.format_exc())

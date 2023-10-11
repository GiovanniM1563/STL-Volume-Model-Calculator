import math
import struct
import sys
import traceback
from io import BufferedReader

import trimesh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

    # Calculate volume for the 3D mesh using Tetrahedron volume
    # based on: http://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    def volume_of_triangle(self, p1, p2, p3) -> float:
        # 0, 1, 2 = X, Y, Z
        v321 = p3[0] * p2[1] * p1[2]
        v231 = p2[0] * p3[1] * p1[2]
        v312 = p3[0] * p1[1] * p2[2]
        v132 = p1[0] * p3[1] * p2[2]
        v213 = p2[0] * p1[1] * p3[2]
        v123 = p1[0] * p2[1] * p3[2]
        return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)
    
    # Refactor and use trimesh area as it is more accurate
    def area_of_triangle(self, p1, p2, p3) -> float:
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
        return 0.5 * np.sqrt(np.linalg.det(v*v.T))
        # lower_brace = trimesh.load('lower_2.stl')
        # scene = trimesh.Scene(lower_brace)
        # scene.add_geometry(geometry=lower_brace)
        # print(scene.area/1000, "cm^2")
        # print(scene.volume/1000, "cm^3")
        # scene.moment_inertia
        # # print(scene.density)
        # # scene.show()
    
    def find_bottommost_vertex(self, model) -> np.ndarray:
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
        bottommost_vertex = vertices[np.argmin(vertices[:, 2])]

        return bottommost_vertex

    def find_topmost_vertex(self, model) -> np.ndarray:
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
                model.surface_area += self.area_of_triangle(edge[0], edge[1], edge[2])
                model.volume += self.volume_of_triangle(edge[0], edge[1], edge[2])
                # print(model.volume, edge[0], edge[1], edge[2])
            
            model.volume /= 1000
            model.surface_area /= 1000
            model.mass = model.material['Resin'] * model.volume
            model.density = model.mass / model.volume

        except:
            print(traceback.format_exc())

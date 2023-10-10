import math
import struct
import sys
import numpy as np
import trimesh

print('Choose desired print material of STL file below:')
print('1 = ABS')
print('2 = PLA')
print('3 = 3k CFRP')
print('4 = Plexiglass')
print('5 = Alumide')
print('6 = Aluminum')
print('7 = Brass')
print('8 = Bronze')
print('9 = Copper')
print('10 = Gold_14K')
print('11 = Gold_18K')
print('12 = Polyamide_MJF')
print('13 = Polyamide_SLS')
print('14 = Rubber')
print('15 = Silver')
print('16 = Steel')
print('17 = Titanium')
print('18 = Resin')

material = input('Enter the number corresponding to the desired print material: ')

# Validate material input
try:
    material = int(material)
    if material < 1 or material > 18:
        material = 1
except ValueError:
    material = 1

class Model:
    def __init__(self, material) -> None:
        self.normals = []
        self.points = []
        self.triangles = []
        self.byte_count = []
        self.material = material
        self.volume = 0.0
        self.mass = 0.0
        self.surface_area = 0.0
        self.density = 0.0

class STL_Util:
    def __init__(self) -> None:
        self.fb = []  # debug list

    def resetVariables(self) -> None:
        self.normals = []
        self.points = []
        self.triangles = []
        self.bytecount = []
        self.fb = []  # debug list

    # Calculate volume for the 3D mesh using Tetrahedron volume
    # based on: http://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    def signedVolumeOfTriangle(self, p1, p2, p3) -> float:
        # 0, 1, 2 = X, Y, Z
        v321 = p3[0] * p2[1] * p1[2]
        v231 = p2[0] * p3[1] * p1[2]
        v312 = p3[0] * p1[1] * p2[2]
        v132 = p1[0] * p3[1] * p2[2]
        v213 = p2[0] * p1[1] * p3[2]
        v123 = p1[0] * p2[1] * p3[2]
        return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)
    
    def area_of_triangle(self, p1, p2, p3) -> float:
        # 0, 1, 2 = X, Y, Z
        '''
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
        return 0.5 * math.sqrt(np.linalg.det(v*v.T))
        '''
        # lower_brace = trimesh.load('lower_2.stl')
        # scene = trimesh.Scene(lower_brace)
        # scene.add_geometry(geometry=lower_brace)
        # print(scene.area/1000, "cm^2")
        # print(scene.volume/1000, "cm^3")
        # scene.moment_inertia
        # # print(scene.density)
        # # scene.show()
        pass

    def unpack(self, sig, l) -> tuple:
        s = self.f.read(l)
        self.fb.append(s)
        return struct.unpack(sig, s)

    def read_triangle(self) -> list:
        n = self.unpack("<3f", 12)
        p1 = self.unpack("<3f", 12)
        p2 = self.unpack("<3f", 12)
        p3 = self.unpack("<3f", 12)
        b = self.unpack("<h", 2)
        self.normals.append(n)
        l = len(self.points)
        self.points.append(p1)
        self.points.append(p2)
        self.points.append(p3)
        self.triangles.append((l, l + 1, l + 2))
        self.bytecount.append(b[0])
        return [p1, p2, p3]

    def read_length(self) -> bytes:
        length = struct.unpack("@i", self.f.read(4))
        return length[0]

    def read_header(self) -> int:
        self.f.seek(self.f.tell() + 80)

    def cm3_To_inch3Transform(self, v) -> float:
        return v * 0.0610237441

    def calculateMassCM3(self, total_volume) -> float:
        if self.material in self.materials:
            material_mass = self.materials[self.material]['mass']
            return total_volume * material_mass
        return 0

    def calculateVolume(self, material, file, unit):
        print(file)
        self.resetVariables()
        total_volume = 0
        total_area = 0
        total_mass = 0
        density = 0
        try:
            self.f = open(file, "rb")
            self.read_header()
            l = self.read_length()
            print("total triangles:", l)
            try:
                i = 0
                while True:
                    edge = self.read_triangle()
                    i += 1
                    # print("reading points", edge)
                    total_volume += self.signedVolumeOfTriangle(edge[0], edge[1], edge[2])
                    # print("After Volume Call", edge)
                    total_area += self.area_of_triangle(edge[0], edge[1], edge[2])
                    # print("After Triangle Call", edge)
            except Exception as e:
                print("ERROR: ", e)
                print(i)

            total_volume /= 1000
            total_area /= 1000
            total_mass = self.calculateMassCM3(material=material, total_volume=total_volume)
            density = total_mass / total_volume

            if total_mass <= 0:
                print('Total mass could not be calculated')
            else:
                print('Total mass:', total_mass, 'g')

                if unit == "cm":
                    print("Volume:", total_volume, "cm^3")
                    print("Surface Area:", total_area, "cm^2")
                    print("Density:", density, "g/cm^3")
                else:
                    total_volume = self.cm3_To_inch3Transform(total_volume)
                    total_area = self.cm3_To_inch3Transform(total_area)
                    density = self.cm3_To_inch3Transform(density)
                    print("Volume:", total_volume, "inch^3")
                    print("Surface Area:", total_area, "inch^2")
                    print("Density:", density, "g/in^3")
        except Exception as e:
            print(e)
        return total_volume

def main():
    if len(sys.argv) == 1:
        print("Define model to calculate volume ej: python measure_volume.py torus.stl")
    else:
        stl_util = STL_Util()
        if(len(sys.argv) > 2 and sys.argv[2] == "inch"):
            stl_util.calculateVolume(material=material, file=sys.argv[1], unit="inch")
        else:
            stl_util.calculateVolume(material=material, file=sys.argv[1], unit="cm")

if __name__ == '__main__':
    main()


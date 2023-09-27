#! /usr/bin/env python3
'''
VOLUME CALCULATION STL binary MODELS
Author: Mar Canet (mar.canet@gmail.com) - september 2012
Description: Added input call for print material (ABS or PLA), added print of object mass, made Python3 compatible, changed tabs for spaces
Notes: Material Mass Source is https://www.toybuilderlabs.com/blogs/news/13053117-filament-volume-and-length

Contributors: Saijin_Naib (Synper311@aol.com)
'''

import math
import struct
import sys

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


class STLUtils:
    def resetVariables(self):
        self.normals = []
        self.points = []
        self.triangles = []
        self.bytecount = []
        self.fb = []  # debug list

    # Calculate volume for the 3D mesh using Tetrahedron volume
    # based on: http://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    def signedVolumeOfTriangle(self, p1, p2, p3):
        # 0 = X
        # 1 = Y
        # 2 = Z
        v321 = p3[0] * p2[1] * p1[2]
        v231 = p2[0] * p3[1] * p1[2]
        v312 = p3[0] * p1[1] * p2[2]
        v132 = p1[0] * p3[1] * p2[2]
        v213 = p2[0] * p1[1] * p3[2]
        v123 = p1[0] * p2[1] * p3[2]
        return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)
    
    def area_of_triangle(self, p1, p2, p3):
        # 0 = X
        # 1 = Y
        # 2 = Z
        ax = p2[0] - p1[0]
        ay = p2[1] - p1[1]
        az = p2[2] - p1[2]
        bx = p3[0] - p1[0]
        by = p3[1] - p1[1]
        bz = p3[2] - p1[2]
        cx = ay*bz - az*by
        cy = az*bx - ax*bz
        cz = ax*by - ay*bx
        return 0.5 * math.sqrt(cx*cx + cy*cy + cz*cz)

    def unpack(self, sig, l):
        s = self.f.read(l)
        self.fb.append(s)
        return struct.unpack(sig, s)

    def read_triangle(self):
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

    def read_length(self):
        length = struct.unpack("@i", self.f.read(4))
        return length[0]

    def read_header(self):
        self.f.seek(self.f.tell() + 80)

    def cm3_To_inch3Transform(self, v):
        return v * 0.0610237441

    def calculateMassCM3(self, totalVolume):
        if material in materials:
            material_mass = materials[material]['mass']
            return totalVolume * material_mass
        return 0

    def calculateVolume(self, infilename, unit):
        print(infilename)
        self.resetVariables()
        totalVolume = 0
        totalArea = 0
        totalMass = 0
        density = 0
        try:
            self.f = open(infilename, "rb")
            self.read_header()
            l = self.read_length()
            print("total triangles:", l)
            try:
                while True:
                    edge = self.read_triangle()
                    totalVolume += self.signedVolumeOfTriangle(edge[0], edge[1], edge[2])
                    totalArea += self.area_of_triangle(edge[0], edge[1], edge[2])
            except Exception as e:
                print("End calculate triangles volume")
            totalVolume /= 1000
            totalArea /= 1000
            totalMass = self.calculateMassCM3(totalVolume)
            density = totalMass / totalVolume

            if totalMass <= 0:
                print('Total mass could not be calculated')
            else:
                print('Total mass:', totalMass, 'g')

                if unit == "cm":
                    print("Volume:", totalVolume, "cm^3")
                    print("Surface Area:", totalArea, "cm^3")
                    print("Density:", density, "g/cm^3")
                else:
                    totalVolume = self.cm3_To_inch3Transform(totalVolume)
                    totalArea = self.cm3_To_inch3Transform(totalArea)
                    density = self.cm3_To_inch3Transform(density)
                    print("Volume:", totalVolume, "inch^3")
                    print("Surface Area:", totalArea, "inch^3")
                    print("Density:", density, "g/in^3")
        except Exception as e:
            print(e)
        return totalVolume


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("Define model to calculate volume ej: python measure_volume.py torus.stl")
    else:
        mySTLUtils = STLUtils()
        if(len(sys.argv) > 2 and sys.argv[2] == "inch"):
            mySTLUtils.calculateVolume(sys.argv[1], "inch")
        else:
            mySTLUtils.calculateVolume(sys.argv[1], "cm")

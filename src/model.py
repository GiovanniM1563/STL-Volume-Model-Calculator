from io import BufferedReader

class Model:
    def __init__(self, material) -> None:
        self.normals = []
        self.points = []
        self.triangles = []
        self.triangles_sum = 0
        self.byte_count = []
        self.material = material
        self.mesh = None
        """Contains 3D Model based off STL File"""
        self.sliced_mesh = []
        """Contains per slice {trimesh.Mesh -> [surface_area, volume, mass, density]"""
        self.surface_area = 0.0
        """Returns total surface area of 3D Model"""
        self.volume = 0.0
        """Returns sum total of volume using tetrahedron formula for 3D Model"""
        self.mass = 0.0
        """Returns calculated mass based off material weight * volume for 3D Model"""
        self.density = 0.0
        """Returns D = M/V"""
        self.reference = ''
        """Returns a reference to the filepath of 3D Model"""


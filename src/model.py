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
        self.volume = 0.0
        self.mass = 0.0
        self.surface_area = 0.0
        self.density = 0.0
        self.reference = None

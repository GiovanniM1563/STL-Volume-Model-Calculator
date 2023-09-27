import numpy
import trimesh

# pip install "pyglet<2
# pip install trimesh
# pip install numpy
# 

def main():
    lower_brace = trimesh.load('LowerBase.stl')
    scene = trimesh.Scene(lower_brace)
    scene.add_geometry(geometry=lower_brace)
    print(scene.area)
    print(scene.volume)
    scene.moment_inertia
    # print(scene.density)
    scene.show()

if __name__ == "__main__":
    main()
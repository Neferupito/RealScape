import bpy
import os
import mathutils
import struct
import numpy as np
import subprocess
import ctypes
import os
bl_info = {
    "name": "RealScape",
    "blender": (4, 3, 0),  # Update this to your version
    "category": "Object",
    "author": "François DALL'ASTA",
    "version": (1, 0),
    "description": "Addon to simulate erosion and deposition on a landscape to make it more realistic",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
}

r=128
c=128
iter_max = 10000
rain = 0.00001
power_coef = 0.5
erosion_strength = 0.0001
deposition_strength = 0.001
diffusion_coef = 0.000000000


# This function will be executed when the button is clicked
def my_function():
    print("Button clicked! Running the Python script...")

    export_to_xyz(os .path.join("inout", "in.xyz"))

    with open(os .path.join("inout", "param.txt"), 'w') as file:
        # Write the values to the file, separated by spaces
        file.write(' '.join(map(str, [r, c, iter_max, rain, power_coef, erosion_strength, deposition_strength, diffusion_coef])) + '\n')

    os.system("docker run -v inout:/RealScape/inout -it realscape")

    os.system("docker exec realscape ./simu")

    xyz2mesh(os .path.join("inout", "out.bin"))


def export_to_xyz(file_path):
    # Ensure the active object is a mesh
    obj = bpy.context.active_object
    if obj and obj.type == 'MESH':
        # Get the mesh data
        mesh = obj.data
        
        # Open the file for writing
        with open(file_path, 'w') as file:
            # Iterate over the vertices and write their coordinates to the file
            for vertex in mesh.vertices:
                # Write the x, y, z coordinates of each vertex
                file.write(f"{vertex.co.x:.8f} {vertex.co.y:.8f} {vertex.co.z:.8f}\n")
        
        print(f"Mesh exported to {file_path}")
    else:
        print("Please select a mesh object.")

def xyz2mesh(file_path):
    # Open the binary file
    with open(file_path, "rb") as file:
        # Read the entire content into bytes
        binary_data = file.read()

    # Unpack the data assuming it is a sequence of 4-byte floats
    num_floats = len(binary_data) // 4  # Calculate how many floats are in the file
    floats = struct.unpack(f"{num_floats}f", binary_data)

    if len(floats) % 3 != 0:
        print("Error: The number of floats in the file is not a multiple of 3.")
        exit(1)

    x = floats[0::3]
    y = floats[1::3]
    z = floats[2::3]

    Nx = np.unique(x).shape[0]
    Ny = np.unique(y).shape[0]

    if Nx * Ny != len(x):
        print("Error: x y z points are not gridded")
        exit(1)

    # Create a new mesh and object
    mesh = bpy.data.meshes.new(name="CubeMesh")
    obj = bpy.data.objects.new("Cube", mesh)

    # Link the object to the scene
    scene = bpy.context.scene
    scene.collection.objects.link(obj)

    vertices = []
    faces = []
    # Define vertices
    for i in range(len(x)):
        # Define the vertices (8 points for a cube)
        vertices.append(mathutils.Vector((x[i], y[i], z[i])))  
    # Define faces
    for i in range(Ny-1):
        for j in range(Nx-1):
            faces.append((i*Ny+j, i*Ny+j+1, (i+1)*Ny+j+1, (i+1)*Ny+j))

    # Create the mesh from the data
    mesh.from_pydata(vertices, [], faces)

    # Update the mesh to reflect changes
    mesh.update()

    # Set the object to be active and selected
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)


class SimpleButtonPanel(bpy.types.Panel):
    bl_label = "Simple Button Panel"
    bl_idname = "OBJECT_PT_simple_button"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Tool'  # This defines which tab the panel will be in ("Tool", "Item", etc.)

    def draw(self, context):
        layout = self.layout

        # Button to execute the function
        layout.operator("object.simple_button_operator", text="Click to Run Script")

# Define the operator class that runs when the button is clicked
class SimpleButtonOperator(bpy.types.Operator):
    bl_idname = "object.simple_button_operator"
    bl_label = "Run Script"

    def execute(self, context):
        my_function()  # Call the function when the button is clicked
        return {'FINISHED'}
    


# Register the classes and the operator
def register():
    bpy.utils.register_class(SimpleButtonPanel)
    bpy.utils.register_class(SimpleButtonOperator)

def unregister():
    bpy.utils.unregister_class(SimpleButtonPanel)
    bpy.utils.unregister_class(SimpleButtonOperator)

if __name__ == "__main__":
    register()

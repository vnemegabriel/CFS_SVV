import Nekpy.NekMesh as Nekmesh
from Nekpy.NekMesh import Node
from Nekpy.NekMesh import ElmtConfig, Element
from Nekpy.LibUtilities import ShapeType
# Modules
from NekPy.NekMesh import Module, ModuleType
from Nekpy.NekMesh import InputModule, OutputModule, ProcessModule

import wurlitzer

# Create a new, empty mesh.
mesh = NekMesh.Mesh()

# Now create an output module. Notice that we can pass the output filename
# to the module by specifying the 'outfile' keyword argument.
xml_writer = OutputModule.Create('xml', mesh, outfile='single_tet_mesh.xml')

# Similarly, let's create a module to check each element's Jacobian.
jac_module = ProcessModule.Create('jac', mesh, list=True)

class MyInputModule(InputModule):
    def __init__(self, mesh):
        # Call InputModule's constructor
        super().__init__(mesh)
    
    def Process(self):
        # Set the mesh expansion and space dimensions
        self.mesh.expDim = 3
        self.mesh.spaceDim = 3 # Or more dimensions, never less than expDim

    tet_nodes = [
        Node(0, 0.0, 0.0, 0.0), 
        Node(1, 1.0, 0.0, 0.0), 
        Node(2, 0.5, 1.0, 0.0),
        Node(3, 0.5, 0.5, 1.0)
        ]

    # Construct a config. for a simple linear tetrahedron
    # In this statementm the arguments to Elmt config are:
    # 1. ShapeType (tetrahedron)
    # 2. Order of element
    # 3. Whether the elements has high-order edge nodes and face interior nodes 
    # (both boolean values)

    tetra_config = ElmtConfig(ShapeType.Tetrahedron, 1, False, False)
    # Construct tetrahedron
    tetra = Element(tetra_config, tet_nodes, [0])
    # Puting the tetrahedron into the mesh
    mesh.element[3].append(Element.Create(tetra_config, tet_nodes, [ 0 ]))

    # Adding boundaries to the domain. In this case we need 4 triangles.
    # We can extract the nodes for each triangle from 
    tri_config = ElmtConfig(ShapeType.Triangle, 1, False, False)
    tri_faces = [(0,1,2), (0,1,3), (1,2,3), (0,2,3)]

    # Now constructing each face
    for i in range(4):
        mesh.element[2].append
        (Element.Create(
            tri_config, 
            [ tet_nodes[j] for j in tri_faces[i] ],
            [ i+1] ))

    # Call the Module functions to create all of the edges, faces and
    # composites for the mesh.
    self.ProcessVertices()
    self.ProcessEdges()
    self.ProcessFaces()
    self.ProcessElements()
    self.ProcessComposites()

# Register this module with our factory.
Module.Register(ModuleType.Input, 'myinputmodule', MyInputModule)

mesh.verbose = True
input_module = InputModule.Create('myinputmodule', mesh)

with wurlitzer.sys_pipes():
    # Run each module in succession by calling its Process function.
    input_module.Process()
    jac_module.Process()
    xml_writer.Process()

with open('single_tet_mesh.xml', 'r') as f:
    print(f.read())

    
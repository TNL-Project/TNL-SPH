import numpy as np
import vtk

import sys
import argparse

# ------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument("-input")
parser.add_argument("-output")

args = parser.parse_args()
inputFileName = args.input
outputFileName = args.output

if args.input == None:
    print("Error: Please insert name of input .ptcs file to read!")
    exit()

if args.output == None:
    print("Error: Please insert name of output .vtk file to write!")
    exit()

print("Input file: ", inputFileName)
print("Output file: ", outputFileName)

# ------------------------------------------------

with open(inputFileName, "r") as f:
    #DATA = f.read()
    LINES =  f.readlines()

#print(LINES)
#print(type(LINES))

ptcs = LINES[1:]
#print(ptcs)

DATA = np.genfromtxt(ptcs, delimiter=" ")
#print(DATA)
#print(type(DATA))
#print(DATA.shape)

r = DATA[:,:3]
v = DATA[:,3:6]
rho = DATA[:,6]
p = DATA[:,7]
ptype = DATA[:,8]

def create_pointcloud_polydata(points, velocity=None, density=None, pressure=None, ptype=None):
    """https://github.com/lmb-freiburg/demon
    Creates a vtkPolyData object with the point cloud from numpy arrays

    points: numpy.ndarray
        pointcloud with shape (n,3)

    velocity: numpy.ndarray
        uint8 array with velocity for each point. shape is (n,3)

    Returns vtkPolyData object
    """
    vpoints = vtk.vtkPoints()
    vpoints.SetNumberOfPoints(points.shape[0])
    for i in range(points.shape[0]):
        vpoints.SetPoint(i, points[i])
    vpoly = vtk.vtkPolyData()
    vpoly.SetPoints(vpoints)


    if not velocity is None:
        #vvelocity = vtk.vtkUnsignedCharArray()
        vvelocity = vtk.vtkFloatArray()
        vvelocity.SetNumberOfComponents(3)
        vvelocity.SetName("Velocity")
        vvelocity.SetNumberOfTuples(points.shape[0])
        for i in range(points.shape[0]):
            vvelocity.SetTuple3(i ,velocity[i,0],velocity[i,1], velocity[i,2])
        #vpoly.GetPointData().SetVectors(vvelocity)
        vpoly.GetPointData().AddArray(vvelocity)

    if not density is None:
        vdensity = vtk.vtkFloatArray()
        vdensity.SetNumberOfComponents(1)
        vdensity.SetName("Density")
        vdensity.SetNumberOfTuples(points.shape[0])
        for i in range(points.shape[0]):
            vdensity.SetTuple1(i ,density[i])
        vpoly.GetPointData().AddArray(vdensity)

    if not pressure is None:
        vpressure = vtk.vtkFloatArray()
        vpressure.SetNumberOfComponents(1)
        vpressure.SetName("Pressure")
        vpressure.SetNumberOfTuples(points.shape[0])
        for i in range(points.shape[0]):
            vpressure.SetTuple1(i ,pressure[i])
        vpoly.GetPointData().AddArray(vpressure)

    if not ptype is None:
        vptype = vtk.vtkIntArray()
        vptype.SetNumberOfComponents(1)
        vptype.SetName("Ptype")
        vptype.SetNumberOfTuples(points.shape[0])
        for i in range(points.shape[0]):
            vptype.SetTuple1(i ,ptype[i])
        vpoly.GetPointData().AddArray(vptype)

    vcells = vtk.vtkCellArray()

    for i in range(points.shape[0]):
        vcells.InsertNextCell(1)
        vcells.InsertCellPoint(i)

    vpoly.SetVerts(vcells)

    #print(vpoly)
    return vpoly

def save_polydata(polydata, file_name, binary=False, color_array_name=None):
    # get file extension (type)
    file_extension = file_name.split(".")[-1].lower()

    # todo better generic load
    # todo test all
    if file_extension == "vtk":
        writer = vtk.vtkPolyDataWriter()
    elif file_extension == "vtp":
        writer = vtk.vtkPolyDataWriter()
    elif file_extension == "fib":
        writer = vtk.vtkPolyDataWriter()
    elif file_extension == "ply":
        writer = vtk.vtkPLYWriter()
    elif file_extension == "stl":
        writer = vtk.vtkSTLWriter()
    elif file_extension == "xml":
        writer = vtk.vtkXMLPolyDataWriter()
    elif file_extension == "obj":
        raise "mni obj or Wavefront obj ?"
    #    writer = set_input(vtk.vtkMNIObjectWriter(), polydata)

    #print(writer)
    writer.SetFileName(file_name)
    #writer = set_input(writer, polydata)
    writer.SetInputData(polydata)
    if color_array_name is not None:
        writer.SetArrayName(color_array_name);

    if binary :
        writer.SetFileTypeToBinary()
    writer.Update()
    writer.Write()

#switch
myData = create_pointcloud_polydata(r, v, rho, p, ptype)
save_polydata(myData, outputFileName)

import numpy as np
import vtk

#dp = 0.01 #initial particle distance
dp = 0.2 #initial particle distance
outputFileName = 'blockOfParticles.vtk'
r = []; v = []; rho = []; p = []

for x in range(int(0.9/dp)):
    for y in range(int(0.9/dp)):
        r.append([x*dp, y*dp, 0])

np_r = np.array(r)

def create_pointcloud_polydata(points, velocity=None, density=None, pressure=None):
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


#myData = create_pointcloud_polydata(r, v, rho, p)
myData = create_pointcloud_polydata(np_r)
save_polydata(myData, outputFileName)

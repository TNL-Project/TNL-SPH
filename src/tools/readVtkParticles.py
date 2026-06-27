"""
Thin VTK polydata reader for TNL-SPH particle files.

Reads a .vtk polydata file and returns the points plus all PointData arrays
as plain numpy arrays.  This eliminates the repeated 7-line VTK boilerplate
that appears in every case ``init.py``.

Array names are NOT aliased — different pre-processing tools produce
different names (``Vel`` vs ``Velocity``, ``Normal`` vs ``Normals``,
``Areas`` vs ``ElementSize``).  Callers must use the exact array name
their input file provides.
"""

import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa


def read_vtk_particles(filename):
   """
   Read a VTK polydata file and return (points, point_data).

   Parameters
   ----------
   filename : str
      Path to the .vtk polydata file.

   Returns
   -------
   points : ndarray, shape (n, 3), dtype float64
      Particle coordinates.
   point_data : dict[str, ndarray]
      All PointData arrays, keyed by their VTK array name.  Original
      dtypes are preserved (caller may convert with ``astype``).
   """
   reader = vtk.vtkPolyDataReader()
   reader.SetFileName(filename)
   reader.ReadAllScalarsOn()
   reader.ReadAllVectorsOn()
   reader.Update()
   polydata = reader.GetOutput()

   wrapped = dsa.WrapDataObject(polydata)
   points = np.array(wrapped.Points, dtype=float)

   point_data = {}
   vtk_point_data = polydata.GetPointData()
   for i in range(vtk_point_data.GetNumberOfArrays()):
      name = vtk_point_data.GetArrayName(i)
      point_data[name] = np.array(wrapped.PointData[name])

   return points, point_data

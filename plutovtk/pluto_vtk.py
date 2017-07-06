def file_name(n, data_dir):
    """
    Find file name.
    """
    import os
    import os.path

    if (n < 10):
        filename = "{0}/data.000{1}.vtk".format(data_dir, n)
    elif (n > 9 and n < 100):
        filename = "{0}/data.00{1}.vtk".format(data_dir, n)
    elif (n > 99 and n < 1000):
        filename = "{0}/data.0{1}.vtk".format(data_dir, n)
    elif (n > 999 and n < 10000):
        filename = "{0]/data.{1}.vtk".format(data_dir, n)
    return filename

def setup_vtk(filename):
    """
    Set vtk reader and get vector dimansions.
    """
    from vtk import vtkStructuredGridReader
    from vtk import vtkGenericDataObjectReader
    from vtk import vtkRectilinearGridReader
    from vtk import vtkStructuredPointsReader

    reader = vtkGenericDataObjectReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    dim = data.GetDimensions()
    sca = list(dim)
    sca = [i - 1 for i in dim]
    vec = list(dim)
    vec = [i - 1 for i in dim]
    vec.append(3)
    return data, vec, sca, dim

def read_vtk(data, vec, sca):
    """
    Read in the field quantites from the vtk file.
    """
    import numpy as np
    from vtk.util import numpy_support as VN

    rho = VN.vtk_to_numpy(data.GetCellData().GetArray('rho'))
    prs = VN.vtk_to_numpy(data.GetCellData().GetArray('prs'))
    u = VN.vtk_to_numpy(data.GetCellData().GetArray('3D_Velocity_Field'))
    b = VN.vtk_to_numpy(data.GetCellData().GetArray('3D_Magnetic_Field'))

    # Redo if sim geometry is cartesian.
    #rho = rho.reshape(sca, order='F')
    #prs = prs.reshape(sca, order='F')

    #u = u.reshape(vec, order='F')
    #b = b.reshape(vec, order='F')

    return rho, prs, u, b


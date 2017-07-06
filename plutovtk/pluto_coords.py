def cartesian_coordinates_vtk(data, dim):
    """
    Get coordinates from the data object (Cartesian).
    """
    import numpy as np

    xg = np.zeros(data.GetNumberOfPoints())
    yg = np.zeros(data.GetNumberOfPoints())
    zg = np.zeros(data.GetNumberOfPoints())

    for i in range(data.GetNumberOfPoints()):
        xg[i], yg[i], zg[i] = data.GetPoint(i)

    xg = xg.reshape(dim, order='F')
    yg = yg.reshape(dim, order='F')
    zg = zg.reshape(dim, order='F')

    xmin = min(xg[:, 0, 0])
    xmax = max(xg[:, 0, 0])
    ymin = min(yg[0, :, 0])
    ymax = max(yg[0, :, 0])
    zmin = min(zg[0, 0, :])
    zmax = max(zg[0, 0, :])

    x = np.zeros(len(xg[:, 0, 0]) - 1)
    y = np.zeros(len(yg[0, :, 0]) - 1)
    z = np.zeros(len(zg[0, 0, :]) - 1)

    for k in range(len(zg[0, 0, :]) - 1):
        z[k] = (zg[0, 0, k] + zg[0, 0, k + 1]) / 2.0

    for j in range(len(yg[0, :, 0]) - 1):
        y[j] = (yg[0, j, 0] + yg[0, j + 1, 0]) / 2.0

    for i in range(len(xg[:, 0, 0]) - 1):
        x[i] = (xg[i, 0, 0] + xg[i + 1, 0, 0]) / 2.0

    del xg, yg, zg

    a = 0
    grid = np.zeros((len(x) * len(y) * len(z), 3))
    for k in range(len(z)):
        for j in range(len(y)):
            for i in range(len(x)):
                grid[a, 0] = x[i]
                grid[a, 1] = y[j]
                grid[a, 2] = z[k]
                a += 1

    limits = np.array([[xmin, xmax], [ymin, ymax], [zmin, zmax]])
    return x, y, z, grid, limits

def spherical_coordinates_vtk(data, dim):
    """
    Get coordinates from data object (spherical).
    """
    import numpy as np

    theta1 = np.array([i * np.pi / 120.0 for i in range(121)])
    phi1 = np.array([j * np.pi / 120.0 for j in range(241)])

    r1 = np.zeros(data.GetNumberOfPoints())
    r2 = np.zeros(data.GetNumberOfPoints())
    r3 = np.zeros(data.GetNumberOfPoints())

    for i in range(data.GetNumberOfPoints()):
        r1[i], r2[i], r3[i] = data.GetPoint(i)
    r3 = r3.reshape(dim, order='F')

    r = np.zeros(len(r3[:, 0, 0]) - 1)
    theta = np.zeros(len(theta1) - 1)
    phi = np.zeros(len(phi1) - 1)

    for i in range(len(r3[:, 0, 0]) - 1):
        r[i] = (r3[i, 0, 0] + r3[i + 1, 0, 0]) / 2.0

    for j in range(len(theta1) - 1):
        theta[j] = (theta1[j] + theta1[j + 1]) / 2.0

    for k in range(len(phi1) - 1):
        phi[k] = (phi1[k] + phi1[k + 1]) / 2.0

    a = 0
    grid = np.zeros((len(r) * len(theta) * len(phi), 3))
    for k in range(len(phi)):
        for j in range(len(theta)):
            for i in range(len(r)):
                grid[a, 0] = r[i]
                grid[a, 1] = theta[j]
                grid[a, 2] = phi[k]
                a += 1

    rmin = 1.0
    rmax = np.ceil(max(r3[:, 0, 0]))
    thetamin = 0.0
    thetamax = np.pi
    phimin = 0.0
    phimax = 2.0 * np.pi
    limits = np.array([[rmin, rmax], [thetamin, thetamax], [phimin, phimax]])
    del r1,r2,r3

    return r, theta, phi, grid, limits

def convert_coordinates(grid, data_geometry, plot_geometry):
    """
    Convert between Cartesian and spherical polar coordinates.
    """
    import numpy as np

    x1 = grid[:, 0]
    x2 = grid[:, 1]
    x3 = grid[:, 2]
    gridc = np.zeros((len(x1), 3))

    if (data_geometry == plot_geometry):
        print("Error: plot geometry == data geometry")
        return
    elif (plot_geometry == "spherical"):
        gridc[:, 0] = np.sqrt(x1*x1 + x2*x2 + x3*x3)
        gridc[:, 1] = np.arccos(x3 / gridc[:, 0])
        gridc[:, 2] = np.arctan(x2 / x1)
    elif(plot_geometry == "cartesian"):
        gridc[:, 0] = x1 * np.sin(x2) * np.cos(x3)
        gridc[:, 1] = x1 * np.sin(x2) * np.sin(x3)
        gridc[:, 2] = x1 * np.cos(x2)

    del x1,x2,x3
    x1 = np.unique(gridc[:, 0])
    x2 = np.unique(gridc[:, 1])
    x3 = np.unique(gridc[:, 2])
    limits = np.array(
        [[round(x1.min(),-1), round(x1.max(),-1)], 
         [round(x2.min(),-1), round(x2.max(),-1)], 
         [round(x3.min(),-1), round(x3.max(),-1)]])
    return x1, x2, x3, gridc, limits
    #return gridc, limits

def rotation(grid, rot, plot_geometry):
    """
    Rotate the simulation data around the z-axis.
    """
    import numpy as np

    if (plot_geometry == 'cartesian'):
        x1 = grid[:, 0]
        x2 = grid[:, 1]
        x3 = grid[:, 2]
        # rotation about z
        x1z = x1*np.cos(rot) - x2*np.sin(rot)
        x2z = x1*np.sin(rot) + x2*np.cos(rot)
        x3z = x3
        """
        # rotation about y
        x1y = x1z*np.cos(rot[1]) + x3z*np.sin(rot[1])
        x2y = x2z
        x3y = -x1z*np.sin(rot[1]) + x3z*np.cos(rot[1])
        # rotation about x
        x1x = x1y
        x2x = x2y*np.cos(rot[0]) - x3y*np.sin(rot[0])
        x3x = x2y*np.sin(rot[0]) + x3y*np.cos(rot[0])
        """
        grid[:, 0] = x1z
        grid[:, 1] = x2z
        grid[:, 2] = x3z
    elif (plot_geometry == 'spherical'):
        grid[:, 2] = grid[:, 2] + rot[2]
    return grid



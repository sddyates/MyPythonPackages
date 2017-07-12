
# coding: utf-8

import numpy as np

def cartesian_coordinates_vtk(data, dim):
    """
    Synopsis
    --------
    Get coordinates from the data object (Cartesian).

    Parameters
    ----------
    data: object-like
        VTK data object containing the data read from 
        the simulation output.

    dim: int-like
        dimensions of the simulations.     

    Returns
    -------
    x, y, z: np.arrays of the unique x, y, z coordinates 
        of the simulation data.

    grid: 3 dimensional np.array containing all the 
        coordinates which form the grid.

    limits: 3x2 array contaning the limits of the 
        simulation box.

    TODO
    ----
    - Possibly combine x, y, z and grid into just grid.
    """

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
    Synopsis
    --------
    Get coordinates from data object (spherical).

    Parameters
    ----------
    data: object-like
        VTK data object containing the data read from 
        the simulation output.

    dim: int-like
        dimensions of the simulations.     

    Returns
    -------
    r, theta, phi: np.arrays of the unique r, theta, phi 
        coordinates of the simulation data.

    grid: 3 dimensional np.array containing all the 
        coordinates which form the grid.

    limits: 3x2 array contaning the limits of the 
        simulation box.

    TODO
    ----
    - Possibly combine x, y, z and grid into just grid.
    """

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
    Synopsis
    --------
    Convert between Cartesian and spherical polar coordinates.

    Parameters
    ----------
    grid: array-like
        3 dimensional np.array containing all the 
        coordinates which form the grid.

    data_geometry: string-like
        Geometry of the raw data read from the vtk file.     

    plot_geometry: string-like
        Geometry which is needed by the main program.     

    Returns
    -------
    r, theta, phi: np.arrays of the unique r, theta, phi 
        coordinates of the simulation data.

    gridc: array-like
        3 dimensional np.array containing the converted 
        grid.

    limits: 3x2 array contaning the limits of the 
        simulation box (post conversion).

    TODO
    ----
    - Improve error handeling.
    """

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


def rotation(grid, alpha, beta, gamma, plot_geometry):
    """
    Synopsis
    --------
    perform sucsessive rotations via rotation marixies 
    on input data and return a new grid.

    Parameters
    ----------
    grid: array-like
        3 dimensional np.array containing all the 
        coordinates which form the grid.

    alpha: double-like
        Angle to rotate about the x-axis.     

    beta: double-like
        Angle to rotate about the y-axis.     

    gamma: double-like
        Angle to rotate about the z-axis.     

    plot_geometry: string-like
        Geometry of the data.     

    Returns
    -------
    grid_new: array-like
        3 dimensional np.array containing the converted 
        grid.

    TODO
    ----
    - Improve error handling.
    """

    if (plot_geometry == 'cartesian'):
        x1 = grid[:, 0]
        x2 = grid[:, 1]
        x3 = grid[:, 2]
        # rotation about z
        x1z = x1*np.cos(gamma) - x2*np.sin(gamma)
        x2z = x1*np.sin(gamma) + x2*np.cos(gamma)
        x3z = x3
        # rotation about y
        x1y = x1z*np.cos(beta) + x3z*np.sin(beta)
        x2y = x2z
        x3y = -x1z*np.sin(beta) + x3z*np.cos(beta)
        # rotation about x
        x1 = x1y
        x2 = x2y*np.cos(alpha) - x3y*np.sin(alpha)
        x3 = x2y*np.sin(alpha) + x3y*np.cos(alpha)

        grid_new = np.zeros_like(grid)        
        grid_new[:, 0] = x1
        grid_new[:, 1] = x2
        grid_new[:, 2] = x3

    del x1, x2, x3

    return grid_new

if __name__ == '__main__':

    import unittest

    class rotate_test(unittest.TestCase):

        def test_pp(self):
            input_coords = np.array([[1.0, 1.0, 1.0]])
            output_coords = np.array([[1.0, -1.0, 1.0]])

            result = rotation(
                input_coords, 
                np.pi/2.0, 
                np.pi/2.0, 
                np.pi/2.0, 
                'cartesian'
                )

            self.assertTrue(np.allclose(result, output_coords))

    unittest.main(exit=False)







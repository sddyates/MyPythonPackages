
# coding: utf-8

import numpy as np

class Grid_Cartesian:
    """
    Synopsis
    --------
    Get coordinates from the data object (Cartesian).

    Args
    ----
    data: object-like
        VTK data object containing the data read from 
        the simulation output.

    dim: int-like
        dimensions of the simulations.     

    Attributes
    ----------
    x, y, z: np.arrays of the unique x, y, z coordinates 
        of the simulation data.


    TODO
    ----
    - Possibly combine x, y, z and grid into just grid.
    """
    
    def __init__(self, data, dim):
        
        self._xg = np.zeros(data.GetNumberOfPoints())
        self._yg = np.zeros_like(self._yg)
        self._zg = np.zeros_like(self._zg)
        self._populate(self, data)

        # xg, yg, and zg are the coordinates of the cell edges.
        self.z = (self._zg[0, 0, 1:] + self._zg[0, 0, :-1]) / 2.0
        self.y = (self._yg[0, 1:, 0] + self._yg[0, :-1, 0]) / 2.0
        self.x = (self._xg[1:, 0, 0] + self._xg[:-1, 0, 0]) / 2.0

        # Cast x, y and z into a form that can be used 
        # for interpolation.
        self.coords = self._get_coords(self, data)

    def _populate(self, data):
        for i in range(data.GetNumberOfPoints()):
            self._xg[i], self._yg[i], self._zg[i] = data.GetPoint(i)
        
    def _get_coords(self):
        """
        Returns
        -------
        grid: 3 dimensional np.array containing all the 
            coordinates which form the grid.
        """
        return np.flip(np.vstack(np.meshgrid(z, y, x)).reshape(3, -1, order='C').T, 1)


    def limits(self):
        """
        Returns
        -------
        limits: 3x2 array contaning the limits of the 
            simulation box.
        """
        return np.array([[min(self._xg[:, 0, 0]), max(self._xg[:, 0, 0])], 
                         [min(self._yg[0, :, 0]), max(self._yg[0, :, 0])], 
                         [min(self._zg[0, 0, :]), max(self._zg[0, 0, :])]])

class Grid_spherical:
    """
    Synopsis
    --------
    Get coordinates from data object (spherical).

    Args
    ----
    data: object-like
        VTK data object containing the data read from 
        the simulation output.

    dim: int-like
        dimensions of the simulations.     

    Attributes
    ----------
    r, theta, phi: np.arrays of the unique r, theta, phi 
        coordinates of the simulation data.

    TODO
    ----
    - Possibly combine x, y, z and grid into just grid.
    """
    
    def __init__(self, data, dim):
        
        self._r1 = np.zeros(data.GetNumberOfPoints())
        self._r2 = np.zeros_like(self._r1)
        self._r3 = np.zeros_like(self._r1)

        self._theta1 = np.array([i * np.pi / 120.0 for i in range(121)])
        self._phi1 = np.array([j * np.pi / 120.0 for j in range(241)])
        self._populate(data, dim)

        r = np.zeros(len(self._r3[:, 0, 0]) - 1)
        self.r = (self._r3[1:, 0, 0] + self._r3[:-1, 0, 0]) / 2.0
        self.theta = (self._theta1[1:] + self._theta1[:-1]) / 2.0
        self.phi = (self._phi1[1:] + self._phi1[:-1]) / 2.0

        #self.coords = self._get_coords(data)
        self.coords = self._get_coords_old(data)

    def _populate(self, data, dim):
        for i in range(data.GetNumberOfPoints()):
            self._r1[i], self._r2[i], self._r3[i] = data.GetPoint(i)
        self._r3 = self._r3.reshape(dim, order='F')
        
    def _get_coords(self, data):
        """
        Retruns
        -------
        grid: 3 dimensional np.array containing all the 
            coordinates which form the grid.
        """
        return np.flip(np.vstack(
                   np.meshgrid(self.phi, self.theta, self.r)
                   ).reshape(3, -1, order='C').T, 1)

    def _get_coords_old(self, data):
        """
        Retruns
        -------
        grid: 3 dimensional np.array containing all the 
            coordinates which form the grid.
        """
        a = 0
        grid = np.zeros((len(self.r) * len(self.theta) * len(self.phi), 3))
        for k in range(len(self.phi)):
            for j in range(len(self.theta)):
                for i in range(len(self.r)):
                    grid[a, 0] = self.r[i]
                    grid[a, 1] = self.theta[j]
                    grid[a, 2] = self.phi[k]
                    a += 1
        return grid

    def limits(self):
        """
        Returns
        -------
        limits: 3x2 array contaning the limits of the 
            simulation box.
        """
        return np.array([[1, np.ceil(max(self._r3[:, 0, 0]))], 
                         [0.0, np.pi], 
                         [0.0, 2.0*np.pi]])


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



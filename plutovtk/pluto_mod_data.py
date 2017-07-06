def interpolate(x1, x2, x3, grid, x1_rez, x2_rez, x3_rez, U):
    """
    Interpolate the field data on to a regular grid. 
    yt curretly dosn't take rectalinear grids.
    """
    import numpy as np
    import scipy.interpolate as interp

    x1i, x2i, x3i = np.mgrid[min(x1):max(x1):x1_rez, min(
        x2):max(x2):x2_rez, min(x3):max(x2):x3_rez]
    ri = np.sqrt(x1i*x1i + x2i*x2i + x3i*x3i)

    Ui = []
    for i in range(len(U[:])):
        a = interp.griddata(grid, U[i], (x1i, x2i, x3i), method='nearest', fill_value=0.0)
        #a[ri<1.0] = 0.0 
        a[ri>40.0] = 0.0
        Ui.append(a)

    dx1 = (max(x1) - min(x1)) / len(x1i[:, 0, 0])
    dx2 = (max(x2) - min(x2)) / len(x2i[0, :, 0])
    dx3 = (max(x3) - min(x3)) / len(x3i[0, 0, :])
    return Ui, dx1, dx2, dx3



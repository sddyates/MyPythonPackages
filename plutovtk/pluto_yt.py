def yt_setup(Ui, plot_geometry, limits):
    """
    Plot the interpolated data using yt-python.
    """
    import yt
    import numpy as np
    import matplotlib as plt
    import matplotlib.cm as cm
    from yt.utilities.physical_constants import mp, kb
    from yt.units import dimensions


    result = dict(
        density = (Ui[0], "g/cm**3"),
        pressure = (Ui[1], "dyne/cm**2"))#,
    """
        tau_ff1 = (Ui[1], ""),
        tau_ff2 = (Ui[2], ""),
        tau_ff3 = (Ui[3], ""))
    """
        pressure = (Ui[1], "dyne/cm**2"),
        velocity_x = (Ui[2], "cm/s"),
        velocity_y = (Ui[3], "cm/s"),
        velocity_z = (Ui[4], "cm/s"),
        magnetic_field_x = (Ui[5], "G"),
        magnetic_field_y = (Ui[6], "G"),
        magnetic_field_z = (Ui[7], "G"))
        magnetic_field_z = (Ui[7], "G"),
        tau_ff = (Ui[8], ""))

    ds = yt.load_uniform_grid(
        result,
        Ui[0].shape,
        length_unit="9*Rsun",
        mass_unit=2.4578492774399997e+23,
        time_unit=6.26e6,
        velocity_unit="cm/s",
        bbox=limits,
        geometry=plot_geometry)
    
    def _temperature(field, data):
        mu = 1.01
        return (data["gas", "pressure"]*mu*mp)/(data["gas", "density"]*kb)
    ds.add_field(("gas", "temperature"), function=_temperature, units="auto",
        dimensions=dimensions.temperature)

    def _H_number_density(field, data):
        mu_i = 1.01
        return data["gas", "density"]/(mu_i*mp)
    ds.add_field(("gas", "H_number_density"), function=_H_number_density, units="cm**-3")

    return ds



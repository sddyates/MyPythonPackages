# coding: utf-8
import yt
import numpy as np
from yt.fields.api import ValidateParameter
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.utilities.physical_constants import mp, kb
from yt import derived_field
from yt.units.yt_array import YTQuantity
from scipy.spatial.distance import euclidean


def unit_override():
    return {"length_unit":(1,"Rsun"),
            "time_unit":(6.955e+05 ,"s"),
            "mass_unit":(3.36427433875e+17,"g"),
            "magnetic_unit":(1.121e-02,"G")}

def create_fields(ds):

    def _radialvelocity(field, data):
        return data['velocity_x']*data['x']/data['radius'] + \
               data['velocity_y']*data['y']/data['radius'] + \
               data['velocity_z']*data['z']/data['radius']

    def _temperature(field, data):
        return (data["gas", "pressure"]*1.01*mp)/(data["gas", "density"]*kb)

    def _radius_planet(field, data):
        a = 0.047*1.496e+13/6.955e+10
        shift = data.ds.arr(np.ones_like(data['x']))*a
        x_planet = data['x'] - shift
        return np.sqrt(x_planet*x_planet + data['y']*data['y'] + data['z']*data['z'])

    def _ni(field, data):
        return data["density"]/(1.09*mp)

    def _fc(field, data):
        return YTQuantity(2.8, 'G**-1*MHz')*data["gas", "mag_field_strength"]

    def _fp(field, data):
        return YTQuantity(8.98e-3, 'cm**(3/2)*MHz')*np.sqrt(1.01*data["gas", "ni"])

    def _f_ratio(field, data):
        return data['gas', 'fc']/data['gas', 'fp']

    def _BGx1(field, data):

        B0s = YTQuantity(2.0, "G")
        B0p = YTQuantity(1.0, "G")

        Rs = YTQuantity(6.955e+10, "cm") 
        Rp = YTQuantity(1.5*0.10045*Rs, "cm")
        a = YTQuantity(0.047, "au").in_units("cm")

        center = data.get_field_parameter('center')

        x1 = data["x"].in_units('cm')
        x2 = data["y"].in_units('cm')
        x3 = data["z"].in_units('cm')
        rs = np.sqrt(x1*x1 + x2*x2 + x3*x3)
        rp = np.sqrt((x1-a)*(x1-a) + x2*x2 + x3*x3)

        BGx1 = data.ds.arr(np.zeros_like(data["magnetic_field_x"]), "G")

        BGx1 = 3.0*x1*x3*B0s*Rs**3*rs**(-5) + 3.0*(x1 - a)*x3*B0p*Rp**3*rp**(-5)
        BGx1[rs <= Rs] = 3.0*x1[rs <= Rs]*x3[rs <= Rs]*B0s*Rs**3*rs[rs <= Rs]**(-5)
        BGx1[rs <= 0.5*Rs] = 0.0
        BGx1[rp <= Rp] = 3.0*(x1[rp <= Rp] - a)*x3[rp <= Rp]*B0p*Rp**3*rp[rp <= Rp]**(-5)
        BGx1[rp <= 0.5*Rp] = 0.0

        return BGx1

    def _BGx2(field, data):

        B0s = YTQuantity(2.0, "G")
        B0p = YTQuantity(1.0, "G")

        Rs = YTQuantity(6.955e+10, "cm") 
        Rp = YTQuantity(1.5*0.10045*Rs, "cm")
        a = YTQuantity(0.047, "au").in_units("cm")

        center = data.get_field_parameter('center')

        x1 = data["x"].in_units('cm')
        x2 = data["y"].in_units('cm')
        x3 = data["z"].in_units('cm')
        rs = np.sqrt(x1*x1 + x2*x2 + x3*x3)
        rp = np.sqrt((x1-a)*(x1-a) + x2*x2 + x3*x3)

        BGx2 = data.ds.arr(np.zeros_like(data["magnetic_field_y"]), "G")

        BGx2 = 3.0*x3*x2*B0s*Rs**3*rs**(-5) + 3.0*x3*x2*B0p*Rp**3*rp**(-5)
        BGx2[rs <= Rs] = 3.0*x3[rs <= Rs]*x2[rs <= Rs]*B0s*Rs**3*rs[rs <= Rs]**(-5)
        BGx2[rs <= 0.5*Rs] = 0.0
        BGx2[rp <= Rp] = 3.0*x3[rp <= Rp]*x2[rp <= Rp]*B0p*Rp**3*rp[rp <= Rp]**(-5)
        BGx2[rp <= 0.5*Rp] = 0.0

        return BGx2

    def _BGx3(field, data):

        B0s = YTQuantity(2.0, "G")
        B0p = YTQuantity(1.0, "G")

        Rs = YTQuantity(6.955e+10, "cm") 
        Rp = YTQuantity(1.5*0.10045*Rs, "cm")
        a = YTQuantity(0.047, "au").in_units("cm")

        x1 = data["x"].in_units('cm')
        x2 = data["y"].in_units('cm')
        x3 = data["z"].in_units('cm')

        rs = np.sqrt(x1*x1 + x2*x2 + x3*x3)
        rp = np.sqrt((x1-a)*(x1-a) + x2*x2 + x3*x3)

        BG_z = data.ds.arr(np.zeros_like(data["magnetic_field_z"]), "G")

        BGx3 = (3.0*x3*x3 - rs*rs)*B0s*Rs**3*rs**(-5) + (3.0*x3*x3 - rp*rp)*B0p*Rp**3*rp**(-5)
        BGx3[rs <= Rs] = (3.0*x3[rs <= Rs]*x3[rs <= Rs] - \
                rs[rs <= Rs]*rs[rs <= Rs])*B0s*Rs**3*rs[rs <= Rs]**(-5) 
        BGx3[rs <= 0.5*Rs] = 16.0*B0s
        BGx3[rp <= Rp] = (3.0*x3[rp <= Rp]*x3[rp <= Rp] - \
                rp[rp <= Rp]*rp[rp <= Rp])*B0p*Rp**3*rp[rp <= Rp]**(-5)
        BGx3[rp <= 0.5*Rp] = 16.0*B0p

        return BGx3

    def _Bx1(field, data):
        return data["gas", "magnetic_field_x"] + data["gas", "BGx1"]

    def _Bx2(field, data):
        return data["gas", "magnetic_field_y"] + data["gas", "BGx2"] 

    def _Bx3(field, data):
        return data["gas", "magnetic_field_z"] + data["gas", "BGx3"]   

    def _mag_energy(field,data):
        return (data["Bx1"]**2 +
                data["Bx2"]**2 +
                data["Bx3"]**2)/(8*np.pi)

    def _mag_field_strength(field,data):
        return np.sqrt(8.*np.pi*data["mag_energy"])

    def _mag_field_magnitude(field,data):
        return (data["Bx1"]**2 +
                data["Bx2"]**2 +
                data["Bx3"]**2)

    def _plasma_b(field,data):
        return data['pressure']/data['mag_energy']

    units_override = {"length_unit":(1,"Rsun"),
                      "time_unit":(6.955e+05 ,"s"),
                      "mass_unit":(3.36427433875e+17,"g"),
                      "magnetic_unit":(1.121e-02,"G")}

    ds.add_field(('index', "radius_planet"), function=_radius_planet, units="cm", take_log=False)
    ds.add_field(("gas", "BGx1"), function=_BGx1, units="G", take_log=False)
    ds.add_field(("gas", "BGx2"), function=_BGx2, units="G", take_log=False)
    ds.add_field(("gas", "BGx3"), function=_BGx3, units="G", take_log=False)
    ds.add_field(("gas", "Bx1"), function=_Bx1, units="G", take_log=False)
    ds.add_field(("gas", "Bx2"), function=_Bx2, units="G", take_log=False)
    ds.add_field(("gas", "Bx3"), function=_Bx3, units="G", take_log=False)
    ds.add_field(("gas", "mag_energy"), function=_mag_energy, units="g*cm**-1*s**-2", take_log=True)
    ds.add_field(("gas", "mag_field_strength"), function=_mag_field_strength, units="G", take_log=True)
    ds.add_field(("gas", "mag_field_magnitude"), function=_mag_field_magnitude, units="G", take_log=True)
    ds.add_field(("gas", "plasma_b"), function=_plasma_b, units="")
    ds.add_field(("gas", "ni"), function=_ni, units="cm**-3")
    ds.add_field(('gas', "fc"), function=_fc, units="MHz", take_log=False)
    ds.add_field(('gas', "fp"), function=_fp, units="MHz", take_log=False)
    ds.add_field(('gas', "fc/fp"), function=_f_ratio, units="", take_log=True)
    ds.add_field(('gas', "temperature"), function=_temperature, units="K", take_log=True)
    ds.add_field(('gas', "radialvelocity"), function=_radialvelocity, units="cm/s", take_log=False)



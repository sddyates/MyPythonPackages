# coding: utf-8
import yt
import numpy as np
from yt.fields.api import ValidateParameter
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.utilities.physical_constants import mp, kb
from yt import derived_field
from yt.units.yt_array import YTQuantity
from yt.funcs import just_one
from scipy.spatial.distance import euclidean
from yt.fields.derived_field import \
    ValidateSpatial

def unit_override():
    return {"length_unit":(1,"Rsun"),
            "time_unit":(6.955e+05 ,"s"),
            "mass_unit":(3.36427433875e+17,"g"),
            "magnetic_unit":(1.121e-02,"G")}


def sim_parameters():
    return {"gamma":1.05}


def create_fields(ds):

    units_override = {"length_unit":(1,"Rsun"),
                      "time_unit":(6.955e+05 ,"s"),
                      "mass_unit":(3.36427433875e+17,"g"),
                      "magnetic_unit":(1.121e-02,"G")}

    def _radialvelocity(field, data):
        return data['velocity_x']*data['x']/data['radius'] + \
               data['velocity_y']*data['y']/data['radius'] + \
               data['velocity_z']*data['z']/data['radius']

    ds.add_field(('gas', "radialvelocity"), 
                 function=_radialvelocity, 
                 units="cm/s", 
                 take_log=False)


    def _sound_speed(field, data):
        gamma = 1.05
        ftype = field.name[0]
        tr = gamma * data[ftype, "pressure"] / data[ftype, "density"]
        return np.sqrt(tr)

    ds.add_field(('gas', "sound_speed"), 
                 function=_sound_speed, 
                 units="cm/s", 
                 take_log=False)


    def _mach_number(field, data):
        """ M{|v|/c_sound} """
        ftype = field.name[0]
        return data[ftype, "velocity_magnitude"] / data[ftype, "sound_speed"]

    ds.add_field(('gas', "mach_number"), 
                 function=_mach_number, 
                 units="", 
                 take_log=False)


    def _temperature(field, data):
        return (data["gas", "pressure"]*mp)/(2.0*data["gas", "density"]*kb)

    ds.add_field(('gas', "temperature"), 
                 function=_temperature, 
                 units="K", 
                 take_log=True)


    def _radius_planet(field, data):
        a = 0.047*1.496e+13/6.955e+10
        shift = data.ds.arr(np.ones_like(data['x']))*a
        x_planet = data['x'] - shift
        return np.sqrt(x_planet*x_planet \
                       + data['y']*data['y'] \
                       + data['z']*data['z'])
    ds.add_field(('index', "radius_planet"), 
                 function=_radius_planet, 
                 units="cm", 
                 take_log=False)


    def _ni(field, data):
        return data["density"]/(1.09*mp)

    ds.add_field(("gas", "ni"), 
                 function=_ni, 
                 units="cm**-3")


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
        BGx1[rp <= Rp] = 3.0*(x1[rp <= Rp] - a)*x3[rp <= Rp]\
                         *B0p*Rp**3*rp[rp <= Rp]**(-5)
        BGx1[rp <= 0.5*Rp] = 0.0

        return BGx1

    ds.add_field(("gas", "BGx1"), 
                 function=_BGx1, 
                 units="G", 
                 take_log=False)



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
        BGx2[rs <= Rs] = 3.0*x3[rs <= Rs]*x2[rs <= Rs]\
                         *B0s*Rs**3*rs[rs <= Rs]**(-5)
        BGx2[rs <= 0.5*Rs] = 0.0
        BGx2[rp <= Rp] = 3.0*x3[rp <= Rp]*x2[rp <= Rp]\
                         *B0p*Rp**3*rp[rp <= Rp]**(-5)
        BGx2[rp <= 0.5*Rp] = 0.0

        return BGx2

    ds.add_field(("gas", "BGx2"), 
                 function=_BGx2, 
                 units="G", 
                 take_log=False)


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

        BGx3 = (3.0*x3*x3 - rs*rs)*B0s*Rs**3*rs**(-5) \
               + (3.0*x3*x3 - rp*rp)*B0p*Rp**3*rp**(-5)
        BGx3[rs <= Rs] = (3.0*x3[rs <= Rs]*x3[rs <= Rs] - \
                rs[rs <= Rs]*rs[rs <= Rs])*B0s*Rs**3*rs[rs <= Rs]**(-5) 
        BGx3[rs <= 0.5*Rs] = 16.0*B0s
        BGx3[rp <= Rp] = (3.0*x3[rp <= Rp]*x3[rp <= Rp] - \
                rp[rp <= Rp]*rp[rp <= Rp])*B0p*Rp**3*rp[rp <= Rp]**(-5)
        BGx3[rp <= 0.5*Rp] = 16.0*B0p

        return BGx3

    ds.add_field(("gas", "BGx3"), 
                 function=_BGx3, 
                 units="G", 
                 take_log=False)


    def _Bx1(field, data):
        return data["gas", "magnetic_field_x"] + data["gas", "BGx1"]

    ds.add_field(("gas", "Bx1"), 
                 function=_Bx1, 
                 units="G", 
                 take_log=False)


    def _Bx2(field, data):
        return data["gas", "magnetic_field_y"] + data["gas", "BGx2"] 

    ds.add_field(("gas", "Bx2"), 
                 function=_Bx2, 
                 units="G", 
                 take_log=False)


    def _Bx3(field, data):
        return data["gas", "magnetic_field_z"] + data["gas", "BGx3"]   

    ds.add_field(("gas", "Bx3"), 
                 function=_Bx3, 
                 units="G", 
                 take_log=False)


    def _mag_energy(field,data):
        return (data["Bx1"]**2 +
                data["Bx2"]**2 +
                data["Bx3"]**2)/(8*np.pi)

    ds.add_field(("gas", "mag_energy"), 
                 function=_mag_energy, 
                 units="g*cm**-1*s**-2", 
                 take_log=True)


    def _mag_field_strength(field,data):
        return np.sqrt(8.*np.pi*data["mag_energy"])

    ds.add_field(("gas", "mag_field_strength"), 
                 function=_mag_field_strength, 
                 units="G", 
                 take_log=True)


    def _mag_field_magnitude(field,data):
        return np.sqrt(data["Bx1"]**2 +
                data["Bx2"]**2 +
                data["Bx3"]**2)

    ds.add_field(("gas", "mag_field_magnitude"), 
                 function=_mag_field_magnitude, 
                 units="G", 
                 take_log=True)


    def _plasma_b(field,data):
        return data['pressure']/data['mag_energy']

    ds.add_field(("gas", "plasma_b"), 
                 function=_plasma_b, 
                 units="")


    def _B_divergence(field, data):
        sl_right = slice(None, -2, None)
        sl_left = slice(2, None, None) 
        div_fac = 2.0
        ds = div_fac*just_one(data["index", "dx"])
        f = data["Bx1"][sl_right, 1:-1, 1:-1]/ds
        f -= data["Bx1"][sl_left, 1:-1, 1:-1]/ds
        ds = div_fac * just_one(data["index", "dy"])
        f += data["Bx2"][1:-1, sl_right, 1:-1]/ds
        f -= data["Bx2"][1:-1, sl_left, 1:-1]/ds
        ds = div_fac * just_one(data["index", "dz"])
        f += data["Bx3"][1:-1, 1:-1, sl_right]/ds
        f -= data["Bx3"][1:-1, 1:-1, sl_left ]/ds
        new_field = data.ds.arr(np.zeros(data["Bx1"].shape, dtype=np.float64),
                                f.units)        
        new_field[1:-1, 1:-1, 1:-1] = f
        return np.abs(new_field)

    ds.add_field(("gas", "B_divergence"), 
                 function=_B_divergence, 
                 units="G/code_length", 
                 validators=[ValidateSpatial(1)], 
                 take_log=True)


    def _divB_measure(field, data):
        return data["index", "dx"]*np.abs(data['B_divergence'])\
               /data['mag_field_magnitude']

    ds.add_field(("gas", "divB_measure"), 
                 function=_divB_measure, 
                 units="dimensionless", 
                 take_log=True)


    def _mag_alfven_speed(field,data):
        ftype = field.name[0]
        B = data[ftype,'mag_field_strength']
        return B/np.sqrt(4.0*np.pi*data[ftype,'density'])

    ds.add_field(("gas", "mag_alfven_speed"), 
                 function=_mag_alfven_speed, 
                 units="cm/s", 
                 take_log=True)


    def _mag_mach_alfven(field,data):
        ftype = field.name[0]
        return data[ftype,'velocity_magnitude']/data[ftype,'mag_alfven_speed']

    ds.add_field(("gas", "mag_mach_alfven"), 
                 function=_mag_mach_alfven, 
                 units="", 
                 take_log=True)


    def _fc(field, data):
        return YTQuantity(2.8, 'G**-1*MHz')*data["gas", "mag_field_strength"]

    ds.add_field(('gas', "fc"), 
                 function=_fc, 
                 units="MHz", 
                 take_log=True)


    def _fp(field, data):
        return YTQuantity(8.98e-3, 'cm**(3/2)*MHz')*\
               np.sqrt(1.01*data["gas", "ni"])

    ds.add_field(('gas', "fp"), 
                 function=_fp, 
                 units="MHz", 
                 take_log=True)


    def _f_ratio(field, data):
        return data['gas', 'fc']/data['gas', 'fp']

    ds.add_field(('gas', "fc/fp"), 
                 function=_f_ratio, 
                 units="", 
                 take_log=True)


    def _specific_angular_momentum_density_x(field, data):
        return data["gas", "specific_angular_momentum_x"]*data["gas", "density"]

    ds.add_field(('index', "specific_angular_momentum_density_x"), 
                 function=_specific_angular_momentum_density_x, 
                 units="cm**-1*g*s**-1", 
                 take_log=False,
                 force_override=True)


    def _specific_angular_momentum_density_y(field, data):
        return data["gas", "specific_angular_momentum_y"]*data["gas", "density"]

    ds.add_field(('index', "specific_angular_momentum_density_y"), 
                 function=_specific_angular_momentum_density_y, 
                 units="cm**-1*g*s**-1", 
                 take_log=False,
                 force_override=True)


    def _specific_angular_momentum_density_z(field, data):
        return data["gas", "specific_angular_momentum_z"]*data["gas", "density"]

    ds.add_field(('index', "specific_angular_momentum_density_z"), 
                 function=_specific_angular_momentum_density_z, 
                 units="cm**-1*g*s**-1", 
                 take_log=False,
                 force_override=True)

    ds.periodicity = (True, True, True)


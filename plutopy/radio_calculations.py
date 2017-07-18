import numpy as np
import scipy.integrate as integ
import scipy.interpolate as interp
import matplotlib.pyplot as plt

def radio_emission_temp(rhoi, prsi, nu, Z, D, gamma, dx1, dx2, dx3, R_star, ri, iso):
    """
    Calculate the radio emission from the sumlation results
    in mJy.
    """

    h = 6.62606957e-27 # [erg s]
    c = 2.99e+10 # [cm s-1]
    k_B = 1.3806488e-16 # [erg K-1]
    kpc = 3.08567758e+21 # [cm]
    D = D * kpc # [cm]
    m_H = 1.67262178e-24 # [g]
    mu_i = 1.01 # [Dimensionless]
    mJy = 1E-26 # [egs s-1 cm-2 Hz-1]
    R_sun = 6.963e+10 # [cm]
    R_star *= R_sun # [cm]
    dx1 *= R_star # [cm]
    dx2 *= R_star
    dx3 *= R_star

    rhoi[ri < 1.0] = 1.0
    T = np.longdouble(prsi[:,:,:])*mu_i*m_H/(np.longdouble(rhoi[:,:,:])*k_B)

    if iso:
        T[:,:,:] = 36.3e+3
        T[ri > 40.0] = 0.0
    else:
        T[ri < 1.0] = 36.3e+3
        T[ri > 40.0] = 0.0
        T[np.isnan(T)] = 0.0

    B_nu = 2.0*k_B*T*nu**2/c**2
    B_nu[ri > 40.0] = 0.0
    B_nu[np.isnan(B_nu)] = 0.0

    g_ff = 9.77 + 1.27 * np.log10((T**1.5) / (nu * Z))
    g_ff[ri > 40.0] = 0.0
    g_ff[np.isnan(g_ff)] = 0.0

    K = 0.0178 * ((Z**2 * g_ff) / ((T**1.5) * nu**2))
    K[ri > 40.0] = 0.0
    K[np.isnan(K)] = 0.0

    n_i2 = np.power(np.longdouble(rhoi) / (mu_i * m_H), 2)
    n_i2[ri > 40.0] = 0.0
    n_i2[np.isnan(n_i2)] = 0.0

    A = K*n_i2
    A[ri > 40.0] = 0.0
    A[np.isnan(A)] = 0.0

    tau = gamma*integ.cumtrapz(A, dx=dx1, axis=0, initial=0.0)
    tau[ri > 40.0] = 0.0
    tau[np.isnan(tau)] = 0.0

    tmax = np.amax(tau, axis=0)

    x, y = np.indices(tau.argmax(axis=0).shape)
    B_nu_max = B_nu[tau.argmax(axis=0), x, y]

    '''
    var_strings = ['T', 'K', 'g_ff', 'B_nu', 'n_i2', 'n_i', 'tau', 'tmax', 'B_nu_max']
    list_of_vars = [T, K, g_ff, B_nu, n_i2, np.sqrt(n_i2), tau, tmax, B_nu_max]
    for i, var in enumerate(list_of_vars):
        if len(np.shape(var)) > 2: 
            slice_point = int(len(var[:, 0, 0])/2.0)
            fig = plt.figure()
            plt.imshow(np.float64(np.log10(var[:, slice_point, :])), cmap=plt.cm.Reds, interpolation='none', extent=[-40.0, 40.0, -40.0, 40.0], origin="lower")
            plt.colorbar()
            plt.savefig('frames/test/{0}.png'.format(var_strings[i]))
        else:
            pass
            fig = plt.figure()
            plt.imshow(np.float64(np.log10(var)), cmap=plt.cm.Reds, interpolation='none', extent=[-40.0, 40.0, -40.0, 40.0])
            plt.colorbar()
            plt.savefig('frames/test/{0}.png'.format(var_strings[i]))
    '''       
     
    I_nu = B_nu_max*(1.0 - np.exp(-tmax))
    B = integ.simps(I_nu, dx=dx2, axis=0)
    S_nu = integ.simps(B, dx=dx3, axis=0)/(D**2 * mJy)

    #I_nu = integ.simps(B_nu*np.exp(-tau)*K*n_i2, dx=dx1, axis=0)
    #S_nu = integ.simps(integ.simps(I_nu, dx=dx2, axis=0), dx=dx3, axis=0)/(D**2 * mJy)

    return S_nu, I_nu

def radio_emission_temp_old(rhoi, prsi, nu, Z, D, gamma, dx1, dx2, dx3, R_star, ri):
     """
     Calculate the radio emission from the sumlation results
     in mJy.
     """
     import numpy as np
     import scipy.integrate as integ
     import scipy.interpolate as interp
 
     h = 6.62606957e-27 # [erg s]
     c = 2.99e+10 # [cm s-1]
     k_B = 1.3806488e-16 # [erg K-1]
     kpc = 3.08567758e+21 # [cm]
     D = D * kpc # [cm]
     m_H = 1.67262178e-24 # [g]
     mu_i = 1.01 # [Dimensionless]
     mJy = 1E-26 # [egs s-1 cm-2 Hz-1]
     R_sun = 6.963e+10 # [cm]
     R_star *= R_sun # [cm]
     dx1 *= R_star # [cm]
     dx2 *= R_star
     dx3 *= R_star
 
     rhoi[ri < 1.0] = 1.0
     T = np.longdouble(prsi)*mu_i*m_H/(np.longdouble(rhoi)*k_B)
     T[ri < 1.0] = 36.0e+3
     T[ri > 40.0] = 300.0
 
     B_nu = 2.0*k_B*T*nu**2/c**2
     B_nu[np.isnan(B_nu)] = 0.0
     B_nu[ri > 40.0] = 0.0
     
     g_ff = 9.77 + 1.27 * np.log10((T**1.5) / (nu * Z))
     g_ff[np.isnan(g_ff)] = 0.0
     g_ff[ri > 40.0] = 0.0
     
     K = 0.0178 * ((Z**2 * g_ff) / ((T**1.5) * nu**2))
     K[np.isnan(K)] = 0.0
     K[ri > 40.0] = 0.0
     
     n_i2 = np.power(np.longdouble(rhoi) / (mu_i * m_H), 2)
     n_i2[np.isnan(n_i2)] = 0.0
     n_i2[ri > 40.0] = 0.0
     
     tau = gamma*integ.cumtrapz(K*n_i2, dx=dx1, axis=0, initial=0.0)
     tau[np.isnan(tau)] = 0.0
     
     I_nu = integ.simps(B_nu*np.exp(-tau)*K*gamma*n_i2, dx=dx1, axis=0)
     
     S_nu = integ.simps(integ.simps(I_nu, dx=dx2, axis=0), dx=dx3, axis=0)/(D**2 * mJy)
     
     return S_nu, I_nu

def radio_emission(rhoi, nu, Z, D, gamma, dx1, dx2, dx3, R_star, ri):
    """
    Calculate the radio emission from the sumlation results
    in mJy.
    """

    h = 6.62606957e-27 # [erg s]
    c = 2.99e+10 # [cm s-1]
    k_B = 1.3806488e-16 # [erg K-1]
    kpc = 3.08567758e+21 # [cm]
    D = D * kpc # [cm]
    m_H = 1.67262178e-24 # [g]
    mu_i = 1.01 # [Dimensionless]
    mJy = 1E-26 # [egs s-1 cm-2 Hz-1]
    R_sun = 6.963e+10 # [cm]
    R_star *= R_sun # [cm]
    dx1 *= R_star # [cm]
    dx2 *= R_star
    dx3 *= R_star

    rhoi[ri < 1.0] = 1.0e10

    T = 36.0e+3
    B_nu = 2.0*k_B*T*nu**2/c**2
    g_ff = 9.77 + 1.27 * np.log10((T**1.5) / (nu * Z))
    K = 0.0178 * ((Z**2 * g_ff) / ((T**1.5) * nu**2))

    n_i2 = np.power(np.longdouble(rhoi) / (mu_i * m_H), 2)
    tmax = gamma * K * integ.simps(n_i2, dx=dx1, axis=0)
    Int = B_nu*(1.0 - np.exp(-tmax))
    B = integ.simps(Int, dx=dx2, axis=0)
    S_nu = integ.simps(B, dx=dx3, axis=0)/(D**2 * mJy)

    return S_nu, Int

def tau(rhoi, prsi, nu, T, Z, D, gamma, dx1, R_star, ri):
    """
    This function calculates the free-free optical depth 
    from the density grid provided. The function assumes 
    the observer is situated at x = -\inf.
    """

    k_B = 1.3806488e-16 # [erg K-1]
    m_H = 1.67262178e-24 # [g]
    mu_i = 1.01 # [Dimensionless]
    R_sun = 6.963e+10 # [cm]
    R_star *= R_sun # [cm]
    dx1 *= R_star # [cm]

    rhoi[ri < 1.0] = 1.0e10

    Temp = prsi*mu_i*m_H/(rhoi*k_B)
    Temp[ri < 1.0] = T
    Temp[ri > 40.0] = 300.0

    g_ff = 9.77 + 1.27 * np.log10((Temp**1.5) / (nu * Z))
    K = 0.0178 * ((Z**2 * g_ff) / ((Temp**1.5) * nu**2))
    n_i2 = (rhoi / (mu_i * m_H))**2

    tau_ff = np.zeros_like(rhoi)
    tau_ff = gamma*integ.cumtrapz(K*n_i2, dx=dx1, axis=0, initial=0.0)

    return tau_ff


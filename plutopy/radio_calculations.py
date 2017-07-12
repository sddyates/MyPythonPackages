
def radio_emission_temp(rhoi, prsi, nu, Z, D, gamma, dx1, dx2, dx3, R_star, ri):
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

    T = prsi*mu_i*m_H/(rhoi*k_B)
    T[ri < 1.0] = 36.0e+3
    T[ri > 40.0] = 300.0

    rhoi[ri < 1.0] = 1.0e10

    B_nu = 2.0*k_B*T*nu**2/c**2
    g_ff = 9.77 + 1.27 * np.log10((T**1.5) / (nu * Z))
    K = 0.0178 * ((Z**2 * g_ff) / ((T**1.5) * nu**2))

    n_i2 = np.power(np.longdouble(rhoi) / (mu_i * m_H), 2)

    tau = gamma*integ.cumtrapz(K*n_i2, dx=dx1, axis=0, initial=0.0)

    I_nu = integ.simps(B_nu*np.exp(-tau)*K*n_i2, dx=dx1, axis=0)

    S_nu = integ.simps(integ.simps(I_nu, dx=dx2, axis=0), dx=dx3, axis=0)/(D**2 * mJy)

    return S_nu, I_nu

def radio_emission(rhoi, nu, Z, D, gamma, dx1, dx2, dx3, R_star, ri):
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

    import numpy as np
    import scipy.integrate as integ

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


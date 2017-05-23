def radio_emission(rhoi, prsi, nu, Z, D, gamma, dx1, dx2, dx3, R_star, ri):
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
    #T = 36.0e+3

    B_nu = 2.0*k_B*T*nu**2/c**2
    #B_nu[np.isnan(B_nu)] = 0.0

    g_ff = 9.77 + 1.27 * np.log10((T**1.5) / (nu * Z))
    #g_ff[np.isnan(g_ff)] = 0.0

    K = 0.0178 * ((Z**2 * g_ff) / ((T**1.5) * nu**2))
    #K[np.isnan(K)] = 0.0
    #K[np.isinf(K)] = 0.0
    #print(K)

    n_i2 = (rhoi / (mu_i * m_H))**2

    """
    tmax = gamma * integ.simps(K*n_i2, dx=dx1, axis=0)
    B = 1.0 - np.exp(-tmax)
    C = integ.simps(B_nu*B, dx=dx2, axis=0)
    S_nu = integ.simps(C, dx=dx3, axis=0)/(D**2 * mJy)
    print(S_nu)
    """
    
    #interp.interp1d()
    x = dx1*np.linspace(-39.4562817, 39.4562817, len(n_i2[:,0,0]))
    print(x, len(x))

    tau_ff = gamma*integ.cumtrapz(K*n_i2, dx=dx1, axis=0, initial=0.0)
    tau_ff[np.isnan(tau_ff)] = 0.0
    #print(tau_ff)
    I_nu = integ.simps(B_nu*np.exp(-tau_ff)*K*n_i2, dx=dx1, axis=0)
    I_nu[np.isnan(I_nu)] = 0.0
    S_nu = integ.simps(integ.simps(I_nu, dx=dx2, axis=0), dx=dx3, axis=0)/(D**2 * mJy)
    print(S_nu)
    
    return S_nu

def tau(rhoi, nu, T, Z, D, gamma, dx1, R_star):
    """
    This function calculates the free-free optical depth 
    from the density grid provided. The function assumes 
    the observer is situated at x = -\inf.
    """

    import numpy as np
    import scipy.integrate as integ

    m_H = 1.67262178e-24 # [g]
    mu_i = 1.01 # [Dimensionless]
    R_sun = 6.963e+10 # [cm]
    R_star *= R_sun # [cm]
    dx1 *= R_star # [cm]

    g_ff = 9.77 + 1.27 * np.log10((T**1.5) / (nu * Z))
    K = 0.0178 * ((Z**2 * g_ff) / ((T**1.5) * nu**2))
    n_i2 = (rhoi / (mu_i * m_H))**2

    tau_ff = np.zeros_like(rhoi)
    tau_ff = K*gamma*integ.cumtrapz(n_i2, dx=dx1, axis=0, initial=0.0)

    return tau_ff

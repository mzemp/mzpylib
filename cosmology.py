##############################
#
# Some cosmological functions
#
# Written by Marcel Zemp
#
##############################

import numpy
import scipy

# Some constants

# Exact definitions

c_light = 299792458. # m s^{-1}
au = 149597870700. # m
yr = 31557600. # s (Julian year = 365.25 d)

# Fundamental measured constants: central values

G_Newton = 6.67384e-11 # m^3 s^{-2} kg^{-1} (NIST)
k_Boltzmann = 1.3806488e-23 # J K^{-1} (NIST)
uamu = 1.660538921e-27 # kg (NIST)
m_proton = 1.672621777e-27 # kg (NIST)
m_neutron =  1.674927351e-27 # kg (NIST)
m_electron = 9.10938291e-31 # kg (NIST)
eV = 1.602176565e-19 # J (NIST) 
mu_sun = 1.32712440018e20 # m^3 s^{-2} (JPL)
Mo = 1.98855e30 # kg (Wikipedia)
Lo = 3.845e26 # W (Binney & Tremaine)

# Fundamental measured constants: standard uncertainty

sigma_G_Newton = 8e-15 # m^3 s^{-2} kg^{-1} (NIST)
sigma_k_Boltzmann = 1.3e-29 # J K^{-1} (NIST)
sigma_uamu = 7.3e-35 # kg (NIST)
sigma_m_proton = 7.4e-35 # kg (NIST)
sigma_m_neutron = 7.4e-35 # kg (NIST)
sigma_m_electron = 4e-38 # kg (NIST)
sigma_eV = 3.5e-27 # J (NIST)
sigma_mu_sun = 8e9 # m^3 s^{-2} (JPL)
sigma_Mo = 2.5e26 # kg (Wikipedia)
sigma_Lo = 8e23 # W (Binney & Tremaine)

# Derived constants: central values

pc = au/numpy.tan(numpy.pi/(180*3600)) # m
G_Newton_cosmology = mu_sun*pow(1e9*yr,2)/pow(1e3*pc,3) # kpc^3 Gyr^{-2} Mo^{-1}
rho_crit = (3*pow(1e-1/pc,2))/(8*numpy.pi*G_Newton) # h^2 kg m^{-3}
rho_crit_cosmology = (3*pow(1e8*yr/pc,2))/(8*numpy.pi*G_Newton_cosmology) # h^2 Mo kpc^{-3}
km_per_s_2_kpc_per_Gyr = 1e9*yr/pc
kpc_per_Gyr_2_km_per_s = pc/(1e9*yr)

# Derived constants: standard uncertainty

sigma_pc = 0 # m
sigma_G_Newton_cosmology = sigma_mu_sun*pow(1e9*yr,2)/pow(1e3*pc,3) # kpc^3 Gyr^{-2} Mo^{-1}
sigma_rho_crit = rho_crit*sigma_G_Newton/G_Newton # h^2 kg m^{-3}
sigma_rho_crit_cosmology = rho_crit_cosmology*sigma_G_Newton_cosmology/G_Newton_cosmology # h^2 Mo kpc^{-3}
sigma_km_per_s_2_kpc_per_Gyr = 0
sigma_kpc_per_Gyr_2_km_per_s = 0

# Functions for setting and printing cosmologial parameters

def set_cosmological_parameters(H0=70,OmegaB0=0.05,OmegaD0=0.25,OmegaL0=0.7,OmegaK0=0.0,OmegaR0=0.0,Delta_DC0_NG=0.0):

    cp = {}

    # Fundamental parameters

    cp['H0'] = H0
    cp['OmegaB0'] = OmegaB0
    cp['OmegaD0'] = OmegaD0
    cp['OmegaL0'] = OmegaL0
    cp['OmegaK0'] = OmegaK0
    cp['OmegaR0'] = OmegaR0
    cp['Delta_DC0_NG'] = Delta_DC0_NG

    # Derived parameters

    cp['OmegaM0'] = OmegaD0+OmegaB0
    cp['h0_100'] = H0/100.0
    cp['rho_crit0'] = rho_crit_cosmology*pow(H0/100.0,2)
    cp['rho_mean0'] = cp['OmegaM0']*cp['rho_crit0']

    return cp

def print_cosmological_parameters(cp):

    print
    print 'Cosmological parameters today at z=0 / a=1:'
    print
    print 'H0 = %g km s^-1 Mpc^-1'%(cp['H0'])
    print 'OmegaB0 = %g'%(cp['OmegaB0'])
    print 'OmegaD0 = %g'%(cp['OmegaD0'])
    print 'OmegaM0 = %g'%(cp['OmegaM0'])
    print 'OmegaL0 = %g'%(cp['OmegaL0'])
    print 'OmegaK0 = %g'%(cp['OmegaK0'])
    print 'OmegaR0 = %g'%(cp['OmegaR0'])
    print 'Delta_DC0_NG = %g'%(cp['Delta_DC0_NG'])
    print 'h0_100 = %g'%(cp['h0_100'])
    print 'rho_crit0 = %g Mo kpc^-3'%(cp['rho_crit0'])
    print 'rho_mean0 = %g Mo kpc^-3'%(cp['rho_mean0'])
    print

# Cosmologial E(a) function

def Ecosmo(a,cp):

    return numpy.sqrt(cp['OmegaR0']*pow(a,-4) + cp['OmegaM0']*pow(a,-3) + cp['OmegaK0']*pow(a,-2) + cp['OmegaL0'])

# Cosmic time

def tcosmic(a,cp):

    integrand = lambda x: 1/(x*Ecosmo(x,cp))
    return scipy.integrate.quad(integrand,0,a)[0]*1e3/cp['H0']/km_per_s_2_kpc_per_Gyr

# Scale factor

def ascale(t,cp,a0=1):

    f = lambda x: tcosmic(x,cp)-t
    return scipy.optimize.newton(f,a0)

# Growth factor

def D(a,cp):

    integrand = lambda x: 1/pow(x*Ecosmo(x,cp),3)
    return 2.5*cp['OmegaM0']*Ecosmo(a,cp)*scipy.integrate.quad(integrand,0,a)[0]

# DC mode

def Delta_DC(a,cp):

    return D(a,cp)*cp['Delta_DC0_NG']

# abox in Eulerian view

def abox_EV(a,cp):

    return a*pow(1+Delta_DC(a,cp),-1/3.0)

# alpha-beta-gamma density profiles with 1. & 2. derivative

def rho_abc(r,rs=1,rho0=1,alpha=1,beta=3,gamma=1,rcutoff=numpy.inf,fdecay=0.3):

    rs = float(rs)
    rho0 = float(rho0)
    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)
    rcutoff = float(rcutoff)
    rdecay = float(fdecay)*rcutoff
    delta = rcutoff/rdecay - (gamma+beta*pow(rcutoff/rs,alpha))/(1+pow(rcutoff/rs,alpha))

    def tau(r):
        return pow(r/rs,gamma)*pow(1+pow(r/rs,alpha),(beta-gamma)/alpha)

    def eta(r):
        return gamma/r + ((beta-gamma)/rs)*(pow(r/rs,alpha-1)/(1+pow(r/rs,alpha)))

    def detadr(r):
        d1 = (alpha-1)*(1+pow(r/rs,alpha))-alpha*pow(r/rs,alpha)
        d2 = pow(1+pow(r/rs,alpha),2)
        return -gamma/pow(r,2) + ((beta-gamma)/pow(rs,2))*pow(r/rs,alpha-2)*d1/d2

    if (beta > 3 or r[-1] <= rcutoff):
        rho = rho0/tau(r)
        drhodr = -rho*eta(r)
        d2rhodr2 = rho*(pow(eta(r),2)-detadr(r))
	return rho, drhodr, d2rhodr2
    else:
        si = -1
        for i in range(len(r)):
            if (r[i] >= rcutoff):
                si = i
                break
        assert(si > -1)
        rho = rho0/tau(r[:si])
        drhodr = -rho*eta(r[:si])
        d2rhodr2 = rho*(pow(eta(r[:si]),2)-detadr(r[:si]))
        rhocutoff = rho0/tau(rcutoff)*pow(r[si:]/rcutoff,delta)*numpy.exp(-(r[si:]-rcutoff)/rdecay)
        drhocutoffdr = rhocutoff*(delta/r[si:]-1/rdecay)
        d2rhocutoffdr2 = rhocutoff*(pow(delta/r[si:]-1/rdecay,2)-delta/pow(r[si:],2))
        return numpy.append(rho,rhocutoff), numpy.append(drhodr,drhocutoffdr), numpy.append(d2rhodr2,d2rhocutoffdr2)

# Einasto density profile with 1. & 2. derivative

def rho_einasto(r,rs=1,rho0=1,alpha=0.18):

    # alpha=0.18 from Springel et al. http://adsabs.harvard.edu/abs/2008MNRAS.391.1685S

    rs = float(rs)
    rho0 = float(rho0)
    alpha = float(alpha)

    rho = rho0*numpy.exp(-(2/alpha)*(pow(r/rs,alpha)-1))
    drhodr = rho*(-2*pow(r/rs,alpha-1)/rs)
    d2rhodr2 = drhodr*(-2*pow(r/rs,alpha-1)/rs)+rho*(-2*(alpha-1)*pow(r/rs,alpha-2)/pow(rs,2))

    return rho, drhodr, d2rhodr2

# Helper function to calculate logarithmic derivatives

def rho_logarithmic_derivatives(r,rho,drhodr,d2rhodr2):

    assert(len(r) == len(rho))
    assert(len(r) == len(drhodr))
    assert(len(r) == len(d2rhodr2))

    dlogrhodlogr = r*drhodr/rho
    d2logrhodlogr2 = dlogrhodlogr-pow(r*drhodr/rho,2)+pow(r,2)*d2rhodr2/rho

    return dlogrhodlogr, d2logrhodlogr2

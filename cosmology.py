##############################
#
# Some cosmological functions
#
# Written by Marcel Zemp
#
##############################

from numpy import sqrt
from scipy import integrate

# Some constants

km_per_s_2_kpc_per_Gyr = 1.02271216511
G = 4.4985021530e-6 # kpc^3 Gyr^-2 Mo^-1

# Cosmology functions

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
    cp['rho_crit0'] = 277.53662719*pow(H0/100.0,2)
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

def Ecosmo(a,cp):
    return sqrt(cp['OmegaR0']*pow(a,-4) + cp['OmegaM0']*pow(a,-3) + cp['OmegaK0']*pow(a,-2) + cp['OmegaL0'])

def tcosmic(a,cp):
    integrand = lambda x: 1/(x*Ecosmo(x,cp))
    return integrate.quad(integrand,0,a)[0]*1e3/cp['H0']/km_per_s_2_kpc_per_Gyr

def D(a,cp):
    integrand = lambda x: 1/pow(x*Ecosmo(x,cp),3)
    return 2.5*cp['OmegaM0']*Ecosmo(a,cp)*integrate.quad(integrand,0,a)[0]

def Delta_DC(a,cp):
    return D(a,cp)*cp['Delta_DC0_NG']

def abox_EV(a,cp):
    return a*pow(1+Delta_DC(a,cp),-1/3.0)

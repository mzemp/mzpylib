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

km_per_s_2_kpc_per_Gyr = 1.02271216511
G = 4.4985021530e-6 # kpc^3 Gyr^-2 Mo^-1

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

# Cosmologial E(a) function

def Ecosmo(a,cp):
    return numpy.sqrt(cp['OmegaR0']*pow(a,-4) + cp['OmegaM0']*pow(a,-3) + cp['OmegaK0']*pow(a,-2) + cp['OmegaL0'])

# Cosmic time

def tcosmic(a,cp):
    integrand = lambda x: 1/(x*Ecosmo(x,cp))
    return scipy.integrate.quad(integrand,0,a)[0]*1e3/cp['H0']/km_per_s_2_kpc_per_Gyr

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

def rho_abc(r,alpha=1,beta=3,gamma=1,rs=1,rho0=1,rcutoff=numpy.inf,fdecay=0.3):

    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)
    rs = float(rs)
    rho0 = float(rho0)
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

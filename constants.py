#################################
#
# Some useful physical constants
#
# Written by Marcel Zemp
#
#################################

import numpy
import uncertainties

cl = {}

# Exact definitions

c_light = uncertainties.ufloat(299792458,0)
cl['c_light'] = [c_light,'m s^{-1}','Speed of light in vacuum','Definition']

au = uncertainties.ufloat(149597870700,0)
cl['au'] = [au,'m','Astronomical unit','Definition']

yr = uncertainties.ufloat(31557600,0)
cl['yr'] = [yr,'s','Julian year','Definition']

# Fundamental measured constants

G_Newton = uncertainties.ufloat_fromstr('6.67384(80)e-11')
cl['G_Newton'] = [G_Newton,'m^{3} s^{-2} kg^{-1}','Newton gravitational constant','NIST']

k_Boltzmann = uncertainties.ufloat_fromstr('1.3806488(13)e-23')
cl['k_Boltzmann'] = [k_Boltzmann,'J K^{-1}','Boltzmann constant','NIST']

uamu =  uncertainties.ufloat_fromstr('1.660538921(73)e-27')
cl['uamu'] = [uamu,'kg','Unified atomic mass unit','NIST']

m_proton = uncertainties.ufloat_fromstr('1.672621777(74)e-27')
cl['m_proton'] = [m_proton,'kg','Proton mass','NIST']

m_neutron = uncertainties.ufloat_fromstr('1.674927351(74)e-27')
cl['m_neutron'] = [m_neutron,'kg','Neutron mass','NIST']

m_electron = uncertainties.ufloat_fromstr('9.10938291(40)e-31')
cl['m_electron'] = [m_electron,'kg','Electron mass','NIST']

eV = uncertainties.ufloat_fromstr('1.602176565(35)e-19')
cl['eV'] = [eV,'J','Electron volt','NIST']

mu_sun = uncertainties.ufloat(1.32712440018e20,8e9)
cl['mu_sun'] = [mu_sun,'m^{3} s^{-2}','Solar standard gravitational parameter','JPL']

Mo = uncertainties.ufloat(1.98855e30,2.5e26)
cl['Mo'] = [Mo,'kg','Solar mass','Wikipedia']

Lo = uncertainties.ufloat(3.845e26,8e23)
cl['Lo'] = [Lo,'W','Solar luminosity','Binney & Tremaine']

# Derived constants

pc = au/numpy.tan(numpy.pi/(180*3600))
cl['pc'] = [pc,'m','Parsec','Derived']

lyr = c_light*yr
cl['lyr'] = [lyr,'m','Lightyear','Derived']

km_per_s_2_kpc_per_Gyr = 1e9*yr/pc
cl['km_per_s_2_kpc_per_Gyr'] = [km_per_s_2_kpc_per_Gyr,'','Conversion factor km s^{-1} => kpc Gyr^{-1}','Derived']

kpc_per_Gyr_2_km_per_s = pc/(1e9*yr)
cl['kpc_per_Gyr_2_km_per_s'] = [kpc_per_Gyr_2_km_per_s,'','Conversion factor kpc Gyr^{-1} => km s^{-1}','Derived']

G_Newton_cosmology = mu_sun*pow(1e9*yr,2)/pow(1e3*pc,3)
cl['G_Newton_cosmology'] = [G_Newton_cosmology,'kpc^{3} Gyr^{-2} Mo^{-1}','Newton gravitational constant for cosmology','Derived']

G_Newton_cosmology_v2 = G_Newton_cosmology*pow(kpc_per_Gyr_2_km_per_s,2)
cl['G_Newton_cosmology_v2'] = [G_Newton_cosmology_v2,'kpc km^{2} s^{-2} Mo^{-1}','Newton gravitational constant for cosmology (2nd version)','Derived']

rho_crit = (3*pow(1e-1/pc,2))/(8*numpy.pi*G_Newton)
cl['rho_crit'] = [rho_crit,'h^{2} kg m^{-3}','Critical density of the universe','Derived']

rho_crit_cosmology = (3*pow(1e8*yr/pc,2))/(8*numpy.pi*G_Newton_cosmology)
cl['rho_crit_cosmology'] = [rho_crit_cosmology,'h^{2} Mo kpc^{-3}','Critical density of the universe for cosmology','Derived']

GeV_per_c2_2_kg = 1e9*eV/pow(c_light,2)
cl['Gev_per_c2_2_kg'] = [GeV_per_c2_2_kg,'','Conversion factor GeV/c^{2} => kg','Derived']

GeV_per_c2_2_Mo = GeV_per_c2_2_kg/Mo
cl['Gev_per_c2_2_Mo'] = [GeV_per_c2_2_Mo,'','Conversion factor GeV/c^{2} => Mo','Derived']

GeV_per_c2_per_cm3_2_kg_per_m3 = GeV_per_c2_2_kg*1e6
cl['GeV_per_c2_per_cm3_2_kg_per_m3'] = [GeV_per_c2_per_cm3_2_kg_per_m3,'','Conversion factor GeV/c^{2} cm^{-3} => kg m^{-3}','Derived']

GeV_per_c2_per_cm3_2_Mo_per_pc3 = GeV_per_c2_2_Mo*pow(pc,3)*1e6
cl['GeV_per_c2_per_cm3_2_Mo_per_pc3'] = [GeV_per_c2_per_cm3_2_Mo_per_pc3,'','Conversion factor GeV/c^{2} cm^{-3} => Mo pc^{-3}','Derived']

# Give list of constants as output if run as main

if __name__ == "__main__":

    print
    print 'List of constants:'
    print
    for key, values in cl.items():
        print '%s: %s = %s %s = %.5e +/- %5e %s (%s)'%(values[2],key,values[0],values[1],values[0].n,values[0].s,values[1],values[3])
        print

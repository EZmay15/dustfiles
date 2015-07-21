
# coding: utf-8

# In[2]:

get_ipython().magic(u'matplotlib inline')
import dust
from scatmodelscopy import *
from cmindex import *
import matplotlib.pyplot as plt
from astropy.io import ascii
data = Mie()


# In[ ]:

#Silicate dust grains with radius 0.1
rad = 0.1

#Wavelengths from 0.01 - 1 microns for ^
tabl = ascii.read("silicate0.1")

wavelengths = tabl["w(micron)"]

#Conversion from wavelength(micron) to energy(keV)
print "Energy conversion done"
energies = [1.2398/(w * 1.e3) for w in wavelengths]

#Qsca
print "Qsca done"
scattering = []
for num in energies:
    scattering.append(data.Qsca(E = num, a = rad, cm = cmi.CmSilicate()))
    
#Qext
print "Qext done"
extinction = []
for num in energies:
    extinction.append(data.Qext(E = num, a = rad, cm = cmi.CmSilicate()))
    
#Qabs
print "Qabs done"
absorption = []
for num in range(len(energies)):
    absorption.append(float(extinction[num]-scattering[num]))


# In[ ]:

plt.xlabel("Energy", size = 12)
plt.ylabel("Efficiency", size = 12)
plt.plot(energies, scattering, "r")
plt.plot(energies, extinction, linewidth = 2.0)
plt.plot(energies, absorption, "g")
plt.loglog()

plt.show()


# In[ ]:




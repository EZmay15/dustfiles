
# coding: utf-8

# In[2]:

get_ipython().magic(u'matplotlib inline')
import dust
from scatmodelscopy import *
import matplotlib.pyplot as plt
import numpy as np
data = Mie()


# In[ ]:

#Silicate dust grains with radius 1
rad = 1.0

#Wavelengths from 0.1 - 1 microns for ^
tabl = ascii.read("silicate1.0")

wavelengths = tabl["w(micron)"]

#Conversion from wavelength(micron) to energy(keV)
energies = []
print "Energies"
for w in wavelengths:
    energies.append(1.2398/(w * 1000))

print "Qsca"
scattering = []
for num in energies:
    scattering.append(data.getQs(rad, num))
    
print "Qext"
extinction = []
for num in energies:
    extinction.append(data.getQs(rad, num, cm=cmi.CmSilicate(), getQ='ext'))

print "Qabs"
absorption = []
for num in energies:
    absorption.append(data.getQs(rad, num, cm=cmi.CmSilicate(), getQ='ext') - data.getQs(rad, num))


# In[ ]:

plt.xlabel("Energy", size = 12)
plt.ylabel("Efficiency", size = 12)
plt.plot(energies, scattering, linewidth = 2.0)
plt.plot(energies, extinction, "r")

plt.loglog()

plt.show()


# In[ ]:

plt.xlabel("Energy", size = 12)
plt.ylabel("Efficiency", size = 12)
plt.plot(energies, absorption, "g")
plt.show()


# In[ ]:




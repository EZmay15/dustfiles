
# coding: utf-8

# In[2]:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import ascii
from pylab import ylim

font = {'size'   : 15}
matplotlib.rc('font', **font)


# In[3]:

import scatmodels as sm
import cmindex as cm
import constants as c


# In[ ]:

mie_model = sm.Mie()


# In[ ]:

GRAIN_RAD = 0.5
tabl = ascii.read("graphite0.5")


# In[ ]:

WAVEL_UVOPT  = tabl["w(micron)"] # um, Draine's wavelengths
ENERGY_UVOPT = c.kev2lam() / (WAVEL_UVOPT * 1.e-4) # keV

HIGH_ENERGY  = np.linspace(0.1, 20.0, 30) # keV, 30 values between 0.1 and 20 keV 


# Show efficiencies for UV/Optical wavelengths

# In[ ]:

qext_uvopt_gra = mie_model.Qext( E=ENERGY_UVOPT, a=GRAIN_RAD, cm=cm.CmGraphite() )
qsca_uvopt_gra = mie_model.Qsca( E=ENERGY_UVOPT, a=GRAIN_RAD, cm=cm.CmGraphite() )


# In[ ]:

#plt.plot( WAVEL_UVOPT, qext_uvopt_gra, 'b-', label='Q_ext - Mie' )
#plt.plot( WAVEL_UVOPT, qsca_uvopt_gra, 'r-', label='Q_sca - Mie' ) 
#plt.plot( WAVEL_UVOPT, tabl["Q_abs"] + tabl["Q_sca"], "g-", label="Q_ext - Draine")
#plt.plot( WAVEL_UVOPT, tabl["Q_sca"], "y-", label="Q_sca - Draine")

#plt.loglog()
#plt.xlabel('Wavelength [um]')
#plt.ylabel('Q')
#plt.legend(loc='lower left')
#plt.figure(figsize=(20,10))
#plt.show()


# In[ ]:

#For extinction and scattering in all the uv_opt plots, create graphs that 
#quantify the difference between the two models:

#Percent difference = ( (Qext,Mie - Qext,Draine) / Qext,Draine ) x 100

qext_percent_diff = []
for i in range(len(tabl["w(micron)"])):
    value = ((qext_uvopt_gra[i] - (tabl["Q_abs"][i] + tabl["Q_sca"][i]))/
            (tabl["Q_abs"][i] + tabl["Q_sca"][i])) * 100
    qext_percent_diff.append(abs(value))
    
qsca_percent_diff = []
for i in range(len(tabl["w(micron)"])):
    value = ((qsca_uvopt_gra[i] - (tabl["Q_sca"][i]))/
            (tabl["Q_sca"][i])) * 100
    qsca_percent_diff.append(abs(value))

plt.plot( WAVEL_UVOPT, qext_percent_diff, 'b-', label = 'Q_ext')
plt.plot( WAVEL_UVOPT, qsca_percent_diff, 'r-', label = 'Q_sca')

plt.semilogx()
pylab.ylim([-10,100])
plt.xlabel('Wavelength [um]')
plt.ylabel('Difference [%]')
plt.legend(loc='upper center')
plt.show()


# Now do the high energy

# In[ ]:

qext_high_sil = mie_model.Qext( E=HIGH_ENERGY, a=GRAIN_RAD, cm=cm.CmSilicate() )
qsca_high_sil = mie_model.Qsca( E=HIGH_ENERGY, a=GRAIN_RAD, cm=cm.CmSilicate() )


# In[ ]:

plt.plot( HIGH_ENERGY, qext_high_sil, 'b-', label='Q_ext - Mie' )
plt.plot( HIGH_ENERGY, qsca_high_sil, 'r-', label='Q_sca - Mie' )

plt.loglog()
plt.xlabel('Energy [keV]')
plt.ylabel('Q')
plt.legend(loc='lower left')

plt.show()


# ## What's up with Drude approximation?
# 
# Let's peek at the **real part** (responsible for scattering) and the **imaginary part** (responsible for absorption) of the **complex index of refraction** to find out.

# In[ ]:

drude = cm.CmDrude()


# In[ ]:

plt.plot( ENERGY_UVOPT, drude.rp(ENERGY_UVOPT), 'r-', label='Real part' )
plt.plot( HIGH_ENERGY, drude.rp(HIGH_ENERGY), 'r-', label='' )

plt.plot( ENERGY_UVOPT, drude.ip(ENERGY_UVOPT)+1.0, 'b-', label='Imaginary part + 1.0' )
plt.plot( HIGH_ENERGY, drude.ip(HIGH_ENERGY)+1.0, 'b-', label='' )

plt.loglog()
plt.ylim(0.1, 1.e5)
plt.xlabel('Energy [keV]')


# Since the Drude approximation is for **scattering only**, the imaginary part is zero. (That's why the blue line, which is imaginary part plus one, is a flat line at $10^0 = 1.0$)

# In[ ]:

print drude.ip(HIGH_ENERGY)



# coding: utf-8

# In[2]:

get_ipython().magic(u'matplotlib inline')

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii


# In[3]:

tabl = ascii.read("silicate0.1")


# In[4]:

plt.xlabel("Wavelengths (microns)", size = 12)
plt.ylabel("Qabs", size = 12)
plt.plot(tabl["w(micron)"],tabl["Q_abs"], linewidth = 2.0)
plt.loglog()

plt.show()


# In[ ]:

plt.xlabel("Wavelengths (microns)", size = 12)
plt.ylabel("Qsca", size = 12)
plt.plot(tabl["w(micron)"],tabl["Q_sca"], linewidth = 2.0)
plt.loglog()

plt.show()


# In[ ]:

plt.xlabel("Wavelengths (microns)", size = 12)
plt.ylabel("Qext", size = 12)

#Adds Qsca and Qabs to get Qext
def findQext(Qsca, Qabs):
    total = []
    for num in range(0, len(Qsca)):
        total.append(float(Qsca[num]) + float(Qabs[num]))
    return total

Q_ext = findQext(tabl["Q_sca"], tabl["Q_abs"])

plt.plot(tabl["w(micron)"], Q_ext, linewidth = 2.0)
plt.plot(tabl["w(micron)"], tabl["Q_sca"], "g--")
plt.plot(tabl["w(micron)"],tabl["Q_abs"], "r--")
plt.loglog()

plt.show()


# In[ ]:




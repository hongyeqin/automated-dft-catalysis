#!/usr/bin/env python
# coding: utf-8

# In[24]:


from Vasp_default import adsorbates_urea
from ase.build import fcc111
from ase.build import add_adsorbate
from ase.io import write


# In[37]:


slab = fcc111('Al', size=(2,2,3))
slab.center(vacuum=10.0, axis=2)
mol = adsorbates_urea()
add_adsorbate(slab,mol['COOH'],height=2,position=(0,0))
write('POSCAR',slab,format='vasp')


# In[13]:





# In[ ]:





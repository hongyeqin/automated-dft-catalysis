#!/usr/bin/env python
# coding: utf-8

# In[1]:


from collections import OrderedDict


# In[17]:


def xc_settings(xc):
    xc_settings = {'lda': OrderedDict(pp='LDA'),
                   # GGAs
                   'gga': OrderedDict(pp='GGA'),
                   'pbe': OrderedDict(pp='PBE'),
                   'revpbe': OrderedDict(pp='LDA', gga='RE'),
                   'rpbe': OrderedDict(gga='RP', pp='PBE'),
                   'am05': OrderedDict(pp='LDA', gga='AM'),
                   'pbesol': OrderedDict(gga='PS', pp='PBE'),
                   # Meta-GGAs
                   'tpss': OrderedDict(pp='PBE', metagga='TPSS'),
                   'revtpss': OrderedDict(pp='PBE', metagga='RTPSS'),
                   'm06l': OrderedDict(pp='PBE', metagga='M06L'),
                   # vdW-DFs
                   'optpbe_vdw': OrderedDict(pp='LDA', gga='OR', luse_vdw=True,
                                             aggac=0.0),
                   'optb88_vdw': OrderedDict(pp='LDA', gga='BO', luse_vdw=True,
                                             aggac=0.0, param1=1.1 / 6.0,
                                             param2=0.22),
                   'optb86b_vdw': OrderedDict(pp='LDA', gga='MK',
                                              luse_vdw=True, aggac=0.0,
                                              param1=0.1234, param2=1.0),
                   'vdw_df2': OrderedDict(pp='LDA', gga='ML', luse_vdw=True,
                                          aggac=0.0, zab_vdw=-1.8867),
                   'beef_vdw': OrderedDict(pp='PBE', gga='BF', luse_vdw=True,
                                           zab_vdw=-1.8867, lbeefens=True),
                   # hybrids
                   'pbe0': OrderedDict(pp='LDA', gga='PE', lhfcalc=True),
                   'hse03': OrderedDict(pp='LDA', gga='PE', lhfcalc=True,
                                        hfscreen=0.3),
                   'hse06': OrderedDict(pp='LDA', gga='PE', lhfcalc=True,
                                        hfscreen=0.2),
                   'b3lyp': OrderedDict(pp='LDA', gga='B3', lhfcalc=True,
                                        aexx=0.2, aggax=0.72, aggac=0.81,
                                        aldac=0.19),
                   'hf': OrderedDict(pp='PBE', lhfcalc=True, aexx=1.0,
                                     aldac=0.0, aggac=0.0)} 
    return xc_settings[xc]
    


# In[27]:


def molecule_settings():
    """
    The default settings used to do DFT calculatins of molecules
    
    """
    molecule_settings = OrderedDict(max_atoms=80,
                                    vasp=OrderedDict(
                                                     ibrion=2,
                                                     nsw=200,
                                                     isif=2,
                                                     kpts=(1, 1, 1,),
                                                     ediff=0.000001,
                                                     encut=600,
                                                     ismear=0,
                                                     sigma=0.01,
                                                     dipole=True,
                                                     idipole=4,
                                                     prec='Accurate',
                                                     **xc_settings()))
    return molecule_settings


# In[30]:


def bulk_settings():
    ''' The default settings we use to do DFT calculations of bulks '''
    bulk_settings = OrderedDict(max_atoms=80,
                                vasp=OrderedDict(ibrion=2,
                                                 nsw=600,
                                                 isif=3,
                                                 isym=0,
                                                 ediff=1e-6,
                                                 kpts=(10, 10, 10),
                                                 prec='Accurate',
                                                 encut=500.,
                                                 **xc_settings('pbesol')))
    return bulk_settings

def surface_slab_settings():
    '''
    The default settings we use to do DFT calculations of bulks
    spefically for surface energy calculations.
    '''
    Surface_slab_settings = OrderedDict(max_atoms=80,
                                        max_miller=2,
                                   vasp=OrderedDict(ibrion=2,
                                                    nsw=600,
                                                    isif=2,
                                                    isym=0,
                                                    ediff=1e-8,
                                                    kpts=(3, 3, 1),
                                                    prec='Accurate',
                                                    encut=500.,
                                                    **xc_settings('pbe')))
    return Surface_slab_settings


# In[44]:


from ase.collections import g2
from ase.data.pubchem import pubchem_atoms_search
from ase import Atoms
def adsorbates_urea():
    Urea = pubchem_atoms_search(name='urea')
    Methanol  = g2['CH3OH']
    adsorbates = {}
    adsorbates['Urea'] = Urea
    adsorbates['CON2H3'] = Atoms(numbers=[6,8,7,7,1,1,1],positions=[(-0.0749,-0.0001,-0.0003),
                                                                    (-1.3042,-0.0008,0.0001),
                                                                    (0.6903,-1.1479,0.0001),
                                                                    (0.6888,1.1489,0.0001),
                                                                    (1.7041,-1.1111,-0.0002),
                                                                    (0.2605,-2.0669,0.0001),
                                                                    (0.2578,2.0672,0.0002)])
    adsorbates['CON2H2'] = Atoms(numbers=[6,8,7,7,1,1],positions=[(-0.0749,-0.0001,-0.0003),
                                                                    (-1.3042,-0.0008,0.0001),
                                                                    (0.6903,-1.1479,0.0001),
                                                                    (0.6888,1.1489,0.0001),
                                                                    (1.7041,-1.1111,-0.0002),
                                                                    (0.2605,-2.0669,0.0001)])
    adsorbates['CON2H'] = Atoms(numbers=[6,8,7,7,1],positions=[(-0.0749,-0.0001,-0.0003),
                                                                    (-1.3042,-0.0008,0.0001),
                                                                    (0.6903,-1.1479,0.0001),
                                                                    (0.6888,1.1489,0.0001),
                                                                    (1.7041,-1.1111,-0.0002)])
    adsorbates['CON2'] = Atoms(numbers=[6,8,7,7],positions=[(-0.0749,-0.0001,-0.0003),
                                                                    (-1.3042,-0.0008,0.0001),
                                                                    (0.6903,-1.1479,0.0001),
                                                                    (0.6888,1.1489,0.0001)])
    adsorbates['CO'] = Atoms(numbers=[6,8],positions=[(-0.0749,-0.0001,-0.0003),
                                                      (-1.3042,-0.0008,0.0001)])
    adsorbates['CO2'] = pubchem_atoms_search(name='CO2')
    adsorbates['COOH'] = Atoms(numbers=[6,8,8,1],positions=[(-0.0749,-0.0001,-0.0003),
                                                            (-1.3042,-0.0008,0.0001),
                                                            (-2.6000,-0.0008,0.0001),
                                                            (-3.5400,-0.0008,0.0001)])
    
        
    return adsorbates

def adsorbates_water():
    
    adsorbates = {}
    adsorbates['H'] = Atoms('H', positions=[[0., 0., -0.5]])
    adsorbates['O'] = Atoms('O')
    adsorbates['OH'] = Atoms('OH', positions=[[0., 0., 0.],
                                              [0.92, 0., 0.32]])
    adsorbates['OOH'] = Atoms('OOH', positions=[[0., 0., 0.],
                                                [1.28, 0., 0.67],
                                                [1.44, -0.96, 0.81]])
    return adsorbates


# In[56]:





# In[ ]:





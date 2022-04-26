#!/usr/bin/env python
# coding: utf-8


import os
from operator import pos
import xml.etree.ElementTree as ET
from pprint import pprint
from ase import Atoms, symbols
from ase.atom import Atom
from ase.atoms import default
from ase.units import create_units
from ase.utils import iofunction
import numpy as np



first_title = ['generator','incar','calculation','kpoints','atominfo']
"""        
            
            
            
                                   -------- "partial"-----"array" ----"set"---"set comment="ion 1-n""----'set comment="spin1/2"'
                                  |                                           return           "s,py,pz,px,dxy,dyz,dz2,dxz,dx2"
                                  |
                                  |
                                  |
            ------------ "dos"---- -------- "i@name="efermi""
            |                     |
                                  |                                -----
            |                     |                               |
            |                     |                               |
            |                      --------"total" -------"array"  -----"set"----"set comment="spin1"
                                                                  |
                                                                  |           
                                                                   ------
            |                               --------"array"
            |                              |
            |                              |
                                           |            
calculation -------------"eigenvalues" ---- --------
            |                              |
            |                              | 
            |                              | 
            |                               --------
            |
            |
            |
             ------------"scstep"
  




"""


def xml_parser(xmlfile):
    """
    """
    pass
tree = ET.parse('vasprun.xml')
root = tree.getroot()
path = '/'.join(['atominfo','array[@name="atoms"]','set','rc','c'])
#print(len(tree.findall(path)))
#print(len(tree.findall(path)[1::2]))



"""                                                  ------------ "v[@name='division']"
                                                    |
                                                    |
                                                    |
           |                                         |
            --------"generation[@param="Gamma"]" ----
           |
           |
           |
           |
           |
kpoints----


"""



"""        
            -------"atoms"
           |-------"type"
atominfo---|-------"array [name="atoms"]"--set--rc-c
           |-------"array [name="atomtypes"]"--set--rc-c # atoms number, element, mass, valence, pseudopotntial


"""
def get_atominfo(directory,filename):
    with open(os.path.join(directory,
                           filename)) as f:
        tree = ET.parse(f)
        atom_dic = {}
        key_list = []
        value_list = []
        atom_name_path = '/'.join(['atominfo','array[@name="atoms"]','set','rc','c'])
        total_elements = len(tree.findall(atom_name_path)) #contain elements and atom type  Ni: 1, W:2  , N:3 . 
        key = tree.findall(atom_name_path)[::2]
        value = tree.findall(atom_name_path)[1::2]
        
        for i in range(len(key)):
            atom_dic[key[i].text] = value[i].text
            key_list.append(key[i].text)
            value_list.append(value[i].text)
        file = open('atom_dict.txt','w')
        file.write(str(atom_dic) + '\n'+ str(key_list) + '\n'+ str(value_list))
        file.close()
        
        
    return atom_dic, key_list, value_list
        
        
    
        
#print(get_atominfo('./','vasprun.xml'))    




"""





    
             ---------------------- 
            | "i[@name='version']"                    
            |----------------------                    
            | "i[@name='program']" 
            |-----------------------
            |"i[@name='subversion']"
"generator" |-----------------------
            |"i[@name='platform']"
            |-----------------------
            |"i[@name='date']" 
            |-----------------------   
            |"i[@name='time']" 
             ------------------------                                                   
                                                                            
                                                                   
 
            
"""

def get_version(directory,filename):
    with open(os.path.join(directory,
                           filename)) as f:
        tree = ET.parse(f)
        program_path = '/'.join(['generator',"i[@name='program']"])
        version_path = '/'.join(['generator',"i[@name='version']"])
        subversion_path = '/'.join(['generator',"i[@name='subversion']"])
        platformtype_path = '/'.join(['generator',"i[@name='platform']"])
        date_path = '/'.join(['generator',"i[@name='date']"])
        time_path = '/'.join(['generator',"i[@name='time']"])
        program = tree.find(program_path).text
        version = tree.find(version_path).text
        subversion = tree.find(subversion_path).text
        platformtype = tree.find(platformtype_path).text
        date = tree.find(date_path).text
        time = tree.find(time_path).text
        print("program is {},version {},runing start date {} {}".format(program,version,date,time))
        
        
def get_fermi_level(directory,filename):
    """Return the Fermi level."""
    

    with open(os.path.join(directory,
                           filename)) as f:
        tree = ET.parse(f)
        path = '/'.join(['calculation',
                         'dos',
                         "i[@name='efermi']"
                         ])
        return float(tree.find(path).text)

def get_orbital_order(directory,filename):
    with open(os.path.join(directory,
                           filename)) as f:
         path = '/'.join(['calculation', 'dos',
                     'partial',
                     'array',
                     'field'])
         tree = ET.parse(f)
         orbital_order = [str(el.text).strip()  for el in tree.findall(path)]
    return orbital_order
a = get_orbital_order('./','vasprun.xml' ) 



orbitals  =["s","px","py","pz","dxy","dyz","dxz","dz2","dx2",
                       "f1","f2","f3","f4","f5","f6","f7"]
orbitals_p=["px","py","pz"]
orbitals_d=["dxy","dyz","dxz","dz2","dx2"]
orbitals_f=["f1","f2","f3","f4","f5","f6","f7"]

        
def get_pdos(directory,filename, atom_index, orbital=None, spin=1, efermi=None, *args):
    """ 

                                                                    | total  
                                                   -----------field--energy         
                 --------"dos[@iname="efermi"]"   |                  |integrated
                |                                 |
                |--------'total'----'array'------- ----------set--set[@comment="spin {}"]'
                |  
                |                               
                |                                 
calculations--dos--------'partial'---'array'---'set'---'set[@comment="ion {}"]'
                |                       |               |
                |                       |              'set[@comment="spin {}"]'
                |                       |
                |                       -------'field'--{'energy, 's, 'px','py',...etc}

   orbital type: 's' or ['dxy','dyz','dxz',dx2,'dz2']

    """
    
    
    with open(os.path.join(directory,
                           filename)) as f:
         path = "/".join(['calculation', 'dos',
                     'partial',
                     'array',
                     'set',
                     'set[@comment="ion {}"]'.format(atom_index),
                     'set[@comment="spin {}"]'.format(spin),
                     "r"])
         field_array = get_orbital_order(directory, filename)
         orbital_array = get_orbital_order(directory, filename)[1:]
         if efermi is None:
             efermi = get_fermi_level(directory, filename)
         else:
             efermi = 0.0
             
         tree = ET.parse(f)
         results = [[float(x) for x in el.text.split()] for el in tree.findall(path)]
         
         energy = np.array([x[field_array.index('energy')] for x in results]) - efermi
         
         dos = []
         for i, j in enumerate(orbital_array):
             dos.append([x[field_array.index(orbital_array[i])] for x in results])
             
         total_dos = np.sum(np.array(dos), axis=0) 
         
                        
         if orbital is None:
            dos = []
            for i, j in enumerate(orbital_array):
                dos.append([x[field_array.index(orbital_array[i])] for x in results]) 
            total_dos = np.sum(np.array(dos), axis=0)
            return [energy, total_dos]
        
         if isinstance(orbital, list):
            dos = []
            for i, j in enumerate(orbital):
                dos.append([x[field_array.index(orbital[i])] for x in results]) 
            total_dos = np.sum(np.array(dos), axis=0)
            return [energy, total_dos]
         
         if orbital == 's':
             dos = np.array([x[orbital_array.index('s')] for x in results])
             return [energy, dos]
         if orbital == 'px':
             dos = np.array([x[orbital_array.index('px')] for x in results])
             return [energy, dos]
                   
         if orbital == 'py':
             dos = np.array([x[orbital_array.index('py')] for x in results])
             return [energy, dos]
                         
         if orbital == 'pz':
             dos = np.array([x[orbital_array.index('pz')] for x in results])
             return [energy, dos]
         
         if orbital == 'dxy':
             dos = np.array([x[orbital_array.index('dxy')] for x in results])
             return [energy, dos]
         
         if orbital == 'dyz':
             dos = np.array([x[orbital_array.index('dyz')] for x in results])
             return [energy, dos]
         
         if orbital == 'dxz':
             dos = np.array([x[orbital_array.index('dxz')] for x in results]) 
             return [energy, dos]
         
         if orbital == 'dz2':
             dos = np.array([x[orbital_array.index('dz2')] for x in results]) 
             return [energy, dos]
         
         if orbital == 'dx2':
             dos = np.array([x[orbital_array.index('x2-y2')] for x in results]) 
             return [energy, dos]
         
         if orbital == 'f1':
             dos = np.array([x[orbital_array.index('f1')] for x in results])
             return [energy, dos]
         
         if orbital == 'f2':
             dos = np.array([x[orbital_array.index('f2')] for x in results]) 
             return [energy, dos]

         if orbital == 'f3':
             dos = np.array([x[orbital_array.index('f3')] for x in results])
             return [energy, dos]
         
         if orbital == 'f4':
             dos = np.array([x[orbital_array.index('f4')] for x in results])
             return [energy, dos]
         
         if orbital == 'f5':
             dos = np.array([x[orbital_array.index('f5')] for x in results])
             return [energy, dos]
         
         if orbital == 'f6':
             dos = np.array([x[orbital_array.index('f6')] for x in results])
             return [energy, dos]
         
         if orbital == 'f7':
             dos = np.array([x[orbital_array.index('f7')] for x in results])
             return [energy, dos]
         #return energy
            

        
#a = get_pdos('./','vasprun.xml',1, orbital=['s','px'])   
#b = get_orbital_order('./','vasprun.xml') 
#print(a)












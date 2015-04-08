#!/usr/bin/python

# Needed libraries

import urllib
import urllib.request
import math
import matplotlib.pyplot as plt
import numpy as np
import io
from bs4 import BeautifulSoup as BS
from scipy.interpolate import interp1d

# Expressions needed later (easy typing)

def col(name, zahl):
    return name[name.dtype.names[zahl]]

# Get the data for the desired gas species

url = urllib.request.urlopen('http://webbook.nist.gov/chemistry/fluid/')
open_url = url.read().decode('UTF-8')
html = BS(open_url)


# Define Nozzle properties and gas species

Gas = 'H2'

P = 40
T_desired = 33

THigh = T_desired + 7
TLow = T_desired - 13

Species = {'H2': 'Hydrogen', 'N2': 'Nitrogen', 'He': 'Helium', 'Ne': 'Neon', 'Ar': 'Argon', 'Kr': 'Krypton',
           'Xe': 'Xenon'}

# Find the value for the desired gas

gas_species = Species[Gas]
ID = html.find(text=gas_species).find_parent('option')['value']
print(ID)

# Define the desired type of data

type_of_data = 'Isobaric properties'

TYPE = html.find(text=type_of_data).find_parent('td').find_previous_sibling()('input')[0]['value']
print(TYPE)

# Defining the default URL values

DEFAULT = 'http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on'
UNITS = '&Digits=5&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm'

# URL where the data can be found

URL = DEFAULT + '&ID=' + ID + '&Type=' + TYPE + '&P=' + str(P) + '&THigh=' + str(THigh) + '&TLow=' + str(
    TLow) + '&TInc=0.001' + UNITS
print(URL)

# Open the web source

source = urllib.request.urlopen(URL)
data = source.read()  # .decode('UTF-8')
source.close()

# print(type(data))
# print(data)

# Data manipulation in order to be readable

matrixdata = np.genfromtxt(io.BytesIO(data), delimiter='\t', names=True, dtype=None, deletechars='', replace_space='')

# print(matrixdata)

# Adiabatic index kappa = c_p / c_V for gas species

kappa = {'H2': 1.41, 'N2': 1.4, 'He': 1.67, 'Ne': 1.67, 'Ar': 1.67, 'Kr': 1.67, 'Xe': 1.67}

# Unit mass of the gas species [kg]

m_u = 1.660538e-27

# Boltzman constant [J/K]

k_C = 1.3806488e-23

# Mass numbers of the gas species

m_Nr = {'H2': 2 * 1.008, 'N2': 2 * 14.0067, 'He': 4.002602, 'Ne': 20.1797, 'Ar': 39.948, 'Kr': 83.798, 'Xe': 131.293}

# Calculate the supersonic gas expansion velocities

v_gas = list(
    math.sqrt(((2 * kappa[Gas]) / (kappa[Gas] - 1)) * (k_C / (m_u * m_Nr[Gas])) * col(matrixdata, 0)[i]) for i in
    range(len(col(matrixdata, 0))))

# Calculate the alternative supersonic gas expansion velocities

kappa_exact = list(col(matrixdata, 8) / col(matrixdata, 7))
v_gas_alt = list(
    math.sqrt(((2 * kappa_exact[i]) / (kappa_exact[i] - 1)) * (k_C / (m_u * m_Nr[Gas])) * col(matrixdata, 0)[i]) for i
    in range(len(col(matrixdata, 0))))

# Define the triple point values of the different gases required for the calculation of the supercritical velocities
# Enthalpy in [kJ/kg]

h_l = {'H2': -53.923, 'N2': -150.75, 'He': -7.5015, 'Ne': -4.9151, 'Ar': -121.44, 'Kr': -2.0639, 'Xe': -1.2418}
h_v = {'H2': 399.82, 'N2': 46.211, 'He': 15.726, 'Ne': 83.204, 'Ar': 42.281, 'Kr': 106.34, 'Xe': 95.165}

# Calculate supercritical velocities

h_mean = 0.5 * (h_l[Gas] + h_v[Gas])

i = 0
v_Sup1 = list()
while h_mean > col(matrixdata, 5)[i]:
    v_Sup1.append(0)
    i = i + 1
    if i == len(col(matrixdata, 0)):
        v_Sup2 = []
        break
else:
    v_Sup2 = list((math.sqrt(2000 * (col(matrixdata, 5)[i] - h_mean))) for i in range(i, len(col(matrixdata, 0))))

v_Sup = v_Sup1 + v_Sup2
# v_Sup=list((math.sqrt(2000*(col(matrixdata,5)[i]-h_mean))) for i in range(len(col(matrixdata,0))))
# print(len(col(matrixdata,0)))
# print(v_Sup)
# len(v_Sup)

# Calculate liquid velocities

v_Bernoulli = list(
    (math.sqrt((2 * 1e5 * col(matrixdata, 1)[i]) / col(matrixdata, 2)[i])) for i in range(len(col(matrixdata, 0))))
T0 = list(col(matrixdata, 0))
# len(T0)

# Plot velocity curves

plt.plot(T0, v_Sup, 'r')
plt.plot(T0, v_Bernoulli, 'b')
plt.plot(T0, v_gas, 'y')
plt.plot(T0, v_gas_alt, 'g')
plt.xlabel('Temperature [K]')
plt.ylabel('Velocity [m/s]')
plt.show()

vgas = interp1d(T0, v_gas_alt)
vSup = interp1d(T0, v_Sup)
vBer = interp1d(T0, v_Bernoulli)
vgas2 = interp1d(T0, v_gas)

float(vgas(T_desired))

float(vSup(T_desired))

float(vgas2(T_desired))

float(vBer(T_desired))

# plt.plot(T0, v_Sup)

Tnew = np.linspace(20, 70, 100)
plt.plot(T0, vBer(T0), 'b')
plt.plot(T0, vgas(T0), 'g')
plt.plot(T0, vgas2(T0), 'y')
plt.plot(T0, vSup(T0), 'r')
plt.xlabel('Temperature [K]')
plt.ylabel('Velocity [m/s]')
plt.show()




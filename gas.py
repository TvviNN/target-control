#!/usr/bin/python

# Needed libraries
from __future__ import division
from urllib.request import urlopen
from numpy import sqrt, genfromtxt
import numpy as np
from io import BytesIO
import os.path
from os import makedirs
from numpy.lib.scimath import sqrt as csqrt
import matplotlib.pyplot as plt

# Needed constants
# Unit mass of the gas species [kg]
M_u = 1.660538e-27

# Boltzmann constant [J/k]
K_b = 1.3806488e-23


class Gas:
    def __init__(self, temp=0, press=0, iso_type='IsoBar', t_low=0, t_high=0, p_low=0, p_high=0):
        self.temp = temp
        self.press = press
        self.iso_type = iso_type
        self.t_low = t_low
        self.t_high = t_high
        self.p_low = p_low
        self.p_high = p_high
        self.long_name = ''
        self.short_name = ''
        self.kappa = 0
        self.mass = 0
        self.liquid_enthalpy_tp = 0
        self.vapour_enthalpy_tp = 0
        self.gas_id = ''
        self.local_directory = os.path.expanduser('~') + '/NIST-gas_database/'
        if self.iso_type == 'IsoBar':
            pass
        elif self.iso_type == 'IsoTherm':
            pass
        else:
            self.iso_type = 'IsoBar'

    def get_file_name(self):
        if self.iso_type == 'IsoBar':
            file_name = self.short_name + '_IsoBar_' + str(self.temp) + '_K_' + str(self.t_low) + '-' + str(
                self.t_high) + '_K_at_' + str(self.press) + '_bar'
        elif self.iso_type == 'IsoTherm':
            file_name = self.short_name + '_IsoTherm_' + str(self.press) + '_bar_' + str(self.p_low) + '-' + str(
                self.p_high) + '_bar_at_' + str(self.temp) + '_K'
        else:
            self.iso_type = 'IsoBar'
            file_name = self.short_name + '_IsoBar_' + str(self.temp) + '_K_' + str(
                self.t_low) + '-' + str(self.t_high) + '_K_at_' + str(self.press) + '_bar'
        return file_name

    def get_long_name(self):
        """The long name of the chosen gas species."""
        return self.long_name

    def get_short_name(self):
        """The short name of the chosen gas species."""
        return self.short_name

    def get_kappa(self):
        """The adiabatic index [c_p/c_V] of the chosen gas species."""
        return self.kappa

    def get_mass(self):
        """The mass number of the gas particles (atoms / molecules)."""
        return self.mass

    def get_liquid_enthalpy_tp(self):
        """The liquid enthalpy [kJ/kg] of the chosen gas species at its triple point."""
        return self.liquid_enthalpy_tp

    def get_vapour_enthalpy_tp(self):
        """The vapour enthalpy [kJ/kg] of the chosen gas species at its triple point."""
        return self.vapour_enthalpy_tp

    def get_gas_id(self):
        """The individual gas ID number, found in http://webbook.nist.gov/chemistry/fluid/"""
        return self.gas_id

    def get_nist_data_url(self):
        """Define the URL where the raw data can be found."""
        url_default = 'http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on'
        url_units = '&Digits=5&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm'
        if self.iso_type == 'IsoBar':
            url_nist = url_default + '&ID=' + self.gas_id + '&Type=' + self.iso_type + '&P=' + str(
                self.press) + '&THigh=' + str(self.t_high) + '&TLow=' + str(
                self.t_low) + '&TInc=0.001' + url_units
        elif self.iso_type == 'IsoTherm':
            url_nist = url_default + '&ID=' + self.gas_id + '&Type=' + self.iso_type + '&T=' + str(
                self.temp) + '&PHigh=' + str(self.p_high) + '&PLow=' + str(self.p_low) + '&PInc=0.001' + url_units
        else:
            self.iso_type = 'IsoBar'
            url_nist = url_default + '&ID=' + self.gas_id + '&Type=' + self.iso_type + '&P=' + str(
                self.press) + '&THigh=' + str(self.t_high) + '&TLow=' + str(
                self.t_low) + '&TInc=0.001' + url_units
        return url_nist

    def download_raw_nist_data(self):
        """Download the raw gas table data from the NIST database."""
        with urlopen(self.get_nist_data_url()) as source_url:
            raw_nist_data = source_url.read()
        return raw_nist_data

    def save_file(self):
        """Save the file to the defined location (self.local_directory)."""
        if os.path.isfile(self.local_directory + self.get_file_name()):
            pass
        else:
            with open(self.local_directory + self.get_file_name(), 'w') as data_file:
                data_file.write(self.download_raw_nist_data().decode(encoding='UTF-8'))
        return

    def save_gas_data(self):
        """Check the path and save the file."""
        if os.path.exists(self.local_directory):
            self.save_file()
        else:
            makedirs(self.local_directory)
            self.save_file()
        return

    def get_raw_nist_data(self):
        """Get the raw gas table data from the NIST database (local or: from the web - then save it!)."""
        if os.path.isfile(self.local_directory + self.get_file_name()):
            with open(self.local_directory + self.get_file_name(), 'r') as data_file:
                raw_nist_data = data_file.read().encode(encoding='UTF-8')
        else:
            raw_nist_data = self.download_raw_nist_data()
            self.save_gas_data()
        return raw_nist_data

    def get_table_nist_data(self):
        """Transform the raw NIST data into a readable format."""
        table_nist_data = genfromtxt(BytesIO(self.get_raw_nist_data()), skip_header=1, delimiter='\t', names=None,
            dtype=None, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
            deletechars='', replace_space='')
        return table_nist_data

    def get_column_data(self, column_number):
        """Get only the data of one particular column of the NIST gas data."""
        if self.get_table_nist_data().ndim >= 2:
            column_data = self.get_table_nist_data()[:, column_number]
        else:
            column_data = self.get_table_nist_data()[column_number]
        return column_data

    def get_column_data_old(self, column_number):
        """Get only the data of one particular column of the NIST gas data [[deprecated!]]."""
        return self.get_table_nist_data()[self.get_table_nist_data().dtype.names[column_number]]

    def get_temperature_column(self):
        return self.get_column_data(0)

    def get_pressure_column(self):
        return self.get_column_data(1)

    def get_density_column(self):
        return self.get_column_data(2)

    def get_enthalpy_column(self):
        return self.get_column_data(5)

    def get_cv_column(self):
        return self.get_column_data(7)

    def get_cp_column(self):
        return self.get_column_data(8)

    def get_simple_supersonic_velocity_value(self):
        """Calculate the supersonic velocity for the given temperature (simple calculation w/o exact kappa-values)."""
        v = sqrt(((2 * self.kappa) / (self.kappa - 1)) * (K_b / (M_u * self.mass)) * self.temp)
        return v

    def get_simple_supersonic_velocity_array(self):
        """Calculate the supersonic velocity for the given temperature range (simple calculation w/o exact kappa-values)."""
        v_array = sqrt(
            ((2 * self.kappa) / (self.kappa - 1)) * (K_b / (M_u * self.mass)) * self.get_temperature_column())
        return v_array

    def get_exact_kappa(self):
        """Calculate exact kappa value from NIST data."""
        kappa_exact = self.get_cp_column() / self.get_cv_column()
        return kappa_exact

    def get_exact_supersonic_velocity_value(self):
        """Calculate the supersonic velocity for the given temperature (with exact kappa values)."""
        self.t_low = self.temp
        self.t_high = self.temp
        v = sqrt(((2 * self.get_exact_kappa()) / (self.get_exact_kappa() - 1)) * (K_b / (M_u * self.mass)) * self.temp)
        return v

    def get_exact_supersonic_velocity_array(self):
        """Calculate the supersonic velocity for the given temperature range (with exact kappa values)."""
        v_array = sqrt(((2 * self.get_exact_kappa()) / (self.get_exact_kappa() - 1)) * (
            K_b / (M_u * self.mass)) * self.get_temperature_column())
        return v_array

    def get_mean_enthalpy_at_tp(self):
        """Calculate the mean enthalpy at the triple point of the gas."""
        mean_enthalpy = 0.5 * (self.get_liquid_enthalpy_tp() + self.get_vapour_enthalpy_tp())
        return mean_enthalpy

    def get_supercritical_velocity_array(self):
        """Calculate the supercritical velocity for the given temperature or pressure range."""
        v_array_calculation = csqrt(2 * 1000 * (self.get_enthalpy_column() - self.get_mean_enthalpy_at_tp()))
        v_array = np.where(np.iscomplex(v_array_calculation), 0, v_array_calculation).astype(float)
        return v_array

    def get_liquid_velocity_array(self):
        """Calculate the velocity of an expanding liquid with the Bernoulli formula."""
        v_array = sqrt((2 * 1e5 * self.get_pressure_column()) / self.get_density_column())
        return v_array


class N2(Gas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.temp == 0:
            self.temp = 140
        if self.iso_type == 'IsoBar':
            if self.t_low == 0:
                self.t_low = self.temp - 1
            if self.t_high == 0:
                self.t_high = self.temp + 1
        if self.press == 0:
            self.press = 10
        if self.iso_type == 'IsoTherm':
            if self.p_low == 0:
                self.p_low = self.press - 1
            if self.p_high == 0:
                self.p_high = self.press + 1
        self.long_name = 'Nitrogen'
        self.short_name = 'N2'
        self.kappa = 1.4
        self.mass = 2 * 14.0067
        self.liquid_enthalpy_tp = -150.75
        self.vapour_enthalpy_tp = 46.211
        self.gas_id = 'C7727379'


class H2(Gas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.temp == 0:
            self.temp = 40
        if self.iso_type == 'IsoBar':
            if self.t_low == 0:
                self.t_low = self.temp - 1
            if self.t_high == 0:
                self.t_high = self.temp + 1
        if self.press == 0:
            self.press = 40
        if self.iso_type == 'IsoTherm':
            if self.p_low == 0:
                self.p_low = self.press - 1
            if self.p_high == 0:
                self.p_high = self.press + 1
        self.long_name = 'Hydrogen'
        self.short_name = 'H2'
        self.kappa = 1.41
        self.mass = 2 * 1.008
        self.liquid_enthalpy_tp = -53.923
        self.vapour_enthalpy_tp = 399.82
        self.gas_id = 'C1333740'


class He(Gas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.temp == 0:
            self.temp = 10
        if self.iso_type == 'IsoBar':
            if self.t_low == 0:
                self.t_low = self.temp - 1
            if self.t_high == 0:
                self.t_high = self.temp + 1
        if self.press == 0:
            self.press = 5
        if self.iso_type == 'IsoTherm':
            if self.p_low == 0:
                self.p_low = self.press - 1
            if self.p_high == 0:
                self.p_high = self.press + 1
        self.long_name = 'Helium'
        self.short_name = 'He'
        self.kappa = 1.67
        self.mass = 4.002602
        self.liquid_enthalpy_tp = -7.5015
        self.vapour_enthalpy_tp = 15.726
        self.gas_id = 'C7440597'


class Ne(Gas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.temp == 0:
            self.temp = 120
        if self.iso_type == 'IsoBar':
            if self.t_low == 0:
                self.t_low = self.temp - 1
            if self.t_high == 0:
                self.t_high = self.temp + 1
        if self.press == 0:
            self.press = 10
        if self.iso_type == 'IsoTherm':
            if self.p_low == 0:
                self.p_low = self.press - 1
            if self.p_high == 0:
                self.p_high = self.press + 1
        self.long_name = 'Neon'
        self.short_name = 'Ne'
        self.kappa = 1.67
        self.mass = 20.1797
        self.liquid_enthalpy_tp = -4.9151
        self.vapour_enthalpy_tp = 83.204
        self.gas_id = 'C7440019'


class Ar(Gas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.temp == 0:
            self.temp = 180
        if self.iso_type == 'IsoBar':
            if self.t_low == 0:
                self.t_low = self.temp - 1
            if self.t_high == 0:
                self.t_high = self.temp + 1
        if self.press == 0:
            self.press = 10
        if self.iso_type == 'IsoTherm':
            if self.p_low == 0:
                self.p_low = self.press - 1
            if self.p_high == 0:
                self.p_high = self.press + 1
        self.long_name = 'Argon'
        self.short_name = 'Ar'
        self.kappa = 1.67
        self.mass = 39.948
        self.liquid_enthalpy_tp = -121.44
        self.vapour_enthalpy_tp = 42.281
        self.gas_id = 'C7440371'


class Kr(Gas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.temp == 0:
            self.temp = 300
        if self.iso_type == 'IsoBar':
            if self.t_low == 0:
                self.t_low = self.temp - 1
            if self.t_high == 0:
                self.t_high = self.temp + 1
        if self.press == 0:
            self.press = 10
        if self.iso_type == 'IsoTherm':
            if self.p_low == 0:
                self.p_low = self.press - 1
            if self.p_high == 0:
                self.p_high = self.press + 1
        self.long_name = 'Krypton'
        self.short_name = 'Kr'
        self.kappa = 1.67
        self.mass = 83.798
        self.liquid_enthalpy_tp = -2.0639
        self.vapour_enthalpy_tp = 106.43
        self.gas_id = 'C7439909'


class Xe(Gas):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.temp == 0:
            self.temp = 300
        if self.iso_type == 'IsoBar':
            if self.t_low == 0:
                self.t_low = self.temp - 1
            if self.t_high == 0:
                self.t_high = self.temp + 1
        if self.press == 0:
            self.press = 10
        if self.iso_type == 'IsoTherm':
            if self.p_low == 0:
                self.p_low = self.press - 1
            if self.p_high == 0:
                self.p_high = self.press + 1
        self.long_name = 'Xenon'
        self.short_name = 'Xe'
        self.kappa = 1.67
        self.mass = 131.293
        self.liquid_enthalpy_tp = -1.2418
        self.vapour_enthalpy_tp = 95.165
        self.gas_id = 'C7440633'


class Ubertest:
    def __init__(self):
        self.long_name = ''
        self.kappa = 0
        pass

    def get_long_name(self):
        return self.long_name

    def get_kappa(self):
        return self.kappa


class Test1(Ubertest):
    def __init__(self):
        super().__init__()
        self.long_name = 'Helium'
        self.kappa = 23

    @staticmethod
    def get_kapp(x):
        return x.kappa


class Test2(Ubertest):
    def __init__(self):
        super().__init__()
        self.long_name = 'Hydrogen'
        self.kappa = 66


def get_k(x):
    return x.kappa


def main():
    a = H2(iso_type='IsoBar', temp=33, t_low=20, t_high=40, press=40)
    b = a.get_liquid_velocity_array()
    c = a.get_exact_supersonic_velocity_array()
    d = a.get_supercritical_velocity_array()
    e = a.get_simple_supersonic_velocity_array()
    T = a.get_temperature_column()
    plt.plot(T, b, 'b')
    plt.plot(T, c, 'y')
    plt.plot(T, d, 'r')
    plt.plot(T, e, 'g')
    plt.show()

# ------------------------

if __name__ == "__main__":
    main()


"""
Velocity calculation tool for internal target applications.

Mar 2015 by TwiNN

"""

# Needed libraries
from gas import *
from PyQt5 import QtWidgets
from velocity_window_ui import Ui_VelocityWindow


class VelocityWindow(QtWidgets.QWidget, Ui_VelocityWindow):
    def __init__(self):
        super().__init__()

        # Initial set up of the GUI

        self.setupUi(self)

        self.combo_box_decoding = {'Hydrogen': H2, 'Nitrogen': N2, 'Helium': He, 'Neon': Ne, 'Argon': Ar, 'Krypton': Kr, 'Xenon': Xe}
        self.selected_gas_species = self.combo_box_decoding[self.combo_box_select_gas.currentText()]()
        self.radio_select_data_type_isobaric.click()
        self.set_default_gas_values()

        # Connect signals
        self.connect_signals()

    def connect_signals(self):
        self.combo_box_select_gas.currentTextChanged.connect(self.set_default_gas_values)
        self.input_t_low.textChanged.connect(self.plot_velocity)

    def set_default_gas_values(self):
        self.selected_gas_species = self.combo_box_decoding[self.combo_box_select_gas.currentText()]
        self.input_temp.setText(str(self.selected_gas_species().temp))
        self.input_press.setText(str(self.selected_gas_species().press))
        self.input_t_low.setText(str(self.selected_gas_species().t_low))
        self.input_t_high.setText(str(self.selected_gas_species().t_high))
        self.input_p_high.setText(str(self.selected_gas_species().p_high))
        self.input_p_low.setText(str(self.selected_gas_species().p_low))

    def plot_velocity(self):
        gui_temp = float(self.input_temp.text())
        gui_press = float(self.input_press.text())
        gui_t_high = float(self.input_t_high.text())
        gui_t_low = float(self.input_t_low.text())
        gui_p_high = float(self.input_p_high.text())
        gui_p_low = float(self.input_p_low.text())
        data = self.selected_gas_species(iso_type='IsoBar', temp=gui_temp, t_low=gui_t_low, t_high=gui_t_high, press=gui_press)
        liquid = data.get_liquid_velocity_array()
        ss_exact = data.get_exact_supersonic_velocity_array()
        sc = data.get_supercritical_velocity_array()
        ss_simple = data.get_simple_supersonic_velocity_array()
        T = data.get_temperature_column()
        self.mpl_widget.canvas.ax.clear()
        self.mpl_widget.canvas.ax.plot(T, liquid, 'b')
        self.mpl_widget.canvas.ax.plot(T, ss_exact, 'y')
        self.mpl_widget.canvas.ax.plot(T, sc, 'r')
        self.mpl_widget.canvas.ax.plot(T, ss_simple, 'g')
        self.mpl_widget.canvas.draw()

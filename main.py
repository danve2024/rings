from visualization import Selection, Model, app
from measure import Measure
from units import *
import sys

# Runs the main application with given parameter values

"""
Automatic values:
    Asteroid:
        R: radius 1-21 km
        D: density 1.3-5.3 g/cm3
        A: semi-major axis 2.5-50 au
    Ring:
        d: density 0.01-0.03 g/cm3
        e: eccentricity 0-0.8
        i: inclination 0-90°
        a: semi-major axis (dependent on other parameters)
        m: mass (dependent on other parameters)
    Star:
        T: temperature
        log(g): logarithmic surface gravity
        λ: wavelength
! To change any parameters ranges use the window opened by the program or the Measure class
"""


defaults = {
            'radius': Measure(1, 21, km, label='R'),
            'density': Measure(1.3, 5.3, gcm3, label='D'),
            'ring_density': Measure(0.01, 0.03, gcm3, label='d'),
            'asteroid_sma': Measure(2.5, 50, au, label='A'),
            'sma': Measure(10, 1000, km, label='a'),  # Dependent slider
            'width': Measure(0.1, 20, km, label='w'), # Dependent slider
            'ring_mass': Measure(0.1, 10, kg, label='m'), # Dependent slider
            'eccentricity': Measure(0, 0.8, label='e'),
            'inclination': Measure(0, 90, deg, label='i'),
        }

# Graphics
selector = Selection(defaults)
app.exec()
parameters = selector.user_values

window = Model(parameters)
window.show()
sys.exit(app.exec())

from visualization import Selection, Model, app
from measure import Measure
from units import *
import sys


defaults = {
            'radius': Measure(1, 21, km, label='R'),
            'density': Measure(1.3, 5.3, gcm3, label='D'),
            'ring_density': Measure(0.01, 0.03, gcm3, label='d'),
            'asteroid_sma': Measure(2.5, 50, au, label='A'),
            'sma': Measure(10, 1000, km, label='?a'),  # Dependent slider
            'ring_mass': Measure(0.1, 10, kg, label='?m'),           # Dependent slider
            'eccentricity': Measure(0, 0.8, label='e'),
            'inclination': Measure(0, 90, deg, label='i'),
            'std_dev': Measure(0.1, 1, label='s')
        }

selector = Selection(defaults)
app.exec()
parameters = selector.user_values

window = Model(parameters)
window.show()
sys.exit(app.exec())

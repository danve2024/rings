import sys
import time

from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel
from PyQt6.QtCore import Qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import numpy as np

# Imports of user-defined modules
from measure import Measure
from units import *
from formulas import *
from space import Asteroid, Rings, Star, hill_sphere

# Initialize the star
star = Star(6.0 * mag, 0.8 * arcsec, None)

class DynamicPlot(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Asteroid with Rings Occultation Simulation")
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.sliders = {}
        self.slider_labels = {}


        # Default parameters
        self.defaults = {
            'radius': Measure(1, 21, km, label='R'),
            'density': Measure(1.3, 5.3, gcm3, label='D'),
            'ring_density': Measure(0.01, 0.03, gcm3, label='d'),
            'asteroid_sma': Measure(2.5, 50, au, label='A'),
            'sma': Measure(10, 1000, km, label='?a'),  # Dependent slider
            'ring_mass': Measure(0.1, 10, kg, label='?m'),           # Dependent slider
            'eccentricity': Measure(0, 0.8, label='e'),
            'inclination': Measure(0, 90, deg, label='i')
        }

        # Create sliders
        for key, measure in self.defaults.items():
            self.create_slider(key, measure)

        # Plot
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)

        # Time spent label
        self.time_label = QLabel("Time spent: 0.0 sec")
        self.layout.addWidget(self.time_label)

        # Occultation duration label
        self.eclipse_label = QLabel("Occultation duration: 0.0 sec")
        self.layout.addWidget(self.eclipse_label)

        self.update_plot()

    def create_slider(self, name, measure):
        container = QHBoxLayout()
        unit = self.get_unit_label(name)
        label = QLabel(f"{name.capitalize()} ({unit}):")
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(100)
        slider.setValue(50)
        slider.valueChanged.connect(self.update_plot)
        value_label = QLabel(str(measure.min))

        container.addWidget(label)
        container.addWidget(slider)
        container.addWidget(value_label)

        self.layout.addLayout(container)
        self.sliders[name] = (slider, measure)
        self.slider_labels[name] = value_label

    @staticmethod
    def get_unit_label(name: str) -> str:
        units = {
            'radius': 'km',
            'density': 'g/cm³',
            'ring_density': 'g/cm³',
            'asteroid_sma': 'AU',
            'sma': 'km',
            'ring_mass': 'kg',
            'eccentricity': '',
            'inclination': 'deg'
        }
        return units.get(name, '')

    def update_plot(self) -> None:
        start_time = time.time()

        params = {}
        for key, (slider, measure) in self.sliders.items():
            real_value = measure.slider(slider.value())
            params[key] = real_value
            self.slider_labels[key].setText(f"{real_value:.2f}")

        self.update_dependent_sliders(params)

        data, occultation_duration = self.calculate_data(**params)

        self.ax.clear()
        if len(data[1]) > 2:
            x, y = zip(*data[1])
            self.ax.plot(x, y)
        else:
            self.ax.plot([0, 1], [0, 0])
        self.ax.set_title("Simulation Result")
        self.canvas.draw()

        elapsed_time = time.time() - start_time
        self.time_label.setText(f"Time spent: {elapsed_time:.2f} sec")
        self.eclipse_label.setText(f"Occultation duration: {occultation_duration:.2f} sec")

    def update_dependent_sliders(self, params: dict) -> None:
        radius = params['radius'].set(km)
        density = params['density'].set(gcm3)
        ring_density = params['ring_density'].set(gcm3)
        asteroid_sma = params['asteroid_sma'].set(au)
        ring_mass = params['ring_mass'].set(kg)
        sma = params['sma'].set(km)

        V = volume(radius) # asteroid volume
        M = V * density    # asteroid mass
        a_min = max(roche_limit(radius, density, ring_density), radius) # semi-major axis minimum
        a_max = hill_sphere(asteroid_sma, M) # semi-major axis maximum
        self.defaults['sma'] = Measure(a_min, a_max, km, label='a')

        m_max = maximum_ring_mass(M, radius, sma)  # ring mass maximum
        m_min = 0.5 * m_max # ring mass minimum
        self.defaults['ring_mass'] = Measure(m_min, m_max, kg, label='m')

    @staticmethod
    def calculate_data(radius: Measure.Unit, density: Measure.Unit, ring_density: Measure.Unit, asteroid_sma: Measure.Unit, sma: Measure.Unit, ring_mass: Measure.Unit, eccentricity: Measure.Unit, inclination: Measure.Unit) -> tuple:
        radius = radius.set(km)
        density = density.set(gcm3)
        ring_density = ring_density.set(gcm3)
        asteroid_sma = asteroid_sma.set(au)
        sma = sma.set(km)
        ring_mass = ring_mass.set(kg)
        inclination = inclination.set(deg)

        V = volume(radius) # asteroid volume
        M = V * density    # asteroid mass
        rings = Rings(ring_density, sma, ring_mass, eccentricity, inclination)  # create rings
        asteroid = Asteroid(rings, radius, density, asteroid_sma, V, M)         # create asteroid
        data = star.occultation(asteroid)
        occultation_duration = data[0]
        return data, occultation_duration

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DynamicPlot()
    window.show()
    sys.exit(app.exec())

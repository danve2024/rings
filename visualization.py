import time
import sys
import traceback

from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel, QPushButton, QCheckBox, \
    QLineEdit, QDialog, QMessageBox
from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QPixmap
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import numpy as np

from measure import Measure
from units import *
from formulas import *
from space import Asteroid, Rings, Star, hill_sphere
from models import cover_animation

from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QGridLayout,
    QLabel, QLineEdit, QCheckBox, QPushButton, QDialog
)


class Selection(QWidget):
    def __init__(self, default_values):
        super().__init__()

        self.default_values = default_values
        self.user_values = {}

        self.checkboxes = {}
        self.mins = {}
        self.maxs = {}

        self.init_ui()

    def init_ui(self):
        layout = QGridLayout()  # Use a grid layout for column alignment

        # Set column headers
        layout.addWidget(QLabel("Parameter"), 0, 0)
        layout.addWidget(QLabel("Auto"), 0, 1)
        layout.addWidget(QLabel("Manual"), 0, 2)
        layout.addWidget(QLabel("Min"), 0, 3)
        layout.addWidget(QLabel("Max"), 0, 4)

        row = 1  # Start from the second row for parameters

        for param, measure in self.default_values.items():
            # Parameter Name Label
            param_label = QLabel(param)
            layout.addWidget(param_label, row, 0)

            # Auto Checkbox (Default)
            auto_cb = QCheckBox("Auto")
            auto_cb.setChecked(True)  # Start with Auto checked
            auto_cb.stateChanged.connect(lambda state, p=param: self.auto(state, p))
            layout.addWidget(auto_cb, row, 1)

            # Manual Checkbox
            manual_cb = QCheckBox("Manual")
            manual_cb.stateChanged.connect(lambda state, p=param: self.manual(state, p))
            layout.addWidget(manual_cb, row, 2)

            # Min and Max Input Fields
            min_input = QLineEdit(str(measure.min(measure.unit)))
            max_input = QLineEdit(str(measure.max(measure.unit)))
            min_input.setEnabled(False)
            max_input.setEnabled(False)
            layout.addWidget(min_input, row, 3)
            layout.addWidget(max_input, row, 4)

            # Store references
            self.checkboxes[param] = (auto_cb, manual_cb)
            self.mins[param] = min_input
            self.maxs[param] = max_input

            row += 1  # Move to the next row for the next parameter

        # Done Button
        done_btn = QPushButton("Done")
        done_btn.clicked.connect(self.click)  # Connect to close function
        layout.addWidget(done_btn, row, 0, 1, 5)  # Span across all columns

        self.setLayout(layout)
        self.setWindowTitle("Parameter Selector")
        self.show()

    def auto(self, state, param):
        if state == 2:  # Auto checked
            self.checkboxes[param][1].setChecked(False)  # Uncheck Manual
            self.mins[param].setEnabled(False)
            self.maxs[param].setEnabled(False)
        else:
            # Prevent unchecking both checkboxes
            if not self.checkboxes[param][1].isChecked():
                self.checkboxes[param][0].setChecked(True)  # Re-check Auto

    def manual(self, state, param):
        if state == 2:  # Manual checked
            self.checkboxes[param][0].setChecked(False)  # Uncheck Auto
            self.mins[param].setEnabled(True)
            self.maxs[param].setEnabled(True)
        else:
            # Prevent unchecking both checkboxes
            if not self.checkboxes[param][0].isChecked():
                self.checkboxes[param][1].setChecked(True)  # Re-check Manual

    def click(self):
        for param in self.default_values.keys():
            measure = self.default_values[param]

            if self.checkboxes[param][0].isChecked():  # Auto selected
                # Use existing measure values without changes
                self.user_values[param] = Measure(measure.min(measure.unit), measure.max(measure.unit), measure.unit,
                                                  measure.label)

            else:  # Manual selected
                min_val = float(self.mins[param].text())
                max_val = float(self.maxs[param].text())
                # Create a new Measure object with user-defined min/max but keep unit and label from defaults
                self.user_values[param] = Measure(min_val, max_val, measure.unit, measure.label)

        # Close the window after clicking Done
        self.close()


class Model(QWidget):
    def __init__(self, parameters: dict) -> None:
        super().__init__()
        self.setWindowTitle("Asteroid with Rings Occultation Simulation")
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.sliders = {}
        self.slider_labels = {}

        # Default parameters
        self.defaults = parameters

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

        # Show animation button
        self.animation_button = QPushButton("Simulate Occultation Animation")
        self.animation_button.clicked.connect(self.show_animation_window)
        self.layout.addWidget(self.animation_button, alignment=Qt.AlignmentFlag.AlignRight)

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
            'inclination': 'deg',
            'std_dev': ''
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

        data, occultation_duration, self.star, self.asteroid = self.calculate_data(**params)

        self.ax.clear()
        self.ax.set_title("Simulation Result")

        if len(data[1]) > 2:
            x, y = zip(*data[1])
            self.ax.plot(x, y)
        else:
            self.ax.plot([0, 1], [0, 0])

        self.ax.invert_yaxis()
        self.ax.set_xlabel("Phase")
        self.ax.set_ylabel("Magnitude Change")
        self.canvas.draw()

        elapsed_time = time.time() - start_time
        self.time_label.setText(f"Time spent: {elapsed_time:.2f} sec")
        self.eclipse_label.setText(f"Occultation duration: {occultation_duration:.2f} sec")

    def update_dependent_sliders(self, params: dict) -> None:
        radius = params['radius'].set(km)
        density = params['density'].set(gcm3)
        ring_density = params['ring_density'].set(gcm3)
        asteroid_sma = params['asteroid_sma'].set(au)
        sma = params['sma'].set(km)

        V = volume(radius)  # asteroid volume
        M = V * density  # asteroid mass
        a_min = max(roche_limit(radius, density, ring_density), radius)  # semi-major axis minimum
        a_max = hill_sphere(asteroid_sma, M)  # semi-major axis maximum
        self.defaults['sma'].update(a_min, a_max)

        m_max = maximum_ring_mass(M, radius, sma)  # ring mass maximum
        m_min = 0.5 * m_max  # ring mass minimum
        self.defaults['ring_mass'].update(m_min, m_max)

    @staticmethod
    def calculate_data(radius: Measure.Unit, density: Measure.Unit, ring_density: Measure.Unit,
                       asteroid_sma: Measure.Unit, sma: Measure.Unit, ring_mass: Measure.Unit,
                       eccentricity: Measure.Unit, inclination: Measure.Unit, std_dev: Measure.Unit) -> tuple:
        # Star initialization
        magnitude = 6.0 * mag
        angular_size = 0.8 * arcsec
        star = Star(magnitude, angular_size, std_dev)

        radius = radius.set(km)
        density = density.set(gcm3)
        ring_density = ring_density.set(gcm3)
        asteroid_sma = asteroid_sma.set(au)
        sma = sma.set(km)
        ring_mass = ring_mass.set(kg)
        inclination = inclination.set(deg)

        V = volume(radius)  # asteroid volume
        M = V * density  # asteroid mass
        rings = Rings(ring_density, sma, ring_mass, eccentricity, inclination)  # create rings
        asteroid = Asteroid(rings, radius, density, asteroid_sma, V, M)  # create asteroid
        data = star.occultation(asteroid)
        occultation_duration = data[0]

        return data, occultation_duration, star, asteroid

    def show_animation_window(self):
        """Open the animation window with error handling."""
        try:
            if hasattr(self, 'star') and hasattr(self, 'asteroid'):
                animation_window = AnimationWindow(self.star, self.asteroid, self)
                animation_window.exec()
            else:
                QMessageBox.warning(self, "Warning", "Please update the plot before running the animation.")
        except Exception as e:
            print("Error during animation:", e)
            traceback.print_exc()


class AnimationWindow(QDialog):
    def __init__(self, star, asteroid, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Asteroid Ring Occultation Simulation")
        self.setModal(True)
        self.resize(600, 600)

        self.layout = QVBoxLayout()
        self.label = QLabel()
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.layout.addWidget(self.label)
        self.setLayout(self.layout)

        self.frames = cover_animation(star.model, asteroid.model)  # Generate frames
        self.current_frame_index = 0

        self.timer = QTimer()
        self.timer.timeout.connect(self.update_frame)
        self.timer.start(66)  # ~15 FPS

    def update_frame(self):
        """Update the animation frame safely."""
        if not self.frames:
            return  # Avoid update if frames are missing

        frame = self.frames[self.current_frame_index]
        if frame.isNull():
            return  # Skip invalid frames

        pixmap = QPixmap.fromImage(frame)
        self.label.setPixmap(pixmap)
        self.current_frame_index = (self.current_frame_index + 1) % len(self.frames)


app = QApplication(sys.argv)

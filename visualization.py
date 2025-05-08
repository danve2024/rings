import time
import sys
import traceback

from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel, QPushButton, QCheckBox, \
    QLineEdit, QDialog, QMessageBox, QFileDialog
from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QPixmap
import matplotlib.pyplot as plt
from astropy.units import temperature
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import numpy as np

from measure import Measure
from units import *
from formulas import *
from space import Asteroid, Rings, Star, hill_sphere
from models import cover_animation
from observations import Observations

from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QGridLayout,
    QLabel, QLineEdit, QCheckBox, QPushButton, QDialog
)


class Selection(QWidget):
    """
    A widget to select and configure parameters for the simulation.
    """
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
            if self.get_unit_label(param):
                param_label = QLabel(f'{param.replace("_", " ").replace("sma", "semi-major axis")} ({self.get_unit_label(param)})')
            else:
                param_label = QLabel(param.replace("_", " ").replace("sma", "semi-major axis"))
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

    @staticmethod
    def get_unit_label(name: str) -> str:
        """Get the unit label for a given parameter"""
        units = {
            'radius': 'km',
            'density': 'g/cm³',
            'ring_density': 'g/cm³',
            'asteroid_sma': 'au',
            'sma': 'km',
            'width': 'km',
            'ring_mass': 'kg',
            'eccentricity': '',
            'inclination': 'deg',
        }
        return units.get(name, '')

    def auto(self, state, param):
        """Auto select parameter"""
        if state == 2:  # Auto checked
            self.checkboxes[param][1].setChecked(False)  # Uncheck Manual
            self.mins[param].setEnabled(False)
            self.maxs[param].setEnabled(False)
        else:
            # Prevent unchecking both checkboxes
            if not self.checkboxes[param][1].isChecked():
                self.checkboxes[param][0].setChecked(True)  # Re-check Auto

    def manual(self, state, param):
        """Manually set the parameter values."""
        if state == 2:  # Manual checked
            self.checkboxes[param][0].setChecked(False)  # Uncheck Auto
            self.mins[param].setEnabled(True)
            self.maxs[param].setEnabled(True)
        else:
            # Prevent unchecking both checkboxes
            if not self.checkboxes[param][0].isChecked():
                self.checkboxes[param][1].setChecked(True)  # Re-check Manual

    def click(self):
        """Close the window and store user-defined parameter values."""
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

class LoadFile(QWidget):
    """
    A widget to load an observations file for comparing observations with the model.
    """
    def __init__(self):
        super().__init__()
        self.button = None
        self.proceed = None
        self.filename = None
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Load observations file")
        self.setGeometry(500, 500, 400, 300)
        self.button = QPushButton("Open File", self)
        self.button.clicked.connect(self.openFileDialog)
        self.button.setGeometry(150, 130, 100, 30)

        self.proceed = QPushButton("No File", self)
        self.proceed.clicked.connect(self.close)
        self.proceed.setGeometry(150, 170, 100, 30)
        self.show()

    def openFileDialog(self):
        """
        A function for opening the file dialog and pre-processing the selected file.
        """
        file_dialog = QFileDialog(self)
        file_dialog.setNameFilter("CSV Files (*.csv)") # Search for CSV files only
        file_dialog.setWindowTitle("Select File")
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)
        file_dialog.setViewMode(QFileDialog.ViewMode.Detail)

        if file_dialog.exec():
            selected_files = file_dialog.selectedFiles()
            self.filename = selected_files[0]
            self.close()


class Model(QWidget):
    """
    A widget to display the simulation model and control parameters sliders.
    """
    def __init__(self, parameters: dict, observations_file: str = None) -> None:
        super().__init__()

        if observations_file is None:
            self.observations = None
        else:
            self.observations = Observations(observations_file) # Working with the observation data if it is available

        self.magnitude_shift = None
        self.magnitude_calibrating = None
        self.phase_shift = None
        self.magnitude_shift_label = None
        self.phase_shift_label = None
        self.magnitude_calibrating_label = None

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

        # Sliders for working with the observed data
        if self.observations is not None:
            self.magnitude_shift_slider()
            self.magnitude_calibrating_slider()
            self.phase_shift_slider()

        # Show animation button
        self.animation_button = QPushButton("Simulate Occultation Animation")
        self.animation_button.clicked.connect(self.show_animation_window)
        self.layout.addWidget(self.animation_button, alignment=Qt.AlignmentFlag.AlignRight)

        for i in range(3):
            self.update_plot()

    def create_slider(self, name, measure):
        """Create a slider"""
        container = QHBoxLayout()
        unit = self.get_unit_label(name)
        if unit:
            label = QLabel(f"{name.replace('_', ' ').replace('sma', 'semi-major axis').capitalize()} ({unit}):")
        else:
            label = QLabel(name.replace('_', ' ').replace('sma', 'semi-major axis').capitalize())
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

    def magnitude_shift_slider(self):
        """
        A slider for making major shifts (0-12 mag) of the magnitude of the by a specific value selected by the slider. It is needed to work with the magnitude change (produced by the model) instead of the magnitude itself (observed).
        """
        container = QHBoxLayout()
        label = QLabel('Observations magnitude shift')
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setMinimum(-20)
        slider.setMaximum(20)
        slider.setValue(0)
        slider.valueChanged.connect(self.update_plot)
        value_label = QLabel('0')

        container.addWidget(label)
        container.addWidget(slider)
        container.addWidget(value_label)

        self.layout.addLayout(container)

        self.magnitude_shift = slider
        self.magnitude_shift_label = value_label

    def magnitude_calibrating_slider(self):
        """
            A slider for making minor shifts (0-1 mag) of the magnitude of the by a specific value selected by the slider. It is needed to work with the magnitude change (produced by the model) instead of the magnitude itself (observed).
        """
        container = QHBoxLayout()
        label = QLabel('Observations magnitude calibrating')
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setMinimum(-1000)
        slider.setMaximum(1000)
        slider.setValue(0)
        slider.valueChanged.connect(self.update_plot)
        value_label = QLabel('0')

        container.addWidget(label)
        container.addWidget(slider)
        container.addWidget(value_label)

        self.layout.addLayout(container)

        self.magnitude_calibrating = slider
        self.magnitude_calibrating_label = value_label

    def phase_shift_slider(self):
        """
        A slider for selecting the moment of the observations to compare with the model. The observation data is cropped by the x-axis and moved for the specific time span selected.
        """
        container = QHBoxLayout()
        label = QLabel('Observations phase shift')
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(100)
        slider.setValue(0)
        slider.valueChanged.connect(self.update_plot)
        value_label = QLabel('0')

        container.addWidget(label)
        container.addWidget(slider)
        container.addWidget(value_label)

        self.layout.addLayout(container)

        self.phase_shift = slider
        self.phase_shift_label = value_label

    @staticmethod
    def get_unit_label(name: str) -> str:
        """Get the unit label for a given parameter"""
        units = {
            'radius': 'km',
            'density': 'g/cm³',
            'ring_density': 'g/cm³',
            'asteroid_sma': 'au',
            'sma': 'km',
            'width': 'km',
            'ring_mass': 'kg',
            'eccentricity': '',
            'inclination': 'deg',
        }
        return units.get(name, '')

    def update_plot(self) -> None:
        """Update the plot based on the current slider values"""
        start_time = time.time()

        params = {}
        for key, (slider, measure) in self.sliders.items():
            real_value = measure.slider(slider.value())
            params[key] = real_value
            self.slider_labels[key].setText(f"{real_value:.2f}")

        self.update_dependent_sliders(params)

        data, occultation_duration, self.star, self.asteroid = self.calculate_data(**params)

        if self.observations is not None: # shifting (magnitude to magnitude change and time to phase) and normalizing (time to phase) the observation data
            phase_shift = self.phase_shift.value() / 100
            magnitude_shift = self.magnitude_shift.value()
            magnitude_calibrating = self.magnitude_calibrating.value() / 1000
            self.phase_shift_label.setText(str(phase_shift))
            self.magnitude_shift_label.setText(str(magnitude_shift))
            self.magnitude_calibrating_label.setText(str(magnitude_calibrating))

            self.observations.normalize_and_shift(occultation_duration, phase_shift, magnitude_shift + magnitude_calibrating)

        self.ax.clear()
        self.ax.set_title("Simulation Result")

        # Plotting the observations data with the selected shifts
        if self.observations is None:
            if len(data[1]) > 2:
                x, y = zip(*data[1])
                self.ax.plot(x, y, label='model')
            else:
                self.ax.plot([0, 1], [0, 0], label='model')
        else:
            if len(data[1]) > 2:
                x, y = zip(*data[1])
                self.ax.plot(x, y, label='model')
            else:
                self.ax.plot([0, 1], [0, 0], label='model')
            if len(self.observations.data) > 1:
                xo, yo = zip(*self.observations.data)
                self.ax.plot(xo, yo, 'g', label='observations')
            else:
                pass

        self.ax.invert_yaxis()
        self.ax.set_xlabel("Phase")
        self.ax.set_ylabel("Magnitude Change")
        self.ax.legend()
        self.canvas.draw()

        elapsed_time = time.time() - start_time
        self.time_label.setText(f"Time spent: {elapsed_time:.2f} sec")
        self.eclipse_label.setText(f"Occultation duration: {occultation_duration:.2f} sec")

    def update_dependent_sliders(self, params: dict) -> None:
        """Update the dependent sliders (ring mass and semi-major axis)"""
        radius = params['radius'].set(km)
        density = params['density'].set(gcm3)
        ring_density = params['ring_density'].set(gcm3)
        asteroid_sma = params['asteroid_sma'].set(au)

        V = volume(radius)  # asteroid volume
        M = V * density  # asteroid mass
        a_min = max(roche_limit(radius, density, ring_density), radius)  # semi-major axis minimum
        # a_max = a_min * 10
        a_max = hill_sphere(asteroid_sma, M)
        self.defaults['sma'].update(a_min / km, a_max / km)
        sma = params['sma'].set(km)

        w_max = min(a_max - a_min, sma / 2) # ring width maximum
        self.defaults['width'].update(10, w_max / km)

        m_max = maximum_ring_mass(M, radius, sma)  # ring mass maximum
        m_min = 0.5 * m_max  # ring mass minimum
        self.defaults['ring_mass'].update(m_min / kg, m_max / kg)

    @staticmethod
    def calculate_data(radius: Measure.Unit, density: Measure.Unit, ring_density: Measure.Unit,
                       asteroid_sma: Measure.Unit, sma: Measure.Unit, width: Measure.Unit, ring_mass: Measure.Unit,
                       eccentricity: Measure.Unit, inclination: Measure.Unit) -> tuple:
        """Calculate the simulation data"""
        # Parameters
        radius = radius.set(km)
        density = density.set(gcm3)
        ring_density = ring_density.set(gcm3)
        asteroid_sma = asteroid_sma.set(au)
        sma = sma.set(km)
        width = width.set(km)
        ring_mass = ring_mass.set(kg)
        inclination = inclination.set(deg)

        # Create asteroid and calculate simulation data
        V = volume(radius)  # asteroid volume
        M = V * density  # asteroid mass
        rings = Rings(ring_density, sma, width, ring_mass, eccentricity, inclination)  # create rings
        asteroid = Asteroid(rings, radius, density, asteroid_sma, V, M)  # create asteroid

        # Star initialization
        angular_size = 0.8 * arcsec
        apsis = asteroid.rings.angular_sma * (1 + asteroid.rings.eccentricity) # apsis: Q = a(1+e)
        star = Star(2 * apsis, angular_size)

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
    """Occultation animation window (opened using the "Show occultation animation" window)"""
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
        self.timer.start(33)  # ~30 FPS

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

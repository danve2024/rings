from units import *
from formulas import *
from measure import Measure
from math import pi, exp, cos, sqrt
from typing import Union
from models import disk, elliptical_ring, normalize, cover, star_model #, gaussian
import numpy as np

# Functions that require other Solar System objects

"""
Parameters used:
    Asteroid:
        R - radius
        D - density
        M - mass
        A - semi-major axis
        V - volume
    Ring:
        r - radius
        d - density
        m - mass
        a - semi-major axis
        w - width
        e - eccentricity
"""

# Parameters for a nonlinear star model
# T -> log(g) -> λ
darkening_coefficients = {5500:
                              {3:
                                   {3437: [0.985, 0.994],
                                    4212: [0.989, 0.994],
                                    4687: [0.989, 0.997],
                                    5475: [0.985, 0.997],
                                    6975: [0.986, 0.999]},
                               4:
                                   {3437: [0.985, 0.994],
                                    4212: [0.989, 0.994],
                                    4687: [0.989, 0.997],
                                    5475: [0.985, 0.997],
                                    6975: [0.986, 0.999]},
                               },
                          6000:
                              {4:
                                   {3437: [0.972, 0.994],
                                    4212: [0.974, 0.993],
                                    4687: [0.977, 0.995],
                                    5475: [0.978, 0.996],
                                    6975: [0.980 , 0.997]},
                               },
                          7000:
                              {4:
                                   {3437: [0.980, 0.998],
                                    4212: [0.966, 0.994],
                                    4687: [0.967, 0.994],
                                    5475: [0.975, 0.997],
                                    6975: [0.978, 0.997]},
                               },
                          8000:
                              {3:
                                   {3437: [0.988, 1.001],
                                    4212: [0.970, 0.998],
                                    4687: [0.971, 0.998],
                                    5475: [0.972, 0.998],
                                    6975: [0.979, 0.999]},
                               4:
                                   {3437: [0.986, 0.999],
                                    4212: [0.956, 0.996],
                                    4687: [0.960, 0.996],
                                    5475: [0.967, 0.997],
                                    6975: [0.977, 0.998]},
                               },
                          10000:
                              {3:
                                   {3437: [0.988, 1.000],
                                    4212: [0.971, 0.999],
                                    4687: [0.973, 0.999],
                                    5475: [0.976, 0.999],
                                    6975: [0.982, 1.000]},
                               4:
                                   {3437: [0.988, 1.000],
                                    4212: [0.969, 0.999],
                                    4687: [0.972, 1.000],
                                    5475: [0.976, 0.999],
                                    6975: [0.982, 1.000]},
                               },
                          15000:
                              {3:
                                   {3437: [0.984, 1.000],
                                    4212: [0.976, 1.000],
                                    4687: [0.978, 1.000],
                                    5475: [0.979, 1.000],
                                    6975: [0.983, 1.001]},
                               4:
                                   {3437: [0.985, 1.001],
                                    4212: [0.975, 1.000],
                                    4687: [0.977, 1.000],
                                    5475: [0.981, 1.000],
                                    6975: [0.985, 1.000]},
                               },
                          20000:
                              {3:
                                   {3437: [0.981, 1.000],
                                    4212: [0.976, 0.999],
                                    4687: [0.978, 1.000],
                                    5475: [0.979, 1.000],
                                    6975: [0.982, 1.001]},
                               4:
                                   {3437: [0.984, 1.000],
                                    4212: [0.977, 1.000],
                                    4687: [0.979, 1.000],
                                    5475: [0.982, 1.000],
                                    6975: [0.986, 1.001]},
                               }
                          }

def hill_sphere(asteroid_semi_major_axis: Union[float, Measure.Unit], asteroid_mass: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
    Computes the hill sphere radius of an asteroid.
    a_hill = max(a_hill_Sun(A, M), a_hill_Jupiter(A, M), a_hill_Neptune(A, M))

    :param asteroid_semi_major_axis: A
    :param asteroid_mass: M
    :return: a_hill
    """
    return max(sun.hill_sphere(asteroid_semi_major_axis, asteroid_mass), jupiter.hill_sphere(asteroid_semi_major_axis, asteroid_mass), neptune.hill_sphere(asteroid_semi_major_axis, asteroid_mass))

def illuminance(magnitude: Union[float, Measure.Unit]) -> Union[float, Measure]:
    """
    Converts apparent magnitude to illuminance.
    E = E_Sun * 10^(0.4(m_Sun - m0))

    :param magnitude: m0 - body apparent magnitude
    :return: E - illuminance
    """
    return sun.illuminance * 10 ** (0.4 * (sun.magnitude - magnitude))

def synodic_period(sma: Union[float, Measure.Unit]) -> Union[float, Measure]:
    """
    Computes the synodic period of a celestial body.
    S = 1/(|1/T(A) - 1/T_Earth|)
    T(A) - period of the body calculated using the 3rd Kepler's law (formulas.py/period)

    :param sma: A
    :return: S - synodic period
    """
    return 1/math.fabs(1/period(sma) - 1/earth.period)


# Classes of the main celestial bodies used for the model and their basic objects

class SSBody:
    """
    A base class for celestial bodies.
    """
    def __init__(self, mass: float, sma: float):
        self.mass = mass # m_i
        self.sma = sma # a_i
        self.period = self.calculate_period() # t_i
    def hill_sphere(self, a_sma: float, a_mass: float) -> float:
        """
        Computes the hill sphere radius of the celestial body.
        a_hill_i = A∛(M/3(M+m_i))

        :param a_sma: A
        :param a_mass: M
        :return: a_hill_i - hill sphere radius
        """
        return a_sma * math.pow(a_mass / (3 * (a_mass + self.mass)), 1/3)
    def calculate_period(self):
        """
        Computes the orbital period of the celestial body using the 3rd Kepler's law (formulas.py/period).

        :return: t_i(a_i) - orbital period
        """
        return period(self.sma)

    def __str__(self):
        return f"ssbody(m: {self.mass/kg}kg, a: {self.sma/au}au)"


# Planet constants
earth = SSBody(6 * 10 ** 24 * kg, 1 * au)
neptune = SSBody(10 ** 26 * kg, 30 * au)
jupiter = SSBody(2 * 10 ** 27 * kg, 5.2 * au)


class Sun(SSBody):
    """
    A class for the Sun.
    """
    def __init__(self, mass: float, sma: float, luminosity: float, magnitude: float, planet: SSBody=earth):
        super().__init__(mass, sma)
        self.luminosity = luminosity # luminosity
        self.magnitude = magnitude # magnitude
        self.planet = planet # planet we observe the Sun from
        self.illuminance = self.luminosity / (4 * pi * self.planet.sma**2) # illuminance

    def __str__(self):
        return f"sun(L: {self.illuminance}W, m: {self.magnitude}, E: {self.illuminance}W/m²)"

# The Sun constant
sun = Sun(2 * 10 ** 30 * kg, 0, 3.8 * 10 ** 26 * W, -26.74 * mag)


class Rings:
    """
    A class for the rings model.
    """
    def __init__(self, density: Measure.Unit, sma: Measure.Unit, width: Measure.Unit, mass: Measure.Unit, eccentricity: Measure.Unit, inclination: Measure.Unit):
        # Ring parameters
        self.density = density # density
        self.sma = sma # semi-major axis
        self.mass = mass # mass
        self.eccentricity = eccentricity # eccentricity
        self.inclination = inclination # inclination
        self.volume = self.mass / self.density # volume
        self.width = width # width

        # For calculating the ring transparency
        self.mass_specific_absorption_coefficient = 2.3 * 10 ** (-3) * (m**2/g) # mass-specific absorption coefficient for silicate dust with quartz dominating
        self.absorption_coefficient = self.mass_specific_absorption_coefficient * self.density # absorption coefficient for silicate
        self.absorption = exp(-self.absorption_coefficient) # light absorption

        # Parameters defined in self.init()
        self.angular_sma = None # angular semi-major axis
        self.angular_width = None # angular width
        self.size = None  # matrix size
        self.model = None  # model of the ring

    def init(self, sma: Union[float, Measure.Unit], size: Union[float, Measure.Unit]) -> None:
        """
        Initializes the ring numpy mask array model using models.py/elliptical_ring

        :param sma: a
        :param size: matrix size (used for matrices concatenation)
        """
        self.angular_sma = to_angle(self.sma, sma)  # angular size of semi-major axis
        self.angular_width = to_angle(self.width, sma)  # angular width
        self.size = size # matrix size
        self.model = elliptical_ring(self.size, to_pixels(self.angular_sma), self.eccentricity, to_pixels(self.angular_width), self.inclination, self.absorption)

    def adjust(self, size) -> None:
        """
        Adjusts the ring model size.

        :param size: matrix size (used for matrices concatenation)
        """
        self.model = elliptical_ring(size, to_pixels(self.angular_sma), self.eccentricity, to_pixels(self.angular_width), self.inclination, self.absorption)

    def __str__(self) -> str:
        return f'rings(d:{self.density(gcm3)}g/cm3, a:{self.sma(km)}km, w: {self.width(km)}km, m:{self.mass(kg)}kg, e:{self.eccentricity}, i:{self.inclination(deg)}°)'


class Asteroid:
    def __init__(self, rings: Rings, radius: Measure.Unit, density: Measure.Unit, sma: Measure.Unit, vol: Measure.Unit, mass: Measure.Unit):
        """
        Sets asteroid parameters and creates its numpy mask array representation with its rings.
        """
        # Asteroid parameters
        self.rings = rings # ring
        self.radius = radius # radius
        self.density = density # density
        self.sma = sma # semi-major axis
        self.volume = vol # volume
        self.mass = mass # mass
        self.rings = rings # rings
        self.angular_diameter = angular_diameter(self.radius, self.sma) # angular diameter
        self.synodic_period = synodic_period(self.sma)  # synodic period
        self.angular_velocity = angular_velocity(self.sma, self.synodic_period) # angular velocity

        # Create asteroid model with its rings
        self.rings.init(self.sma, to_pixels(self.angular_diameter))
        apsis = self.rings.angular_sma * (1 + self.rings.eccentricity) # apsis: Q = a(1+e)
        crop_factor = to_pixels(max(self.angular_diameter, 2 * apsis))
        self.disk = disk(to_pixels(self.angular_diameter/2), crop_factor) # disk
        self.rings.adjust(crop_factor)
        self.adjust(crop_factor)
        self.model = normalize(self.disk + self.rings.model)

    def adjust(self, size: Union[float, Measure.Unit]) -> None:
        """
        Adjusts the asteroid model size and rings model size.

        :param size: matrix size (used for matrices concatenation)
        """
        self.disk = disk(to_pixels(self.angular_diameter/2), size)
        self.rings.adjust(size)
        self.model = normalize(self.disk + self.rings.model)

    def __str__(self):
        return f'asteroid(R:{self.radius(km)}km, D:{self.density(gcm3)}g/cm3, A:{self.sma(au)}, V:{self.volume(m3)}m3, M:{self.mass(kg)}kg) + ' + str(self.rings)

class Star:
    def __init__(self, cropped: Union[float, Measure.Unit], angular_radius: Union[float, Measure.Unit], temperature: Union[float, Measure.Unit] = 8000, log_g: Union[float, Measure.Unit] = 4, wavelength: Union[float, Measure.Unit] = 4687) -> None: # std_dev -> temperature, log_g
        # Basic star parameters (will be replaced)
        self.radius = angular_radius # angular radius
        self.diameter = 2 * self.radius # angular diameter
        self.area = angular_area(self.radius) # angular area
        self.temperature = temperature # temperature
        self.log_g = log_g # log(g)
        self.cropped = cropped # width in pixels (should be equal to the asteroid diameter)
        self.wavelength = wavelength # observation wavelength

        # Intensity calculation
        # Derivation:
        # I(μ) = I(1)[1-c(1-μ)-d(1-√μ)]

        # I(μ) - pixel brightness in relative units (0 - min, 1 - max)
        # μ - cosine of the angle by the emergent radiation and the direct perpendicular to the stellar surface

        # Full brightness: I = ∫[0, 1]2πI(μ)dμ
        # I = π(15-15c-15d+12d+10с)/15 = π(15-5c-3d)/15
        try:
            c, d = darkening_coefficients[self.temperature][self.log_g][self.wavelength]
        except IndexError:
            raise ValueError(f'Wrong star parameter value selected (temperature [{self.temperature/K}K], log_g [{self.log_g}] or wavelength [{self.wavelength}Å]).\nTo see available parameter values check space.py/darkening_coefficients.')

        self.intensity = np.sum(star_model([round(to_pixels(self.radius*2)), round(to_pixels(self.radius*2))], [c, d]))

        # Creating the model
        # self.model = gaussian(to_pixels(self.radius*2), std_dev=std_dev) # disk
        self.model = star_model([round(to_pixels(self.cropped)), round(to_pixels(self.radius*2))],
                                [c, d]) # -> Will use nonlinear limb darkening approximation


    def __str__(self):
        return f"star(cropped: {self.cropped/arcsec}'', angular_radius: {self.radius/arcsec}'', temperature: {self.temperature/K}K, log_g: {self.log_g}, wavelength: {self.wavelength}Å)"

    '''
    def match(self, other: Asteroid) -> None:
        """
        Adjusts the star model size to match the size of the other asteroid.

        :param Asteroid other: the covering asteroid
        """
        size = to_pixels(max(2 * self.radius, other.angular_diameter))
        other.adjust(size)
    '''

    def occultation(self, asteroid: Asteroid):
        """
        Calculates the occultation time and the lightcurve of the magnitude change of the star being observed by the asteroid.

        :param Asteroid asteroid: the covering asteroid
        :return: [t, Δm(Φ)] - occultation time and its lightcurve
        """

        time = (asteroid.angular_diameter + self.diameter) / asteroid.angular_velocity
        # self.match(asteroid)
        return [time, cover(self.model, asteroid.model, self.intensity)]

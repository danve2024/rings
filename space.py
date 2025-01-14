from units import *
from formulas import *
from measure import Measure
from math import pi, exp
from typing import Union
from models import disk, gaussian, elliptical_ring, normalize, cover

# Functions that require other Solar System objects

def hill_sphere(asteroid_semi_major_axis: Union[float, Measure.Unit], asteroid_mass: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return max(sun.hill_sphere(asteroid_semi_major_axis, asteroid_mass), jupiter.hill_sphere(asteroid_semi_major_axis, asteroid_mass), neptune.hill_sphere(asteroid_semi_major_axis, asteroid_mass))

def illuminance(magnitude: Union[float, Measure.Unit]) -> Union[float, Measure]:
    return sun.illuminance * 10 ** (0.4 * (sun.magnitude - magnitude))

def synodic_period(sma: Union[float, Measure.Unit]) -> Union[float, Measure]:
    return 1/math.fabs(1/period(sma) - 1/earth.period)


# Classes and their basic objects

class SSBody:
    def __init__(self, mass: float, sma: float):
        self.mass = mass
        self.sma = sma
        self.period = self.calculate_period()
    def hill_sphere(self, a_sma: float, a_mass: float) -> float:
        return a_sma * math.pow(a_mass / (3 * (a_mass + self.mass)), 1/3)
    def calculate_period(self):
        return period(self.sma)


earth = SSBody(6 * 10 ** 24 * kg, 1 * au)
neptune = SSBody(10 ** 26 * kg, 30 * au)
jupiter = SSBody(2 * 10 ** 27 * kg, 5.2 * au)


class Sun(SSBody):
    def __init__(self, mass: float, sma: float, luminosity: float, magnitude: float, planet: SSBody=earth):
        super().__init__(mass, sma)
        self.luminosity = luminosity # luminosity
        self.magnitude = magnitude # magnitude
        self.planet = planet # planet we observe the Sun from
        self.illuminance = self.luminosity / (4 * pi * self.planet.sma**2) # illuminance

sun = Sun(2 * 10 ** 30 * kg, 0, 3.8 * 10 ** 26 * W, -26.74 * mag)


class Rings:
    def __init__(self, density: Measure.Unit, sma: Measure.Unit, mass: Measure.Unit, eccentricity: Measure.Unit, inclination: Measure.Unit):
        # Ring parameters
        self.density = density # density
        self.sma = sma # semi-major axis
        self.mass = mass # mass
        self.eccentricity = eccentricity # eccentricity
        self.inclination = inclination # inclination
        self.volume = self.mass / self.density # volume
        self.width = ring_width(self.volume, self.sma, self.eccentricity)

        self.mass_specific_absorption_coefficient = 2.3 * 10 ** (-3) * (m**2/g) # mass-specific absorption coefficient for silicate dust with quartz dominating
        self.absorption_coefficient = self.mass_specific_absorption_coefficient * self.density # absorption coefficient for silicate
        self.absorption = exp(-self.absorption_coefficient) # light absorption

        self.angular_sma = None # angular semi-major axis
        self.angular_width = None # angular width
        self.size = None  # matrix size
        self.model = None  # model of the ring

    def init(self, sma: Union[float, Measure.Unit], size: Union[float, Measure.Unit]) -> None:
        self.angular_sma = to_angle(self.sma, sma)  # angular size of semi-major axis
        self.angular_width = to_angle(self.width, sma)  # angular width
        self.size = size # matrix size
        self.model = elliptical_ring(self.size, to_pixels(self.angular_sma), self.eccentricity, to_pixels(self.angular_width), self.inclination, self.absorption)

    def adjust(self, size) -> None:
        self.model = elliptical_ring(size, to_pixels(self.angular_sma), self.eccentricity, to_pixels(self.width), self.inclination, self.absorption)

    def __str__(self) -> str:
        return f'd:{self.density(gcm3)}g/cm3 a:{self.sma(km)}km m:{self.mass(kg)}kg e:{self.eccentricity} i:{self.inclination(deg)}°\n'


class Asteroid:
    def __init__(self, rings: Rings, radius: Measure.Unit, density: Measure.Unit, sma: Measure.Unit, vol: Measure.Unit, mass: Measure.Unit):
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
        self.disk = disk(to_pixels(self.angular_diameter/2)) # disk
        self.rings.init(self.sma, self.disk.shape[0])
        self.model = normalize(self.disk + self.rings.model)

    def adjust(self, size: Union[float, Measure.Unit]) -> None:
        self.disk = disk(to_pixels(self.angular_diameter/2), size)
        self.rings.adjust(size)
        self.model = normalize(self.disk + self.rings.model)

    def __str__(self):
        return f'R:{self.radius(km)}km D:{self.density(gcm3)}g/cm3 A:{self.sma(au)} V:{self.volume(m3)}m3 M:{self.mass(kg)}kg' + ' ' + str(self.rings)


class Star:
    def __init__(self, magnitude: float, angular_radius: float, std_dev: Union[float, Measure.Unit]) -> None:
        # Basic star parameters
        self.magnitude = magnitude # magnitude
        self.radius = angular_radius # angular radius
        self.diameter = 2 * self.radius # angular diameter
        self.area = angular_area(self.radius) # angular area
        self.surface_brightness = self.magnitude / self.area # surface brightness
        self.illuminance = illuminance(magnitude) # illuminance
        self.brightness = self.illuminance / self.area # brightness
        self.model = gaussian(to_pixels(self.radius*2), std_dev=std_dev) # disk

    def adjust(self, size: Union[float, Measure.Unit]) -> None:
        self.model = gaussian(size)

    def match(self, other: Asteroid) -> None:
        size = to_pixels(max(2 * self.radius, other.angular_diameter))
        self.adjust(size)
        other.adjust(size)

    def occultation(self, asteroid: Asteroid):
        time = (asteroid.angular_diameter + self.diameter) / asteroid.angular_velocity
        self.match(asteroid)
        return [time, cover(self.model, asteroid.model)]



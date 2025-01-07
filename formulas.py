import math
from typing import Union
from measure import Measure
from units import au, yrs

def volume(radius: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return 4/3 * math.pi * radius**3

def roche_limit(radius: Union[float, Measure.Unit], asteroid_density: Union[float, Measure.Unit], ring_density: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return radius * (2 * asteroid_density/ring_density) ** (1/3)

def maximum_ring_mass(asteroid_mass: Union[float, Measure.Unit], asteroid_radius: Union[float, Measure.Unit], ring_sma: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return (asteroid_mass * asteroid_radius)/(2 * ring_sma)

def angular_area(angular_radius: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return math.pi * angular_radius**2

def angular_diameter(radius: Union[float, Measure.Unit], distance: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return 2 * math.degrees(radius / distance)

def to_angle(size: Union[float, Measure.Unit], distance: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return 2 * math.degrees(size/2 / distance)

def period(sma: Union[float, Measure.Unit]) -> Union[float, Measure]:
    return math.sqrt((sma/au)**3) * yrs

def velocity(sma: Union[float, Measure.Unit], time: Union[float, Measure.Unit]) -> Union[float, Measure]:
    return (2 * math.pi * sma)/time

def angular_velocity(sma: Union[float, Measure.Unit], time: Union[float, Measure.Unit]) -> Union[float, Measure]:
    return math.degrees(velocity(sma, time) / sma)

def format_data(data: list) -> list:
    ans = []
    initial_illuminance = data[0]
    for i in range(len(data)):
        illuminance = data[i]
        magnitude_change = 2.5 * math.log10(initial_illuminance/illuminance)
        phase = i / (len(data) - 1)
        ans.append((phase, magnitude_change))

    return ans

def to_pixels(angle: Union[float, Measure.Unit], pixel: Union[float, Measure.Unit] = 1/90000) -> Union[float, Measure.Unit]:
    return angle / pixel

def ring_width(vol: Union[float, Measure.Unit], sma: Union[float, Measure.Unit], eccentricity: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    return vol / (math.pi ** 2 * sma ** 2 * math.sqrt(1 - eccentricity ** 2))
    # Derivation:
        # The shape of the ring can be approximated by a torus with a circular cross-section.
        # The volume of torus is: 2π²Rr², R - torus radius, r - half the torus width
        # To account for the eccentricity of the ring, we need to stretch the torus of unit radius by a and b times:
        # The volume of stretched torus: V = 2π²abr, a - semi-major axis, b - semi-minor axis
        # b = a√(1-e²), e - eccentricity
        # V = 2π²a²r√(1-e²)
        # Ring width: w = 2r = 2(V/(2π²a²√(1-e²)))
        # w = V/(π²a²√(1-e²))


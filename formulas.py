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

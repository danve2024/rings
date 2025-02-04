import math
from typing import Union
from measure import Measure
from units import au, yrs

# Independent formulas used for model calculation

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

def volume(radius: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
    Computes the volume of a sphere.
    V = 4/3*πR³

    :param radius: R
    :return: V
    """
    return 4/3 * math.pi * radius**3

def roche_limit(radius: Union[float, Measure.Unit], asteroid_density: Union[float, Measure.Unit], ring_density: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
    Computes the Roche limit for an asteroid.
    a_Roche = R∛(2D/d)

    :param radius: R
    :param asteroid_density: D
    :param ring_density: d
    :return: a_Roche
    """
    return radius * (2 * asteroid_density/ring_density) ** (1/3)

def maximum_ring_mass(asteroid_mass: Union[float, Measure.Unit], asteroid_radius: Union[float, Measure.Unit], ring_sma: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
    Computes the maximum mass of a ring for the asteroid so that the ring-asteroid system is not binary.
    m_max = MR/2a

    :param asteroid_mass: M
    :param asteroid_radius: R
    :param ring_sma: a
    :return: m_max
    """
    return (asteroid_mass * asteroid_radius)/(2 * ring_sma)

def angular_area(angular_radius: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
    Computes the area of a circle on a sphere.
    Ω = πρ²

    :param angular_radius: ρ - angular radius
    :return: Ω - angular area
    """
    return math.pi * angular_radius**2

def angular_diameter(radius: Union[float, Measure.Unit], distance: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
    Computes the angular diameter of a celestial body.
    2ρ = 2arctan(r/l) ≈ 2r/l

    :param radius: r - radius of the body
    :param distance: l - distance
    :return: 2ρ - angular diameter
    """
    return 2 * math.degrees(radius / distance)

def to_angle(size: Union[float, Measure.Unit], distance: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
        Computes the angular size of a celestial body.
        α = 2arctan(L/2l) ≈ L/l

        :param size: L - linear size of the body
        :param distance: l - distance
        :return: α - angular size
        """
    return 2 * math.degrees(size/2 / distance)

def period(sma: Union[float, Measure.Unit]) -> Union[float, Measure]:
    """
    Computes the orbital period of a celestial body using the 2nd Kepler's Law for the Solar System.
    T[years] = √(a[au]³)

    :param sma: a - semi-major axis
    :return: T - orbital period
    """
    return math.sqrt((sma/au)**3) * yrs

def velocity(sma: Union[float, Measure.Unit], time: Union[float, Measure.Unit]) -> Union[float, Measure]:
    """
    Computes the orbital velocity of a celestial body.
    v = 2πA/t

    :param sma: A
    :param time: t - time
    :return: v - orbital velocity
    """
    return (2 * math.pi * sma)/time

def angular_velocity(sma: Union[float, Measure.Unit], time: Union[float, Measure.Unit]) -> Union[float, Measure]:
    """
    Computes the angular velocity of a celestial body.
    ω = v(A, t)/A

    :param sma: A
    :param time: t - time
    :return: ω - angular velocity
    """
    return math.degrees(velocity(sma, time) / sma)

def format_data(data: list) -> list:
    """
    Formats the data created by simulation for creating a lightcurve Δm(Φ)
    [I1, I2, ..., Ii] -> [(Φ1, Δm1), (Φ2, Δm2), ..., (Φi, Δmi)]

    Δmi = 2.5log(I0/Ii)
    Φi = i/(|data| - 1)

    I - illuminance
    Φ - occultation phase
    Δm - magnitude change

    :param list data: [I1, I2, ..., Ii]
    :return: [(Φ1, Δm1), (Φ2, Δm2), ..., (Φi, Δmi)]
    :rtype: list
    """

    ans = []
    initial_illuminance = data[0]
    for i in range(len(data)):
        illuminance = data[i]
        magnitude_change = 2.5 * math.log10(initial_illuminance/illuminance)
        phase = i / (len(data) - 1)
        ans.append((phase, magnitude_change))

    return ans

def to_pixels(angle: Union[float, Measure.Unit], pixel: Union[float, Measure.Unit] = 1/900000) -> Union[float, Measure.Unit]:
    """
    Converts angular units to pixel units for creating mask matrices.

    :param angle: angular unit
    :param pixel: pixel unit
    :return: converted angular unit to pixel unit
    """
    return angle / pixel

def ring_width(vol: Union[float, Measure.Unit], sma: Union[float, Measure.Unit], eccentricity: Union[float, Measure.Unit]) -> Union[float, Measure.Unit]:
    """
    Computes the ring width.
    w = V_ring/(π²a²√(1-e²))

    :param vol: V_ring - ring volume
    :param sma: a
    :param eccentricity: e
    :return: w
    """
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


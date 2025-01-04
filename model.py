from measure import Measure
from units import *
from formulas import *
from running import Progress
from space import Asteroid, Rings, Star, hill_sphere

star = Star(6.0 * mag, 0.8 * arcsec)

def calculate_parameters(radius: Measure=None, density: Measure=None, asteroid_sma: Measure=None, eccentricity: Measure=None, inclination: Measure=None, ring_density: Measure=None, filename='model.txt', show_progress=True, print_values=False):
    #Default parameters
    #Asteroid
    if radius is None:
        radius = Measure(1*km, 21*km, 5*km) # radius
    if density is None:
        density = Measure(1.3*gcm3, 5.3*gcm3, 1*gcm3) # density
    if asteroid_sma is None:
        asteroid_sma = Measure(35*au, 50*au, 5*au)  # semi-major axis

    # Ring
    if eccentricity is None:
        eccentricity = Measure(0, 0.8, 0.2)  # eccentricity
    if inclination is None:
        inclination = Measure(0*deg, 90*deg, 15 * deg)  # inclination
    if ring_density is None:
        ring_density = Measure(0.01 * gcm3, 0.03 * gcm3, 0.02 * gcm3)  # density

    bar = Progress('Progress',[radius, density, asteroid_sma, eccentricity, inclination, ring_density], show_progress)
    s = ''
    print('Calculating...')


    for R in radius:
        for D in density:
            for d in ring_density:
                V = volume(R)  # asteroid volume
                M = V * D  # asteroid mass
                for A in asteroid_sma:
                    a_min = max(roche_limit(R, D, d), R)  # semi-major axis minimum
                    a_max = hill_sphere(A, M)  # semi-major axis maximum
                    a_step = 10 * km # semi-major axis step
                    sma = Measure(a_min, a_max, a_step)
                    for a in sma:
                        m_max = maximum_ring_mass(M, R, a) # ring mass minimum
                        m_min = 0.5 * m_max # ring mass maximum
                        m_step = 0.1 * m_max # ring mass step
                        ring_mass = Measure(m_min, m_max, m_step)
                        for m_ring in ring_mass:
                            for e in eccentricity:
                                for i in inclination:
                                    if a == a_min and m_ring == m_min:
                                        bar.next()
                                    rings = Rings(d, a, m_ring, e, i) # create rings
                                    asteroid = Asteroid(rings, R, D, A, V, M) # create asteroid
                                    s += str(asteroid)
                                    if print_values:
                                        print(asteroid)
                                    star.occultation(asteroid)



    print('Saving...')
    with open(filename, 'w') as f:
        f.write(s)
    print('Done!')
    bar.finish()

    return filename


if __name__ == '__main__':
    calculate_parameters()


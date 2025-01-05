import numpy as np
from measure import Measure
from typing import Union
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from math import sin, radians


def disk(radius: Union[float, Measure.Unit], size: Union[float, Measure.Unit] = None) -> np.array:
    # problem with units or numpy array sizes
    if size is None:
        size = int(radius * 2)

    radius = round(radius)
    size = round(size)

    size |= 1

    center = size // 2
    ans = np.zeros((size, size))
    y, x = np.ogrid[:size, :size]
    mask = (x - center) ** 2 + (y - center) ** 2 <= radius ** 2
    ans[mask] = 1

    return ans


def elliptical_ring(size: Union[float, Measure.Unit], a: Union[float, Measure.Unit], e: Union[float, Measure.Unit], w: Union[float, Measure.Unit], i: Union[float, Measure.Unit], fill: float, focus: tuple = None) -> np.array:
    if focus is None:
        focus = (round(size/2), round(size/2))
    size = round(size)
    a = round(a)
    w = round(w)

    shape = (size, size)
    fy, fx = focus
    c = e * a
    b = np.sqrt(a ** 2 - c ** 2) * sin(radians(i))
    center_x = fx + c
    center_y = fy
    y, x = np.ogrid[:shape[0], :shape[1]]

    outer_mask = ((x - center_x) / a) ** 2 + ((y - center_y) / b) ** 2 <= 1
    inner_mask = ((x - center_x) / (a - w)) ** 2 + ((y - center_y) / (b - w)) ** 2 <= 1
    ring_mask = outer_mask & ~inner_mask

    arr = np.zeros(shape, dtype=float)
    arr[ring_mask] = fill

    return arr


def gaussian(diameter: Union[float, Measure.Unit], size: Union[float, Measure.Unit] = None) -> np.array:
    if size is None:
        size = diameter

    diameter = round(diameter)
    size = round(size)

    grid = size // 2
    x = np.linspace(-grid, grid, size)
    y = np.linspace(-grid, grid, size)
    xs, ys = np.meshgrid(x, y)
    model = Gaussian2D(amplitude=1, x_mean=0, y_mean=0, x_stddev=diameter/3, y_stddev=diameter/3)(xs, ys)

    return model


def crop(array: np.array, rows: int, end: bool = False) -> np.array:
    if array.size == 0:
        return np.array([])
    if rows > array.shape[0]:
        return np.zeros_like(array)

    result = np.zeros_like(array)
    if rows == 0:
        return array
    else:
        if end:
            result[rows:] = array[:-rows]
        else:
            result[:-rows] = array[rows:]
        return result


def cover(star: np.array, asteroid: np.array) -> list:
    data = []
    grid = (len(star), len(star[0]))

    for i in range(len(star) - 1, 0, -1):
        mask = crop(asteroid, i, True)
        result = np.zeros(grid)
        for x in range(len(star)):
            for y in range(len(star[0])):
                result[x][y] = max(0, star[x][y] - mask[x][y])
        area = np.sum(result)
        data.append(float(area))
        # show_model(result)

    for i in range(len(star)):
        mask = crop(asteroid, i)
        result = np.zeros(grid)
        for x in range(len(star)):
            for y in range(len(star[0])):
                result[x][y] = max(0, star[x][y] - mask[x][y])
        area = np.sum(result)
        data.append(float(area))
        # show_model(result)

    return data

def show_model(model: np.array):
    plt.figure()
    plt.imshow(model)
    plt.show()

if __name__ == '__main__':
    show_model(disk(10, 50) + elliptical_ring(51, 20, 0.2, 1, 5, 1))
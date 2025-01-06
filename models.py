import numpy as np

from formulas import format_data
from measure import Measure
from typing import Union
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from math import sin, cos, sqrt, radians


def disk(radius: float, size: float = None) -> np.ndarray:
    if size is None:
        size = 2 * radius + 1
    size = int(round(size))
    center = (size - 1) / 2.0

    ans = np.zeros((size, size), dtype=float)
    y, x = np.ogrid[:size, :size]

    mask = (x - center) ** 2 + (y - center) ** 2 <= radius ** 2
    ans[mask] = 1
    return ans


def elliptical_ring(
    size: Union[float, Measure.Unit], 
    a: Union[float, Measure.Unit], 
    e: Union[float, Measure.Unit], 
    w: Union[float, Measure.Unit], 
    i: Union[float, Measure.Unit],
    fill: float, 
    rotation_angle: Union[float, Measure.Unit] = 90,
    focus: tuple = None) -> np.array:
    size = int(round(size))
    shape = (size, size)
    a = float(a)
    w = float(w)

    if focus is None:
        cxy = (size - 1) / 2.0
        focus = (cxy, cxy)
    fy, fx = focus
    c = e * a
    b0 = sqrt(a ** 2 - c ** 2)
    b = b0 * sin(radians(i))
    
    y, x = np.ogrid[:size, :size]

    # Перемещаем координаты к центру вращения
    x_centered = x - fx
    y_centered = y - fy

    # Применяем поворот
    theta = radians(rotation_angle)
    print("rotation angle: ", rotation_angle)
    x_rot = x_centered * cos(theta) + y_centered * sin(theta)
    y_rot = -x_centered * sin(theta) + y_centered * cos(theta)

    # Смещаем центр эллипса от фокуса
    x_rot = x_rot - c

    if abs(b) < 1e-12:
        outer_mask = np.zeros(shape, dtype=bool)
        inner_mask = np.zeros(shape, dtype=bool)
    else:
        outer_mask = (x_rot / a) ** 2 + (y_rot / b) ** 2 <= 1
        a_inner = max(a - w, 1e-9)
        b_inner = max(b - w, 1e-9)
        inner_mask = (x_rot / a_inner) ** 2 + (y_rot / b_inner) ** 2 <= 1
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

    return format_data(data)

def normalize(array: np.array) -> np.array:
    for x in range(len(array)):
        for y in range(len(array[0])):
            if array[x][y] > 1:
                array[x][y] = 1
    return array

def show_model(model: np.ndarray) -> None:
    plt.figure()
    plt.imshow(model, origin='lower', cmap='viridis')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    # print((cover(gaussian(20), normalize(disk(3, 20) + elliptical_ring(20, 8, 0.5, 1, 90, 0.2)))))
    show_model(disk(2, 51) + elliptical_ring(51, 13, 0.5, 2, 40, 1, 0.0,))
    show_model(disk(2, 51) + elliptical_ring(51, 13, 0.5, 2, 40, 1, 45.0))
    show_model(disk(2, 51) + elliptical_ring(51, 13, 0.5, 2, 40, 1))


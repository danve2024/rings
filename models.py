import numpy as np
from measure import Measure
from typing import Union
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt


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

def gaussian(size: Union[float, Measure.Unit], radius: Union[float, Measure.Unit] = 1) -> np.array:
    grid = size // 2
    x = np.linspace(-grid, grid, size)
    y = np.linspace(-grid, grid, size)
    xs, ys = np.meshgrid(x, y)
    model = Gaussian2D(amplitude=1, x_mean=0, y_mean=0, x_stddev=2*radius, y_stddev=2*radius)(xs, ys)
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
    print(cover(gaussian(10, 2), disk(3, 10)))
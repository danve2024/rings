from PyQt6.QtGui import QImage
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling.models import Gaussian2D
from math import sin, cos, sqrt, radians

from formulas import format_data
from measure import Measure
from typing import Union


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

    x_centered = x - fx
    y_centered = y - fy

    theta = radians(rotation_angle)
    # print("rotation angle: ", rotation_angle)
    x_rot = x_centered * cos(theta) + y_centered * sin(theta)
    y_rot = -x_centered * sin(theta) + y_centered * cos(theta)

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


def gaussian(diameter: Union[float, Measure.Unit], size: Union[float, Measure.Unit] = None,
             std_dev: Union[float, Measure.Unit] = None) -> np.array:
    if size is None:
        size = diameter
    if std_dev is None:
        std_dev = diameter / 3

    diameter = round(diameter)
    size = round(size)

    grid = size // 2
    x = np.linspace(-grid, grid, size)
    y = np.linspace(-grid, grid, size)
    xs, ys = np.meshgrid(x, y)
    model = Gaussian2D(amplitude=1, x_mean=0, y_mean=0, x_stddev=std_dev, y_stddev=std_dev)(xs, ys)

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
    initial_area = np.sum(star)

    for i in range(len(star) - 1, 0, -1):
        mask = crop(asteroid, i, True)
        if np.count_nonzero(mask) != 0:
            result = np.zeros(grid)
            for x in range(len(star)):
                for y in range(len(star[0])):
                    result[x][y] = max(0, star[x][y] - mask[x][y])
            area = np.sum(result)
            data.append(float(area))
            # show_model(result)

    for i in range(len(star)):
        mask = crop(asteroid, i)
        if np.count_nonzero(mask) != 0:
            result = np.zeros(grid)
            for x in range(len(star)):
                for y in range(len(star[0])):
                    result[x][y] = max(0, star[x][y] - mask[x][y])
            area = np.sum(result)
            data.append(float(area))
            # show_model(result)

    data.insert(0, initial_area)
    data.append(initial_area)

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

def cover_animation(star: np.array, asteroid: np.array) -> list:
    """
    Generates frames for the asteroid ring occultation animation, preserving the original logic.

    Args:
        star: 2D NumPy array representing the star's brightness distribution.
        asteroid: 2D NumPy array representing the asteroid's mask.

    Returns:
        List of QImage frames for animation.
    """
    frames = []
    grid = (len(star), len(star[0]))

    # Phase 1: Decreasing mask size
    for i in range(len(star) - 1, 0, -1):
        mask = crop(asteroid, i, True)
        if np.count_nonzero(mask) != 0:
            result = np.zeros(grid)
            for x in range(len(star)):
                for y in range(len(star[0])):
                    result[x][y] = max(0, star[x][y] - mask[x][y])
            frame = _array_to_qimage(result)
            frames.append(frame)

    # Phase 2: Increasing mask size
    for i in range(len(star)):
        mask = crop(asteroid, i)
        if np.count_nonzero(mask) != 0:
            result = np.zeros(grid)
            for x in range(len(star)):
                for y in range(len(star[0])):
                    result[x][y] = max(0, star[x][y] - mask[x][y])
            frame = _array_to_qimage(result)
            frames.append(frame)

    return frames


def _array_to_qimage(array: np.array) -> QImage:
    """
    Converts a 2D NumPy array to a QImage for animation display.

    Args:
        array: 2D NumPy array.

    Returns:
        QImage representing the array as a grayscale image.
    """
    # Normalize array to 0-255
    normalized = (array - np.min(array)) / (np.max(array) - np.min(array) + 1e-8) * 255
    img_array = normalized.astype(np.uint8)
    height, width = img_array.shape

    # Create QImage with a copy to avoid segmentation fault
    return QImage(img_array.data, width, height, width, QImage.Format.Format_Grayscale8).copy()


if __name__ == '__main__':
    print(cover(gaussian(20), normalize(disk(3, 20) + elliptical_ring(20, 9, 0.3, 1, 30, 0.5))))


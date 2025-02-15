import math

from PyQt6.QtGui import QImage
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling.models import Gaussian2D
from math import sin, cos, sqrt, radians

from formulas import format_data, to_angle
from measure import Measure
from typing import Union

"""
Uses numpy arrays to model the event of the occultation

0 - minimum matrix pixel brightness
1 - maximum matrix pixel brightness
"""


def disk(radius: float, size: float = None) -> np.ndarray:
    """
    Draws a disk of specified radius

    :param radius: disk radius (in pixel units)
    :param size: matrix size (used for matrix concatenation)
    :return: numpy array with a disk of specified radius
    """
    if size is None:
        size = 2 * radius + 1
    size = int(round(size))
    if size == 0:
        return np.zeros((1, 1))
    if size % 2 == 0:
        size -= 1

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
    """
    Draws an elliptical ring of specified parameters

    :param size: matrix size (used for matrix concatenation)
    :param a: ring semi-major axis (in pixel units)
    :param e: ring eccentricity
    :param w: ring width (in pixel units)
    :param i: ring inclination (in degrees)
    :param fill: ring fill percentage (depends on the absorption coefficient)
    :param rotation_angle: rotation angle (in degrees)
    :param focus: focus coordinates of the ring (in pixel units)
    :return: numpy array with an elliptical ring of specified parameters
    """

    size = int(round(size))

    if size == 0:
        return np.zeros((1, 1))

    if size % 2 == 0:
        size -= 1

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
    """
    Will be replaced with the square root law model for applying to the star's physical parameters (T, log(g))!

    Draws a Gaussian star model of specified diameter

    :param diameter: disk diameter (in pixel units)
    :param size: matrix size (used for matrix concatenation)
    :param std_dev: standard deviation of the Gaussian distribution (in pixel units)
    :return: numpy array with a Gaussian star model of specified diameter
    """
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


def star_model(shape: list[int], coefficients: list[float]) -> np.array:
    """
    Creates a star model using the square root limb darkening approximation
    I(μ) = I(1)[1-c(1-μ)-d(1-√μ)]

    I - pixel brightness in relative units (0 - min, 1 - max)
    μ - cosine of the angle by the emergent radiation and the direct perpendicular to the stellar surface

    :param shape: shape of the star model (m x n)
    :param coefficients: coefficients of the star model (c, d)
    :param angular_size: angular radius of the star
    :return: star model as a numpy array
    """

    m, n = shape
    c, d = coefficients

    if m == 0:
        return np.zeros((1, 1))
    if n == 0:
        return np.zeros((1, 1))

    if m % 2 == 0:
        m -= 1

    if n % 2 == 0:
        n -= 1

    # Create an empty matrix
    result = np.zeros((m, n))

    # Calculate the center of the matrix
    center_x = (m - 1) / 2
    center_y = (n - 1) / 2

    # Fill the matrix with distances from the center
    for i in range(m):
        for j in range(n):
            # mu = cos((2 * np.sqrt((i - center_x) ** 2 + (j - center_y) ** 2) / m) * math.pi / 2)
            mu = cos(np.sqrt((i - center_x) ** 2 + (j - center_y) ** 2)/(n / 2) * math.pi / 2)
            # result[i, j] = mu
            if mu < 0:
                mu = 0
            result[i, j] = 1 - c * (1 - mu) - d * (1 - math.sqrt(mu))
            if result[i, j] < 0:
                result[i, j] = 0

    return np.rot90(result)


def crop(array: np.array, rows: int, shape: tuple[int] = None, end: bool = False) -> np.array:
    """
    Crops the given array to the specified number of rows. Used for modeling the occultation stages.

    :param np.array array: numpy array to be cropped
    :param int rows: number of rows to crop
    :param tuple shape: shape of the cropped array (optional, default is the original shape)
    :param bool end: if True, the cropped array will start from the end
    :return: cropped array
    """

    if not shape:
        shape = array.shape

    if array.size == 0:
        return np.array([])
    if rows > shape[0]:
        return np.zeros_like(shape)

    result = np.zeros_like(shape)
    if rows == 0:
        return array
    else:
        if array.shape[0] >= shape[0]:
            if end:
                    result[rows:] = array[:-rows]
            else:
                result[:-rows] = array[rows:]
        return result


def cover_old(star: np.array, asteroid: np.array, initial: float) -> list:
    """
    Models the occultation by masking the star array with the array of the asteroid with rings

    :param np.array star: the array of the star covered
    :param np.array asteroid: the array of the covering asteroid with its rings
    :param float initial: the initial intensity of the star radiation
    :return: the list of data used for drawing a lightcurve
    """

    m, n = star.shape

    data = []
    initial_intensity = initial # calculating initial illuminance
    cutout_intensity = np.sum(star)

    # Phase 1: Decreasing mask size
    for i in range(m - 1, 0, -1):
        print(i)
        mask = crop(asteroid, i, star.shape,True)
        show_model(mask)
        if np.count_nonzero(mask) != 0:
            result = np.zeros(star.shape)
            for x in range(m):
                for y in range(n):
                    result[x][y] = max(0, star[x][y] - mask[x][y])
            intensity = np.sum(result) # calculating the sum of intensities from each pixel - Ii
            data.append(float(initial_intensity - cutout_intensity + intensity)) # adding Ii to the output array
            # show_model(result)

    # Phase 2: Increasing mask size
    for i in range(m):
        mask = crop(asteroid, i, star.shape)
        if np.count_nonzero(mask) != 0:
            result = np.zeros(star.shape)
            for x in range(m):
                for y in range(n):
                    result[x][y] = max(0, star[x][y] - mask[x][y])
            intensity = np.sum(result)  # calculating the sum of intensities from each pixel - Ii
            data.append(float(initial_intensity - cutout_intensity + intensity))  # adding Ii to the output array
            # show_model(result)

    data.insert(0, initial_intensity)
    data.append(initial_intensity)

    return format_data(data) # [I1, I2, ..., Ii] -> [(Φ1, Δm1), (Φ2, Δm2), ..., (Φi, Δmi)]


def cover(star: np.array, mask: np.array, initial: float) -> list:
    """
        Models the occultation by masking the star array with the array of the asteroid with rings

        :param np.array star: the array of the star covered
        :param np.array mask: the array of the covering asteroid with its rings
        :param float initial: the initial intensity of the star radiation
        :return: the list of data used for drawing a lightcurve
    """

    m, n = star.shape
    n_cover = mask.shape[0]

    data = []
    initial_intensity = initial  # calculating initial illuminance
    cutout_intensity = np.sum(star)

    if mask.shape[1] != n:
        raise ValueError(f"Mask array must have {n} columns to match the base array. But it has {mask.shape[1]}.")

    for start_row in range(m - 1, -n_cover, -1):
        result = star.copy()

        overlap_start = max(start_row, 0)
        overlap_end = min(start_row + n_cover, m)

        cover_start = overlap_start - start_row
        cover_end = cover_start + (overlap_end - overlap_start)

        # Apply the mask with max(0, star[x][y] - mask[x][y])
        result[overlap_start:overlap_end, :] = np.maximum(
            0, star[overlap_start:overlap_end, :] - mask[cover_start:cover_end, :]
        )

        # show_model(result)

        intensity = np.sum(result)  # calculating the sum of intensities from each pixel - Ii
        data.append(float(initial_intensity - cutout_intensity + intensity))  # adding Ii to the output array

    data.insert(0, initial_intensity)
    data.append(initial_intensity)

    return format_data(data)  # [I1, I2, ..., Ii] -> [(Φ1, Δm1), (Φ2, Δm2), ..., (Φi, Δmi)]

def normalize(array: np.array) -> np.array:
    """
    Normalizes the given array to the range [0, 1]

    :param np.array array: the array to normalize
    :return: normalized array
    """
    for x in range(len(array)):
        for y in range(len(array[0])):
            if array[x][y] > 1:
                array[x][y] = 1
    return array

def show_model(model: np.ndarray) -> None:
    """
    Shows the given model as an image using matplotlib

    :param np.ndarray model: the model to show
    """
    plt.figure()
    plt.imshow(model, origin='lower', cmap='viridis')
    plt.colorbar()
    plt.show()

def cover_animation_old(star: np.array, asteroid: np.array) -> list:
    """
    Generates frames for the asteroid ring occultation animation, preserving the original logic.

    :param np.array star: the array of the star covered
    :param np.array asteroid: the array of the covering asteroid with its rings
    :return: list of QImage frames for animation
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

def cover_animation(star: np.array, mask: np.array) -> list:
    """
        Models the occultation by masking the star array with the array of the asteroid with rings

        :param np.array star: the array of the star covered
        :param np.array mask: the array of the covering asteroid with its rings
        :param float initial: the initial intensity of the star radiation
        :return: the list of data used for drawing a lightcurve
    """

    m, n = star.shape
    n_cover = mask.shape[0]

    frames = []

    if mask.shape[1] != n:
        raise ValueError(f"Mask array must have {n} columns to match the base array.")

    for start_row in range(m - 1, -n_cover, -1):
        result = star.copy()

        overlap_start = max(start_row, 0)
        overlap_end = min(start_row + n_cover, m)

        cover_start = overlap_start - start_row
        cover_end = cover_start + (overlap_end - overlap_start)

        # Apply the mask with max(0, star[x][y] - mask[x][y])
        result[overlap_start:overlap_end, :] = np.maximum(
            0, star[overlap_start:overlap_end, :] - mask[cover_start:cover_end, :]
        )

        frames.append(_array_to_qimage(result))  # adding a new frame

    return frames

def _array_to_qimage(array: np.array) -> QImage:
    """
    Converts a 2D numpy array to a QImage for animation display.

    :param np.array array: array to convert
    :return: QImage representing the array as a grayscale image.
    """
    # Normalize array to 0-255
    normalized = (array - np.min(array)) / (np.max(array) - np.min(array) + 1e-8) * 255
    img_array = normalized.astype(np.uint8)
    height, width = img_array.shape

    # Create QImage with a copy to avoid segmentation fault
    return QImage(img_array.data, width, height, width, QImage.Format.Format_Grayscale8).copy()


if __name__ == '__main__':
    cover(star_model([10, 100], [0.933, 0.922]), disk(5), 500)


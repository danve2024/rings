import numpy as np
from measure import Measure
from typing import Union


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


def crop(array: np.array, rows: int, end: bool = False) -> np.array:
    if array.size == 0:
        return np.array([])

    num_rows = array.shape[0]
    if rows > num_rows:
        new_array = np.zeros_like(array)
        return new_array

    new_array = np.zeros_like(array)

    if end:
        start_index = num_rows - rows
        new_array[:rows] = array[start_index:]
    else:
        start_index = num_rows - rows
        new_array[start_index:] = array[:rows]

    return new_array


def cover(star: np.array, asteroid: np.array) -> list:
    data = []
    initial_area = np.count_nonzero(star)

    for i in range(len(star)):
        mask = crop(asteroid, i)
        result = (star.astype(int) & mask.astype(int)).astype(float)
        area = initial_area - np.count_nonzero(result)
        data.append(area)

    for i in range(len(star) - 1, -1, -1):
        mask = crop(asteroid, i, True)
        result = (star.astype(int) & mask.astype(int)).astype(float)
        area = initial_area - np.count_nonzero(result)
        data.append(area)

    return data

if __name__ == '__main__':
    print(cover(disk(10), disk(1, size=20)))

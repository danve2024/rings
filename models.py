import numpy as np
from measure import Measure
from typing import Union


def disk(radius: Union[float, Measure.Unit], size: Union[float, Measure.Unit] = None) -> np.array:
    # problem with units or numpy array sizes
    if size is None:
        size = int(radius * 2)

    radius = round(radius)
    size = round(size)

    size = size if size % 2 != 0 else size + 1

    center = size // 2
    ans = np.zeros((size, size))
    for x in range(size):
        for y in range(size):
            if (x - center) ** 2 + (y - center) ** 2 <= radius ** 2:
                ans[x, y] = 1
    return ans

def cover(star: np.array, asteroid: np.array) -> list:
    data = []
    initial_area = np.count_nonzero(star)
    area = np.count_nonzero(star)

    for i in range(len(star)):
        pass

    return data

if __name__ == '__main__':
    print(cover(disk(10), disk(1)))

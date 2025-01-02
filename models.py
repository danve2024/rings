import numpy as np
from measure import Measure
from typing import Union


def disk(radius: Union[float, Measure.Unit], size: Union[float, Measure.Unit] = None) -> np.array:
    if size is None:
        size = int(radius * 2)

    size = size if size % 2 != 0 else size + 1

    center = size // 2
    ans = np.zeros((size, size))
    for x in range(size):
        for y in range(size):
            if (x - center) ** 2 + (y - center) ** 2 <= radius ** 2:
                ans[x, y] = 1
    return ans

if __name__ == '__main__':
    print(disk(10))

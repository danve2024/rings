from typing import Union

class Measure:
    def __init__(self, minimum: float, maximum: float, unit: Union[int, float] = 1, label: str = None, print_values=False):
        self.min = self.Unit(minimum).set(unit)
        self.max = self.Unit(maximum).set(unit)
        self.unit = self.Unit(unit)
        self.range = self.max - self.min
        self.label = label
        if self.label:
            self.prefix = self.label + ': '
        else:
            self.prefix = ''
        if print_values:
            print(self)

    def slider(self, value: int):
        return self.min(self.unit) + self.range(self.unit) * (value / 100)

    def update(self, minimum: float, maximum: float):
        self.__init__(minimum, maximum, self.unit, self.label)

    def __str__(self):
        return f'{self.prefix}measure({self.min}, {self.max}, {self.unit})'

    class Unit(float):
        def __add__(self, other):
            return Measure.Unit(super().__add__(float(other)))

        def __sub__(self, other):
            return Measure.Unit(super().__sub__(float(other)))

        def __mul__(self, other):
            return Measure.Unit(super().__mul__(float(other)))

        def __truediv__(self, other):
            return Measure.Unit(super().__truediv__(float(other)))

        def __pow__(self, other):
            return Measure.Unit(super().__pow__(float(other)))

        def __radd__(self, other):
            return Measure.Unit(super().__add__(float(other)))

        def __rsub__(self, other):
            return Measure.Unit(Measure.Unit(other).__sub__(float(self)))

        def __rmul__(self, other):
            return Measure.Unit(super().__mul__(float(other)))

        def __rtruediv__(self, other):
            return Measure.Unit(Measure.Unit(other).__truediv__(float(self)))

        def __rpow__(self, other):
            return Measure.Unit(Measure.Unit(other).__pow__(float(self)))

        def __floordiv__(self, other):
            return round(self/other) + 1

        def __call__(self, n):
            return self/n

        def set(self, n):
            return self*n
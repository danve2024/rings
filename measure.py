class Measure:
    def __init__(self, minimum: float, maximum: float, step: float):
        self.min = self.Unit(minimum)
        self.max = self.Unit(maximum)
        self.delta = self.Unit(step)

    def __iter__(self) -> iter:
        return iter(self.min + k * self.delta for k in range((self.max - self.min) // self.delta))

    def __len__(self) -> int:
        return (self.max - self.min)//self.delta

    def __add__(self, other) -> list:
        return list(self) + list(other)

    def __radd__(self, other) -> list:
        return list(other) + list(self)

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

def single(value: float) -> list[Measure.Unit]:
    return [Measure.Unit(value)]

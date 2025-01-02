from progress.bar import IncrementalBar

class Progress:
    @staticmethod
    def if_show(f):
        def execute(self):
            if self.show_progress:
                f(self)
                return f

        return execute

    def __init__(self, name: str, parameters: list, show_progress: bool, operation_time: float=100):
        self.name = name
        self.parameters = parameters
        self.show_progress = show_progress
        self.bar = None
        print(f'Estimated time: {(operation_time * self.len())//3600000}hrs')
        self.init()

    @if_show
    def init(self):
        self.bar = IncrementalBar(self.name, max=self.len())

    def len(self):
        ans = 1
        for item in self.parameters:
            ans *= len(item)
        return ans

    @if_show
    def next(self, n=1):
        self.bar.next(n)

    @if_show
    def finish(self):
        self.bar.finish()

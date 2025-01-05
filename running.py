import time
from progress.bar import IncrementalBar
from math import prod

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
        self.start_time = time.time()
        self.init()

    @if_show
    def init(self) -> None:
        self.bar = IncrementalBar(self.name, max=self.len())

    def len(self) -> int:
        return prod(len(item) for item in self.parameters)

    @staticmethod
    def format_time(seconds) -> str:
        hrs, rem = divmod(seconds, 3600)
        mins, secs = divmod(rem, 60)
        return f"{int(hrs):02}:{int(mins):02}:{int(secs):02}"

    @if_show
    def next(self, n=1) -> None:
        elapsed_time = time.time() - self.start_time
        completed = self.bar.index + n
        total = self.len()
        if completed > 0:
            estimated_total_time = (elapsed_time / completed) * total
        else:
            estimated_total_time = 0

        elapsed_formatted = self.format_time(elapsed_time)
        total_formatted = self.format_time(estimated_total_time)

        self.bar.suffix = f" {self.bar.index}/{self.bar.max} [{elapsed_formatted}/{total_formatted}]"
        self.bar.next(n)

    @if_show
    def finish(self) -> None:
        self.bar.finish()
        total_time = time.time() - self.start_time
        print(f'Total time elapsed: {self.format_time(total_time)}')

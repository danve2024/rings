import csv

class Observations:
    """
    A class for working with the observation data in CSV format.
    """
    def __init__(self, filename: str):
        """
        Converting the CSV data to a list of points for the graph for comparing with the model
        :param filename: CSV file
        """
        self.source_data = []

        with open(filename, 'r', encoding='utf-8-sig') as data:
            for line in csv.reader(data):
                self.source_data.append((float(line[0]), float(line[1])))

        self.data = self.source_data.copy()

    def __str__(self):
        return str(self.data)

    def normalize_and_shift(self, occultation_duration: float, phase_shift: float, magnitude_shift: float):
        """
        Normalizes and moves the source data along the axes. It is needed for working with the magnitude change (produced by the model) instead of the magnitude itself (observed) and normalizing the time span (using phase instead of time).

        :param occultation_duration: occultation duration (for normalization time -> phase)
        :param phase_shift: shift along the x-axis for selecting a specific time period required (time -> phase)
        :param magnitude_shift: shift of the magnitude (magnitude -> magnitude change)
        :return:
        """
        ending_point = phase_shift + 1
        new_data = []
        for i in self.source_data:
            if ending_point > i[0] / occultation_duration > phase_shift:
                new_data.append((i[0] / occultation_duration - phase_shift, i[1] + magnitude_shift))
        self.data = new_data

    def shift(self, phase_shift: float, magnitude_shift: float):
        """
        Moves the pre-processed data along the axes. It is needed for working with the magnitude change (produced by the model) instead of the magnitude itself (observed) and selecting a specific observation time (trimming the observations, using phase instead of time).
        :param phase_shift: shift along the x-axis for selecting a specific time period required (time -> phase)
        :param magnitude_shift: shift of the magnitude (magnitude -> magnitude change)
        :return:
        """
        ending_point = phase_shift + 1
        new_data = []
        for i in self.data:
            if ending_point > i[0] > phase_shift:
                new_data.append((i[0] - phase_shift, i[1] + magnitude_shift))
        self.data = new_data

    def time_to_phase(self, occultation_duration: int):
        """
        Normalizes the time, converting it into the phase parameter.
        :param occultation_duration: occultation duration
        :return:
        """
        new_data = []
        for i in self.source_data:
            new_data.append((i[0] / occultation_duration, i[1]))
        self.data = new_data

    def convert_hours(self):
        """
        Converts the time values in the source data from hours to seconds.
        """
        self.normalize_x(3600)

    def convert_minutes(self):
        """
        Converts the time values in the source data from minutes to seconds.
        """
        self.normalize_x(60)

    def convert_seconds(self):
        """
        Converts the time values in the source data from seconds to seconds.
        """
        pass

    def convert_milliseconds(self):
        """
        Converts the time values in the source data from milliseconds to seconds.
        """
        self.normalize_x(0.001)

    def convert_microseconds(self):
        """
        Converts the time values in the source data from microseconds to seconds.
        """
        self.normalize_x(10^(-6))

    def convert_kiloseconds(self):
        """
        Converts the time values in the source data from kiloseconds to seconds.
        """
        self.normalize_x(1000)

    def convert_megaseconds(self):
        """
        Converts the time values in the source data from megaseconds to seconds.
        """
        self.normalize_x(1000000)

    def normalize_x(self, time_span: float):
        """
        Normalizes the time, converting it into the phase parameter in the source data array.
        :param time_span:
        :return:
        """
        new_data = []
        for i in self.source_data:
            new_data.append((i[0] / time_span, i[1]))
        self.source_data = new_data

    def crop_x(self, starting_point: float = 0., duration: float = 1.):
        """
        Crops the source data to represent the selected time period.
        :param starting_point: beginning of the trimming time
        :param duration: occultation duration
        :return:
        """
        ending_point = starting_point + duration
        new_data = []
        for i in self.data:
            if ending_point > i[0] > starting_point:
                new_data.append((i[0] - duration, i[1]))
        self.source_data = new_data

    def move_y(self, delta: float):
        """
        Shifts the magnitude by a specific value in the source data.
        :param delta: magnitude shift
        :return:
        """
        new_data = []
        for i in self.data:
            new_data.append((i[0], i[1] + delta))
        self.source_data = new_data


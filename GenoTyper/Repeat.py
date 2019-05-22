'''docstring'''

class Repeat():
    """docstring for Genotype"""
    def __init__(self, last_unit_index, number_of_units=1, repeat_unit="CTG"):
        self.last_unit_index = last_unit_index
        self.number_of_units = number_of_units

    def add_unit(self):
        self.number_of_units +=1

    def change_last_unit_index(self, index):
        self.last_unit_index = index
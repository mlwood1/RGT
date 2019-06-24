'''docstring'''

class MatchingSequence():
    """docstring for Genotype"""
    def __init__(self, sequence_string, repeat_units_count, abundance, order):
        self.sequence_string = sequence_string
        self.repeat_units_count = repeat_units_count
        self.abundance = abundance
        self.order_in_genotable = order+1 #idx starts from 1 not zero

'''docstring'''

class Repeat():
    """docstring for Genotype"""
    def __init__(self, last_unit_index,window, number_of_units=0, repeat_unit="CTG"):
        self.repeat_unit = repeat_unit
        self.last_unit_index = last_unit_index
        self.number_of_units = number_of_units
        self.start_index = last_unit_index - ((number_of_units+1)*len(repeat_unit)) #for future use
        self.single_point_mutation_indexes = []
        self.unconfirmed_units_buffer = 0
        self.non_perfect_units = 0
        self.repeat_sequence = ""
        self.unconfirmed_sequence = "" #sequence of non pure repeates waiting for a confirmed unit to  be added
        self.add_unit(window, last_unit_index)

    def add_unit(self, window, index): #index is the last base index in the sequence
        if window == self.repeat_unit: #perfect match            
            self.number_of_units +=1
            self.number_of_units += self.unconfirmed_units_buffer #add the unconfirmed units (points with SNPs)
            self.non_perfect_units += self.unconfirmed_units_buffer
            self.unconfirmed_units_buffer = 0 #zero the points with SNPs
            self.change_last_unit_index(index)
            self.repeat_sequence += self.unconfirmed_sequence
            self.unconfirmed_sequence = ""
            self.repeat_sequence += window

        else:
            self.unconfirmed_units_buffer += 1
            self.single_point_mutation_indexes.append(self.get_SNP_index(window, index))
            self.unconfirmed_sequence += window
    
    def get_SNP_index(self, window, index):
        for i in range(0,len(window)):
            if (window[i] != self.repeat_unit[i]) :
                return index -len(window)+i

    def change_last_unit_index(self, index):
        self.last_unit_index = index

    def get_non_perfect_units_percentage(self):
    	return self.non_perfect_units/self.number_of_units

    def get_seq(self, read):
    	return read[self.start_index:self.last_unit_index]

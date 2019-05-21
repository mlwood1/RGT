'''docstring'''

class Genotype():
    """docstring for Genotype"""
    def __init__(self, reads, repeat_unit="CTG", min_size_repeate=3):
        self.repeat_unit = repeat_unit
        self.reads = reads
        self.min_size_repeate = min_size_repeate

    def get_repeates(self):
        geno_table = {}
        for line in self.reads:
            self.get_line_repeates(line, geno_table)
        return geno_table

    def get_line_repeates(self, sequence, geno_table):
        window_length = len(self.repeat_unit)
        window_inside_repeates = False
        number_of_current_repeats = 0


        i = window_length
        while i < len(sequence): #sliding window
            window = sequence[i-window_length:i]
            if self.is_window_equals_repeat_unit(window, self.repeat_unit) and (not window_inside_repeates):
                window_inside_repeates = True
                number_of_current_repeats = 1
                i = i+2
            elif self.is_window_equals_repeat_unit(window, self.repeat_unit) and window_inside_repeates:
                number_of_current_repeats += 1
                i = i+2
            elif (not(self.is_window_equals_repeat_unit(window, self.repeat_unit)) and window_inside_repeates):# or (window_inside_repeates and i == len(sequence)):
                '''  first condition: when the repeat ends, second is when the read is finished'''
                window_inside_repeates = False
                if number_of_current_repeats >= self.min_size_repeate: #check that number of repeates is larger than the minimum size repeate
                    self.add_repeat_to_genotable(number_of_current_repeats, geno_table)
            i +=1   
        return geno_table
   
  
    def is_window_equals_repeat_unit(self, window, repeat_unit):
        mismatch = False
        for i in range(0,len(window)):
            if (window[i] != repeat_unit[i]) and (not mismatch):
                mismatch = True
            elif (window[i] != repeat_unit[i]) and mismatch:
                return False
        return True

    def add_repeat_to_genotable(self, number_of_repeat_units, geno_table):

        if(number_of_repeat_units in geno_table):
            geno_table[number_of_repeat_units] += 1
        else:
            geno_table[number_of_repeat_units] = 1

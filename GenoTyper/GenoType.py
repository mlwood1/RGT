'''docstring'''
from .Repeat import Repeat

class Genotype():
    """docstring for Genotype"""
    def __init__(self, reads, repeat_unit="CTG", min_size_repeate=1, max_interrupt_tract=6):
        self.repeat_unit = repeat_unit
        self.reads = reads
        self.min_size_repeate = min_size_repeate
        self.max_interrupt_tract = max_interrupt_tract

    def get_repeates(self):
        geno_table = {}
        for line in self.reads:
            self.get_line_repeates(line, geno_table)
        return geno_table

    def get_line_repeates(self, sequence, geno_table):
        window_length = len(self.repeat_unit)
        window_inside_repeates = False

        i = window_length
        repeat = None
        while i < len(sequence): #sliding window
            window = sequence[i-window_length:i]
            if self.window_enters_repeat_sequence(window, self.repeat_unit, repeat, window_inside_repeates):
                #if window detects a repeat unit, while it is not inside a repeat sequence
                repeat = Repeat(i, window) #creat a repeat object
                window_inside_repeates = True
                i = i+3 #Jumb one window
                continue
            #elif self.detect_repeat_unit_inside_repeat():
            elif self.is_window_equals_repeat_unit(window, self.repeat_unit, repeat) and window_inside_repeates:
                #if it detects a repeat while inside the repeat sequence
                repeat.add_unit(window,i) #add a repeat unit count
                i = i+3 #jumb one window
                continue

            elif (not(self.is_window_equals_repeat_unit(window, self.repeat_unit, repeat)) and window_inside_repeates) or (window_inside_repeates and i == len(sequence)):
                '''  first condition: when the repeat ends, second is when the read is finished'''
                '''  non repeat elements found while in the repeat sequence '''
                if i-repeat.last_unit_index <= self.max_interrupt_tract:
                    #ignore if length is smaller than max interrupt tract
                    i += 1
                    continue
                window_inside_repeates = False #if length is larger than max interrupt tract
                if repeat.number_of_units >= self.min_size_repeate: #check that number of repeates is larger than the minimum size repeate
                    self.add_repeat_to_genotable(repeat, geno_table)
                repeat = None

            i +=1   
        return geno_table
   
    def window_enters_repeat_sequence(self, window, repeat_unit, repeat_object, window_inside_repeates_flag):
        if not window_inside_repeates_flag:
            if self.is_window_equals_repeat_unit(window, repeat_unit, repeat_object):
                return True
        return False

    def is_window_equals_repeat_unit(self, window, repeat_unit, repeat_object):
        mismatch = False
        for i in range(0,len(window)):
            if (window[i] != repeat_unit[i]) and (not mismatch) and (repeat_object != None):
                mismatch = True
            elif (window[i] != repeat_unit[i]) and (mismatch or repeat_object == None):
                return False
        return True

    def add_repeat_to_genotable(self, repeat, geno_table):
        number_of_repeat_units = repeat.number_of_units
       
        if(number_of_repeat_units in geno_table):
            geno_table[number_of_repeat_units] += 1
        else:
            geno_table[number_of_repeat_units] = 1

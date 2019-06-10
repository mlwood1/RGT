'''docstring'''
from .PeakIdentifier import PeakIdentifier

class AllelesDetector():
    """docstring for main"""
    def __init__(self, counts_table, geno_table):

        #sorted_geno_list: sorted by abundance, 2D list => list of lists(key index 0, value 1)
        self.sorted_geno_list = (sorted(geno_table.items(), key=lambda x: x[1][0], reverse=True))

        self.possible_alleles_list = self.get_possible_alleles_list_from_sorted_geno_list()  
        peak_identifier = PeakIdentifier(counts_table)
        self.peak_repeat_counts = peak_identifier.get_peaks()

    def get_alleles(self):
        matching_sequences = self.get_matches_between_peaks_and_possible_alleles_list()

        if len(matching_sequences) == 2:
            return(matching_sequences, "ok hetero")
        elif len(matching_sequences) >2:
            return(matching_sequences[:1], "check, >2")
        elif len(matching_sequences) == 1:
            return([matching_sequences[1], matching_sequences[0]], "ok homo" )
        elif len(matching_sequences) ==0:
            return([self.sorted_geno_list[0][1][0],self.sorted_geno_list[1][1][0]])

    def get_matches_between_peaks_and_possible_alleles_list(self):
        matching_sequences = []
        for possibe_allele in self.possible_alleles_list:
            if (possibe_allele[1][1]) in self.peak_repeat_counts:
                matching_sequences.append(possibe_allele[0])
        return matching_sequences
    
    def get_possible_alleles_list_from_sorted_geno_list(self):
        most_abundant_allele_values = self.sorted_geno_list[0][1]
        
        # minimum threshold for identifing a repeat sequence as possible allele
        min_threshold_abundance = most_abundant_allele_values[0] * 0.2

        for i in range(0,len(self.sorted_geno_list)):
            #first index the repeat, second the values, third the repeat count            
            if self.sorted_geno_list[i][1][0] < min_threshold_abundance: 
                possible_allels_list = self.sorted_geno_list[0:i]
                break
        return possible_allels_list

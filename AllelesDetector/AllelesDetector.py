'''docstring'''
from .PeakIdentifier import PeakIdentifier
from .MatchingSequence import MatchingSequence

class AllelesDetector():
    """docstring for main"""
    def __init__(self, counts_table, geno_table):

        #sorted_geno_list: sorted by abundance, 2D list => list of lists(key index 0, value 1)
        self.sorted_geno_list = (sorted(geno_table.items(), key=lambda x: x[1][0], reverse=True))

        self.counts_table = counts_table
        self.possible_alleles_list = self.get_possible_alleles_list_from_sorted_geno_list()  
        peak_identifier = PeakIdentifier(counts_table)
        self.peak_repeat_counts = peak_identifier.get_peaks()
        self.plateau = 15 #plateau to ignore small peaks in the region of the detected allele


    def get_alleles(self):
        matching_sequences = self.get_matches_between_peaks_and_possible_alleles_list()

        if len(matching_sequences) == 2:
            message = "ok hetero"
            if self.peaks_list_has_peaks_bigger_than_genotyped_alleles(matching_sequences):
                message += " ,expanded may exist, please check "

            if matching_sequences[0].repeat_units_count == matching_sequences[1].repeat_units_count:
                message += " ,two matches with one repeat count, please check manually"

            return([matching_sequences[0].sequence_string,matching_sequences[1].sequence_string, message , str(self.peak_repeat_counts)])
        
        elif len(matching_sequences) >2:
            return([matching_sequences[0].sequence_string, matching_sequences[1].sequence_string, "check, more than two potential alleles >2",
                    str(self.peak_repeat_counts)])
        
        elif len(matching_sequences) == 1:
            new_matching_sequences = self.explore_if_a_close_allele_exists()
            if len(new_matching_sequences)==2 :
                return([new_matching_sequences[0].sequence_string,new_matching_sequences[1].sequence_string,
                       "ok hetero, alleles next to each other", str(self.peak_repeat_counts)])
            
            if self.peaks_list_has_peaks_bigger_than_genotyped_alleles(matching_sequences):
                other_allele = self.get_seq_from_geno_table_by_repeats_count(matching_sequences[0])
                return([matching_sequences[0].sequence_string, other_allele.sequence_string,
                       "expanded allele detected with small count, please check", str(self.peak_repeat_counts)])

            message = "ok homo"
            if matching_sequences[0].repeat_units_count != self.possible_alleles_list[0][1][1]:
                message += " most abundant repeat sequence not selected, please chek manually"
            return([matching_sequences[0].sequence_string, matching_sequences[0].sequence_string, message, str(self.peak_repeat_counts)] )
        
        

        elif len(matching_sequences) == 0:
            new_matching_sequences = self.explore_if_a_close_allele_exists()
            if len(new_matching_sequences)==2 :
                return([new_matching_sequences[0].sequence_string ,new_matching_sequences[1].sequence_string, "ok hetero, better check", str(self.peak_repeat_counts)])
            
            elif len(new_matching_sequences)==1:
                return([new_matching_sequences[0].sequence_string , new_matching_sequences[0].sequence_string, "homo, please check" , str(self.peak_repeat_counts) ] )
           
            return([self.sorted_geno_list[0][1][0],self.sorted_geno_list[1][1][0], "check no matching peaks", str(self.peak_repeat_counts) ])

    def peaks_list_has_peaks_bigger_than_genotyped_alleles(self, matching_sequences):

        largest_detected_peak = matching_sequences[0].repeat_units_count+self.plateau
       
        for matching_sequence in matching_sequences:
            if matching_sequence.repeat_units_count > largest_detected_peak:
                largest_detected_peak = matching_sequence.repeat_units_count+self.plateau
        for peak in self.peak_repeat_counts:
            if peak > largest_detected_peak:
                return True
        
        return False

    def get_seq_from_geno_table_by_repeats_count(self, matching_sequence):
        detected_peak = matching_sequence.repeat_units_count
        repeat_counts_bigger_than_detected_allele = [x for x in self.peak_repeat_counts if x >= detected_peak+self.plateau]
        
        for possibe_allele in self.sorted_geno_list:
            if (possibe_allele[1][1]) in repeat_counts_bigger_than_detected_allele:
                matching_sequence = MatchingSequence(possibe_allele[0], possibe_allele[1][1])
                return matching_sequence
        return MatchingSequence("can't find the other allele", 0) #this line should be impossible to reach


    def explore_if_a_close_allele_exists(self):
        new_matching_sequences = []
        new_peak_repeat_counts = list(self.peak_repeat_counts)

        first_key = list(self.counts_table)[0]
        counts_min_threshold = self.counts_table[first_key]* 0.5

        #adding near by points to peak list
        for number_of_repeat_counts in self.counts_table.keys():
            count = self.counts_table[number_of_repeat_counts]
           
            if number_of_repeat_counts in new_peak_repeat_counts:
                continue
            if count > counts_min_threshold:
                new_peak_repeat_counts.append(number_of_repeat_counts)
      
        #checking if a match happens with the near by points, now identified as new peaks
        for possibe_allele in self.possible_alleles_list:
            if (possibe_allele[1][1]) in new_peak_repeat_counts:
                matching_sequence = MatchingSequence(possibe_allele[0], possibe_allele[1][1])
                new_matching_sequences.append(matching_sequence)
        return new_matching_sequences
        

    def get_matches_between_peaks_and_possible_alleles_list(self):
        matching_sequences = []
        for possibe_allele in self.possible_alleles_list:
            if (possibe_allele[1][1]) in self.peak_repeat_counts:
                matching_sequence = MatchingSequence(possibe_allele[0], possibe_allele[1][1])
                matching_sequences.append(matching_sequence)
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

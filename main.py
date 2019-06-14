from FileReader.ReadFile import ReadFile
from GenoTyper.GenoType import Genotype
from ExcelExporter.ExcelExport import ExcelWriter
from AllelesDetector.AllelesDetector import AllelesDetector
from CountsPlotter.CountsPlotter import CountsPlotter

import glob
from joblib import Parallel, delayed, cpu_count


def RGT(sample):

    output_table = {} #dictionary to export result

    #sample code chopping from file + directory name
    sample_code = sample.split("/")[1]
    sample_code = sample_code.split(".")[0]
    print(sample_code)


    #read file and extract sequence from between flanks
    file = ReadFile(sample ,start_flank="CAACATGGGC", end_flank="ACCAGCCCAGCAGAACCAGTACGTCCACATTT")
    reads = file.reads #extracted reads from between flanks

    #genotype the reads (create counts table and repeat sequence abundance table)
    genotype = Genotype(reads,repeat_units=["CAG", "CAT"], unique_repeat_units=["CAG"])
    geno_table = genotype.get_geno_table() #the repeat sequence abundance table
    counts_table = genotype.get_counts_table() 
    unique_counts_table = genotype.get_unique_counts_table()

    #sort the three tables by abundance
    sorted_geno_table = dict(sorted(geno_table.items(), key=lambda x: x[1], reverse=True))
    sorted_counts_table = dict(sorted(counts_table.items(), key=lambda x: x[1], reverse=True))
    sorted_unique_counts_table = dict(sorted(unique_counts_table.items(), key=lambda x: x[1], reverse=True))

    #write three tabels to excel 
    excel_writer = ExcelWriter()
    excel_writer.add_table_to_sheet(sorted_geno_table,"genotype")
    excel_writer.add_table_to_sheet(sorted_counts_table,"counts")
    excel_writer.add_table_to_sheet(sorted_unique_counts_table,"unique counts")
    excel_writer.save_file("FilesSpecificResults/"+sample_code+".xlsx")

    #export plot
    table = genotype.get_counts_table()
    CountsPlotter.plot_counts_table(table, "FilesSpecificResults/Plots/"+sample_code+".png")

    #Automaticly detect allels from counts table and geom table
    a = AllelesDetector(sorted_counts_table,sorted_geno_table)
    output_table[sample_code] = a.get_alleles()


    return(output_table)
def get_collective_dictionary_from_list_of_output_dictionaries(list_of_output_dictionaries):
    ''' Convert the output of the parallel processing to a single dictionary''' 
    output_dictionary = {}
    for dictionary in list_of_output_dictionaries:
        key = list(dictionary.keys())[0] #get the first key(the only one)
        output_dictionary[key] = dictionary[key]
    return(output_dictionary)
    
def main():
    samples = glob.glob("SCAData/*.fastq")
    list_of_output_dictionaries = Parallel(n_jobs=cpu_count(),verbose=1)(map(delayed(RGT),(samples)))
    output_dictionary = get_collective_dictionary_from_list_of_output_dictionaries(list_of_output_dictionaries)

    collective_excel_writer = ExcelWriter()
    collective_excel_writer.add_table_to_sheet(output_dictionary,"results")
    collective_excel_writer.save_file("results.xlsx")



if __name__== "__main__":
    main()



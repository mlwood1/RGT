from rgt import RGT
import glob
import sys

from joblib import Parallel, delayed, cpu_count
from interface.interface import get_user_inputs
from interface.json_parser import extract_parameters
from ExcelExporter.ExcelExport import ExcelWriter


def get_collective_dictionary_from_list_of_output_dictionaries(list_of_output_dictionaries):
    ''' Convert the output of the parallel processing to a single dictionary''' 
    output_dictionary = {}
    for dictionary in list_of_output_dictionaries:
        key = list(dictionary.keys())[0] #get the first key(the only one)
        output_dictionary[key] = dictionary[key]
    return(output_dictionary)

def main():
    input_directory, output_directory, settings_file, number_of_threads = get_user_inputs(sys.argv[1:])
    settings = extract_parameters(settings_file)
    samples = glob.glob(input_directory + "/*.fastq")
    rgt_ = RGT(settings, input_directory, output_directory)

    if number_of_threads == None:
        number_of_threads = cpu_count()
    list_of_output_dictionaries = Parallel(n_jobs=number_of_threads,verbose=1)(map(delayed(rgt_.rgt),(samples)))
    output_dictionary = get_collective_dictionary_from_list_of_output_dictionaries(list_of_output_dictionaries)

    collective_excel_writer = ExcelWriter()
    results_headers = ["sample ID", "First allele structure", "Second allele structure",
                        "Comments and Flags", "Identified peaks"]
    collective_excel_writer.add_table_to_sheet(output_dictionary,"results", results_headers)
    collective_excel_writer.save_file(output_directory + "/ResultsSummary.xlsx")


if __name__== "__main__":
    main()
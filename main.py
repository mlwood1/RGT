import rgt
import glob
from joblib import Parallel, delayed, cpu_count


def get_collective_dictionary_from_list_of_output_dictionaries(list_of_output_dictionaries):
    ''' Convert the output of the parallel processing to a single dictionary''' 
    output_dictionary = {}
    for dictionary in list_of_output_dictionaries:
        key = list(dictionary.keys())[0] #get the first key(the only one)
        output_dictionary[key] = dictionary[key]
    return(output_dictionary)

def main():
    samples = glob.glob("SCAData/*.fastq")
    list_of_output_dictionaries = Parallel(n_jobs=cpu_count(),verbose=1)(map(delayed(rgt.RGT),(samples)))
    output_dictionary = get_collective_dictionary_from_list_of_output_dictionaries(list_of_output_dictionaries)

    collective_excel_writer = ExcelWriter()
    collective_excel_writer.add_table_to_sheet(output_dictionary,"results")
    collective_excel_writer.save_file("results.xlsx")



if __name__== "__main__":
    main()
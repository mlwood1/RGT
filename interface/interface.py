from . import check_files 
import getopt
import sys, os

def get_user_inputs(argv):
    input_directory = ''
    output_directory = ''
    settings_file = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:", ["help", "input", "output=", "sys"])
    except getopt.GetoptError:
        print ('test.py -i <input directory> -o <output directory> -s <settings json file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h','--help'):
            print('test.py -i <input directory> -o <output directory> -s <settings json file> ')
            sys.exit()
        elif opt in ("-i"):
            input_directory = arg
        elif opt in ("-o"):
            output_directory = arg
        elif opt in ("-s"):
            settings_file = arg
    

    check_files.check_input_directory(input_directory)
    check_files.check_or_create_output_directory(output_directory)
    check_files.check_settings_file(settings_file)
    

    return(input_directory, output_directory, settings_file)

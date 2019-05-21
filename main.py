from FileReader.ReadFile import ReadFile
from GenoTyper.GenoType import Genotype


def main():
    file = ReadFile("MadeUpTest.fastq",start_flank="startflank", end_flank="endflank")
    reads = file.reads
    genotype = Genotype(reads)
    print(genotype.get_repeates())

if __name__== "__main__":
    main()
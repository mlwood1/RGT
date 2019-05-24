from FileReader.ReadFile import ReadFile
from GenoTyper.GenoType import Genotype
import matplotlib.pyplot as plt
from ExcelExporter.ExcelExport import ExcelWriter


def main():
    file = ReadFile("DMGV229C.fastq",start_flank="AACGGGGCTCGAAGGGTCCT", end_flank="CAGGCCTGCAGTTTGCCCATC")
    #file = ReadFile("MadeUpTest.fastq",start_flank="startflank", end_flank="endflank")
    print("read")
    reads = file.reads
    genotype = Genotype(reads)
    table = genotype.get_repeates()
    print("table")
    sorted_table = dict(sorted(table.items()))
    

    ExcelWriter.write_to_excel(sorted_table)

    sorted_table = sorted(table.items())
    x, y = zip(*sorted_table)
    print(x,y)
    graph = plt.plot(x,y)
    plt.show()


if __name__== "__main__":
    main()
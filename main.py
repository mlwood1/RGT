from FileReader.ReadFile import ReadFile
from GenoTyper.GenoType import Genotype
import matplotlib.pyplot as plt
from ExcelExporter.ExcelExport import ExcelWriter


def main():
    #file = ReadFile("DMGV97C_S70_hexamers.fastq",start_flank="AACGGGGCTCGAAGGGTCCT", end_flank="CAGGCCTGCAGTTTGCCCATC")
    file = ReadFile("MS225-N701-A-S506-A_S37_L001_R1_001.fastq",start_flank="CCACAGCCTA", end_flank="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC") #sca
    #file = ReadFile("MadeUpTest.fastq",start_flank="startflank", end_flank="endflank")
    print("read")
    reads = file.reads
    genotype = Genotype(reads)
    table = genotype.get_geno_table()
    print("table")
    sorted_table = dict(sorted(table.items(), key=lambda x: x[1], reverse=True))
    #sorted_table = dict(sorted(table.items()))
    print("sorted")

    ExcelWriter.write_to_excel(sorted_table)

    #sorted_table = sorted(table.items(),reverse=True)
    #x, y = zip(*sorted_table)
    #print(x,y)
    #graph = plt.plot(x,y)
    #plt.show()


if __name__== "__main__":
    main()

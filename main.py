from FileReader.ReadFile import ReadFile
from GenoTyper.GenoType import Genotype
import matplotlib.pyplot as plt
from ExcelExporter.ExcelExport import ExcelWriter


def main():
    file_name = "MS225-N701-A-S503-A_S13_L001_R1_001"
    file = ReadFile("DMGV12C_S58.fastq",start_flank="AACGGGGCTCGAAGGGTCCT", end_flank="CAGGCCTGCAGTTTGCCCATC")
    #file = ReadFile(file_name+".fastq" ,start_flank="GTGTATGGGC", end_flank="ACCAGCCCAGCAGAACCAGTACGTCCACATTT") #sca
    #file = ReadFile(file_name+".fastq") #sca

    #file = ReadFile("MadeUpTest.fastq",start_flank="startflank", end_flank="endflank")
    print("read")
    reads = file.reads
    genotype = Genotype(reads)
    geno_table = genotype.get_geno_table()
    counts_table = genotype.get_counts_table()
    unique_counts_table = genotype.get_unique_counts_table()

    print("table")
    sorted_geno_table = dict(sorted(geno_table.items(), key=lambda x: x[1], reverse=True))
    sorted_counts_table = dict(sorted(counts_table.items(), key=lambda x: x[1], reverse=True))
    sorted_unique_counts_table = dict(sorted(unique_counts_table.items(), key=lambda x: x[1], reverse=True))
    print("sorted")
    excel_writer = ExcelWriter()
    excel_writer.add_table_to_sheet(sorted_geno_table,"genotype")
    excel_writer.add_table_to_sheet(sorted_counts_table,"counts")
    excel_writer.add_table_to_sheet(sorted_unique_counts_table,"unique counts")
    excel_writer.save_file()
   # excel_writer.save_file()
    #ExcelWriter.write_to_excel(sorted_table)

    table = genotype.get_counts_table()
    sorted_table = sorted(table.items(),reverse=True)
    x, y = zip(*sorted_table)
    #print(x,y)
    graph = plt.plot(x,y)
    plt.show()


if __name__== "__main__":
    main()

from FileReader.ReadFile import ReadFile
from GenoTyper.GenoType import Genotype
from ExcelExporter.ExcelExport import ExcelWriter
from AllelesDetector.AllelesDetector import AllelesDetector
import matplotlib.pyplot as plt
import glob


def main():

    output_table = {}
    for sample in glob.glob("HuntingtonData/*.fastq"):
        sample_code = sample.split("/")[1]
        sample_code = sample_code.split(".")[0]
        print(sample_code)

        file = ReadFile(sample ,start_flank="GTCCCTCAAGTCCTTC", end_flank="CAGCTTCCTCAGCCGC") #huntington
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
        excel_writer.save_file("FilesSpecificResults/"+sample_code+".xls")

        '''
        table = genotype.get_counts_table()
        sorted_table = sorted(table.items(),reverse=True)
        x, y = zip(*sorted_table)
        graph = plt.plot(x,y)
        plt.show()
        '''
        a = AllelesDetector(sorted_counts_table,sorted_geno_table)
        output_table[sample_code] = a.get_alleles()

        #print(output_table)
        collective_excel_writer = ExcelWriter()
        collective_excel_writer.add_table_to_sheet(output_table,"results")
        collective_excel_writer.save_file("results.xls")




if __name__== "__main__":
    main()

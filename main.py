from FileReader.ReadFile import ReadFile

def main():
    file = ReadFile("AtypicalAllele_R1.fastq",start_flank="NTGCGACCCTGGAAAAGC", end_flank="CCCACCACCACA")
    print(len(file.reads))
    print(len(file.raw_reads))

if __name__== "__main__":
    main()
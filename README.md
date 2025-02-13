[![DOI](https://zenodo.org/badge/187811748.svg)](https://doi.org/10.5281/zenodo.14865580)

# RGT
Repeats Genotyping Tool

## This is a developers guide, if you are looking for a user manual it's [here](https://github.com/hossam26644/RGT/wiki)

RGT is a software to:
1. Extract structure from fastq files representing targeted sequencing samples of SSRs (simple sequence repeats)
2. Export the identified germline/inherited alleles of the sample (assuming it is a human sample)

## System requirement
To run the RGT executable on Mac you need: macOS 12.7.5 or above plus the macOS RGT executable available [here](https://sourceforge.net/projects/rgt/files/macOS/)

To run the RGT executable on linux you need: Ubuntu 18.10 or above plus the linux RGT executable available [here](https://sourceforge.net/projects/rgt/files/Linux/). We have only tested RGT on Ubuntu but it should also run fine on other linux OS.

If you are planning to run it from the source code, the software dependencies are listed in [requirements.txt](https://github.com/hossam26644/RGT/blob/master/requirements.txt)

no non-standard harware required to run RGT

## The approach taken here is an alignment-free approach
* Because targeted sequencing gives thousands of reads, why do alignment
* We just extract structures and get the most repeating ones
* This of course means you are loosing a percentage of the reads as any sequencing error will result in a mistakenly identified repeat structure, but who cares if you have tens of thousands of reads

## Structures are identified by the expected repeating units defined by the user
* structures are stored with their number of repeating units
* A graph of the number of repeating units count and the abundance of allele structures having the same count is plotted
* Peaks from this graph are used to identify the modal alleles which are assumed to correspond to the germline/inherited alleles

> Every allele structure identified is stored regardless of the sequencing errors and variations and exported in a separate excel file for the user
![Semantic description of image](https://github.com/hossam26644/RGT/blob/enh/pep8/ImagesForReadme/CountsPlot.png "Example of a a counts plot")*Example of a counts plot*
##### The most abundant repeat structures are matched with the peaks identified in the "repeating units count vs structures having this count plot" to get the two human germline alleles

## lets get to the developers guide:
#### first the strucutre of the software:
```
RGT
│─── main.py               
│─── RGT.py                
│─── settings.json    
│─── README.md
│ 
└──────────── interface
│             │───  interface.py
│             │───  check_files.py
│             │───  json_parser.py
│             └───  __init__.py
│ 
└──────────── filereader
│             │─── ReadFile.py 
│             └───  __init__.py
│
└──────────── excelexporter
│             │───  ExcelExport.py 
│             └───  __init__.py
│ 
└──────────── genotyper
│             │───  GenoType.py
│             │───  Repeat.py
│             │───  SmartString.py
│             │───  GroupingString.py
│             │───  revComplementry.py
│             └───  __init__.py
│ 
└──────────── graphsplotter
│             │───  plotter.py
│             │───  table_2d_plotter.py   
│             │───  plot_3D.py
│             └───  __init__.py
│ 
└──────────── allelesdetector
              │───  AllelesDetector.py
              │───  MatchingSequence.py
              │───  PeakIdentifier.py
              └───  __init__.py

```
## how parallel processing is done:
* The main class uses joblib parallel, it does the job for you
* Processes are created by it and instances of the RGT class where the analysis is done are excuted in parallel
* The results of all processes are put in a list of dictionaries, there is a method *get_collective_dictionary_from_list_of_output_dictionaries* (do you like my naming style?) to (as you would imagine) get a collective dictionary from the list of output dictionaries

## Do I have to explain every single detail? sorry I am not going to
* the **interface** gets user inputs, nothing smart in that 
* **excelexporter** exports tables to excel sheets, nothing smart also I guess
* **filereader** parses fastq files and extracts the reads:
  * it searches for the flanking sequences and extract reads from between these flanks
  * flanking sequence match is identified when a sequence has a [hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) of one or zero with the flank
  * when there is no upstream flank (i.e. start_flank in settings file) found (if user specified a start_flank) the whole read is discarded
  * if the user specifies to discard reads when there is no downstream flank (i.e. end_flank in settings file).. and no downstream flank is found: guess what, the read is discarded
  * The percentage of the discarded reads is calculated and the sample is flagged if this percentage is greater than a threshold identified by the user (or default value set by the software)
  
## The genotyper:
It is the package that is responisble for extracting the repeat structure from reads, and to count the repeat units 

  #### first lets discuss other modules in the genotyper pacakage:
  * **revComplementary**: computes the reverse complement of a sequence, it's only useful if the user is working on the reverse strand
  * **GroupingString**: groups a sequence by the user input grouping units 
  e.g *CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTCCT* is grouped to *[CAG]15[CAACAG]1[CCGCCA]1[CCG]10[CCT]2* when the user input grouping units are: CAG CAACAG CCGCCA CCG CCT
    * grouping is done by using a sliding windows that identifies the user input units then slides to count them, eventualy it substitutes the repeating units with one unit between square brackets and the number of copies of this unit 
  * **SmartString** does the grouping when the user identifies no grouping units, 
    >it uses a slow runner fast runner algorithm, where a slow runner sliding window identifies repeat units, and the fast runner window counts the repeating units. This code is executed recursively on ungrouped parts of the sequence with increasing window size to detect larger repeating units. The stopping condition is to have a window size larger than half the sequence length -because there will be no possibility to have a longer repeating unit-.
    
  * **Repeat**: the repeat class, has the attributes of a repeat object, like its start index, the number of repeat units in it,..etc


### Now how the repeat strucutre is identified (in the GenoTyper module):
  * The detection of the first unit creates a repeat object, this object stores the start index of the repeat sequence.
  * Each detected repeat unit is inserted in the repeat object, where the count of repeat units is incremented, and the end index of the repeat sequence is shifted to the index of the last inserted repeat unit.
  * Units having one mismatch with one of the repeat units are considered; to compensate expected sequencing errors, they are inserted in a buffer (temporary storage) that is flushed to the repeat object if another repeat unit exists after them
 
 ![Semantic description of image](https://github.com/hossam26644/RGT/blob/enh/pep8/ImagesForReadme/Repeatunits.png "adding repeat units")*Identifying the start index and the end index of a repeat sequence, including a unit with mismatch within the sequence structure (GTG in a CTG repeat) , green arrow represents the identified repeat sequence start index, blue arrow represents the end index of the identified repeat sequence, red curly bracket represents the sliding window identifying repeat units*
  
  * The software keeps searching for repeat units in the region downstream of the last identified repeat unit in the region downstream of the last identified repeat unit. This region starts with the base succeeding the last identified repeat unit, and the size of the region equals the user identified maximum allowed interruption length 
  ![Semantic description of image](https://github.com/hossam26644/RGT/blob/enh/pep8/ImagesForReadme/Interruption.png "dealing with interruptions")*Identifying repeat sequence with an interruption(TTTT in a CTG repeat), sliding window keeps searching for possible repeats after the final identified one, allowing a maximum interruption of length identified by the user (four or greater in this case), green arrow represents the identified repeat sequence start index, blue arrow represents the end index of the identified repeat sequence, red curly bracket represents the sliding window identifying repeat units*

Now the repeat sequence is identified, and the number of repeat units is counted in every one of them

  * Identified repeat sequences are inserted in a hashtable (a python dictionary), key is the structure and the value is the number of reads having the same structures
* Simultaneously another table is created (count table), where the key is the number of units in a structure, and the value is the number of reads having this number of units
  
  
### Now we have the structures, how the germline alleles are detected:
  * Peaks in the counts table that are above a threshold are identified, this threshold is set to the average value of read count per sequence in the counts table
    ![Semantic description of image](https://github.com/hossam26644/RGT/blob/enh/pep8/ImagesForReadme/Peaks.png "peaks in the counts table")*Counts table plotted for a DM1 sample, green line represents average value, peaks identified above the threshold value in green circles, red circle represents peaks below threshold value, peak in the red circle is discarded in this case*

  * The most abundant structures are extracted from the genotable, they are identified as structures having 30% or more reads of the most abundant structure read count
  
  * This list is matched with the peaks identified from the counts table by the repeat structures units count, then a decision tree is executed depending on the number of matches identified.

  * The decision tree: 
  
  ```


└──────────── Zero matches
│             │
│             └───  Error can't genotype
│ 
│
└──────────── >2 matches
│             │
│             └───  Error can't genotype
│
│
└──────────── 2 matches
│             │
│             └───  is the most abundant structure not matched with a peak?
│             |       │
│             |       └─── flag the sample
│             | 
│             └───  do both sequences have the same repeat units count?
|             |       |   (peak may be a result of overlapping two less abundant 
│             |       │    structures in the sample)
│             |       │
│             |       └─── flag the sample
│             |
|             └───  Export both alleles
│             
│
│
└──────────── one match
              │
              └───  Check expanded allele 
              |       | (there is a peak and there is no structure matching             
              |       │  in the most abundant structures list) 
              |       │              
              |       └─── flag the sample, export the most abundant structure 
              |            outside the list, matching the peak as the second allele
              |    
              └───  Check tandem alleles
              |       | (the n+1 structure is more abundant than the n-1  
              |       |  where n is the identified modal allele)     
              |       |              
              |       | (n+1 represents the first somatic variant, n-1 the first
              |       |  PCR slippage variant, theoretically the first is less abundant
              |       |  unless it is an allele)     
              |       |              
              |       └─── Export the n+1 structure as the second allele          
              |    
              └───  Export the identified allele as the homozygous allele of the sample              
              

```
## Plotting (graphsplotter package):
* The plotter module prepares data to be plotted, methods and variables names should be clear enough that I dont have to explain them
* table_2d_plotter module does what you think it does, it plots a table in a 2d graph, currently it is used twice, for the total counts plot and for the specified units count
* plot_3D module plots the 3d plot of 2 counts vs each other vs the abundance of samples having these counts
  * these counts are in the table_3d in the genotype object

![Semantic description of image](https://github.com/hossam26644/RGT/blob/enh/pep8/ImagesForReadme/3D_Plot.png "3d plot example") *3D plot of CAG counts vs CCG counts in an HD sample*

## Now if you have read the manual, and the code, but still have questions: you are a very unlucky person because you have to contact me
* email: hossam26644@live.com
* linkedin: [https://uk.linkedin.com/in/hossameldin-ali](https://uk.linkedin.com/in/hossameldin-ali)

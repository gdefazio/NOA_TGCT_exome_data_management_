# Exome sequencing data management and variant filtering in azoospermic and testicular germ cell tumor patients

This program executes all the passages contained in the Krausz et al. pipeline for SNV and 
INDEL filtering. It takes in input .xlsx or .csv files obtained after exome sequencing data
analysis and is aimed to find the rarest variants.

### USAGE 

```
SNV_INDEL_filter_PV.py [-h] (-s SINGLE | -d DIRECTORY) 
                       [-p PARALLEL] [-a AUX1] [-b AUX2] [-o OUTPUT]

```
Indicate a single file or a directory containing multiple files.

**OPTIONS**:

&emsp; -h, --help &emsp; show this help message and exit

&emsp; -s SINGLE, --single SINGLE &emsp; Path to a single INDEL or SNV file (debug mode)

&emsp; -d DIRECTORY, --directory DIRECTORY &emsp; Path to directory to elaborate

&emsp; -p PARALLEL, --parallel PARALLEL &emsp; Multiprocessing: Parallel degree number (default 2) unused in single mode
  
&emsp; -a AUX1, --aux1 AUX1  &emsp; Auxiliary data sheet containing SAMPLE as sample name and PHENOTYPE as 			           diagnosis or ethnic labels column. It is a CSV file format.

&emsp; &emsp; &emsp; E.g.<br> 
&emsp; &emsp; &emsp; SAMPLE,PHENOTYPE<br>
&emsp; &emsp; &emsp; A731RV,AZO<br>
&emsp; &emsp; &emsp; B732RV,AZO<br>
&emsp; &emsp; &emsp; C733RV,CRTL<br>

&emsp; -b AUX2, --aux2 AUX2 &emsp;  Auxiliary directory containing OMIM genes classification (e.g. D/r).

&emsp; -o OUTPUT, --output OUTPUT &emsp; Path to output directory. (default ./)

&emsp; -v, --verification &emsp; It verifies whether each sample in sample sheet has files (SNV, INDEL) in input directory


**Warning**:<br>
the file name must have the following format<br>
&emsp; {SAMPLE}.{SNV|INDEL}.{something_else}.{xlsx|xls|csv|csv.gz|csv.gzip}<br>
&emsp; E.g. <br>
&emsp; A731RV.SNV.FINAL.xlsx<br>

## OUTPUT description in output directory

**AF_recalculation.csv**: file containing alleles for AF recalculation
PHENOA/<br>
&emsp;	/SNV/<br>
&emsp;&emsp;		filtered SNV for each sample<br>
&emsp;	/INDEL/<br>
&emsp;&emsp;		filtered INDEL for each sample<br>
**PHENOA_INDEL_filtered.csv**: Filtered INDELs for all the samples belonging to the PHENOA phenotype;<br>
**PHENOA_SNV_filtered.csv**: Filtered SNVs for all the samples belonging to the PHENOA phenotype;<br>
**PHENOA_RECESSIVE.csv**: homozygous, putative composite Heterozygous and X-linked variants crossed against OMIM recessive genes list. HET column has the ‘pC-HET’ value indicating the putative composite Heterozygous;<br>
**PHENOA_DOMINANT.csv**: heterozygous crossed against OMIM dominant genes list;<br>
**PHENOA_noDOM_noREC.csv**: Variants of genes without an OMIM dominant|recessive annotation. HET column has the ‘pC-HET’ value indicating the putative composite Heterozygous;<br>
**PHENOA_genes4GO.txt**: Non redundant list of dominant and recessive genes contained in PHENOA_RECESSIVE.csv and PHENOA_DOMINANT.csv useful for the Gene Ontology analysis;<br>


## Rehearsal
You can test the program by using at least 50 samples.

First, try to check sample sheet with the following command:
```
SNV_INDEL_filter_PV.py\
	-d ./test/esomi_prova\  
	-a ./test/esomi_prova_data_sheet.csv\ 
	-b ./script/OMIM_30_07_22/\ 
	-v
```

Then, if the data sheet is correctly checked you can launch

```
SNV_INDEL_filter_PV.py\
	-d ./test/esomi_prova\  
	-p 20
	-a ./test/esomi_prova_data_sheet.csv\ 
	-b ./script/OMIM_30_07_22/\ 
	-o ./test/esomi_prova_test
```

After rehearsal, you should observe the obtained data to find problems, if present.

Then you can execute the script for the analysis of the entire samples ensemble.


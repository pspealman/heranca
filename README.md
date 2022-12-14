# heranca
A small script for reducing false positive variant calls using lineage information. 
(v0.7)

### Overview
  False positive variant calls can add extra labor to the task of making sense of a genotype. This tool uses lineage information supplied from the user to identify instances where a variant breaks from the infinite site assumption (ISA) and is more likely a false positive than a real variant. By defualt this is if it occurs in three or more unrelated strains.
  ![simple_lineages](https://user-images.githubusercontent.com/32845376/184647554-1ed94eeb-20b7-4978-a88a-99e67f8aa2ac.png)

  Using the figure above a candidate variant that shows up in S3, S4, and S6 is unlikely to be real as these strains lack a common ancestor that also has the same variant. Similarly, if the same variant was found in the ancestor (Anc) then it would be reasonable to assume it was likely to be real. 

  Additionally, heranca, will flag variants that continue a homopolymer run (default of 5 nucelotides of the same base over a window with a default size of 7). 
  
### Requirements
  Heranca requires the pairwise2 package from biopython.
  
### Note for NYU Greene HPC users. 
  A modified main.nf file (and example config file) of Mohammed Khalfan's GATK4 nextflow script ([original here](https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)) have been added to the ```tools``` directory. These have been modified for use on haploid yeast strains with two important changes: 
* gatk HaplotypeCaller --sample-ploidy 1
* gatk VariantFiltration -filter-name "MQ_filter" -filter "MQ < 30.0"
  
  Modified config files can be defined using the ```-C``` argument
```
nextflow -C DGY1657_A_40_nextflow.config run main.nf
```

### Lineage Metadata file format
```
Anc	demo/Anc.vcf	Anc
S1	demo/S1.vcf	Anc
S2	demo/S2.vcf	Anc
S3	demo/S3.vcf	Anc
S4	demo/S4.vcf	Anc
S5	demo/S5.vcf	Anc
S6	demo/S6.vcf	Anc, S5
```
  To import lineage data into heranca the user needs to supply a tab delimited file containing a *Name* column (unique), a *File_path* column to the vcf file, and an *Lineage* column containing the ancestors *Names* of the strain. 

### Commands 
  #### Input / Output 
  Metadata file
  ```
  '-i', '--input_metadata_file' 
  default 'demo/vcf_metadata.txt'
  ```
  Genome reference file 
  ```
  '-fa', '--fasta_file'
  default 'demo/s288c_chr3.fa'
  ```
  Export Annotation - If a VCF file has been previously annotated these annotations can be exported in a summary file ending with '_heranca.anno.tab'
  ```
  '-anno', '--export_annotation', type=bool, default = False
  ```
  
  Verbose - In addition to outputting the filter reasons additional information, such as the number of nucleotides or sizes, will also be reported
  ```
  '-v', '--verbose', type=bool, default = False
  ```
  
  Remove filtered variants - When set to True (default) variants that are filtered by any process will not be written out to the edited VCF.
  ```
  '-filter', '--remove_filtered', type=bool, default = True
  ``` 
  
  Label Inherited - When set to True (default) variants that are inherited from an ancestor will be labelled as INHERITED. NB: If also used with *'--remove_filtered'* this will result in the inherited variants being filtered.
  ```
  '-label', '--label_inherited', type=bool, default = True
  ``` 
  
  Output 
  ```
  '-o',"--output_path"
  default = current working directory
  ```
  
  Rename output - When set to True (not default) the Name in the *Name* column of the metadata file will replace the file stem in the file output name   
  ```
  '-rename',"--rename"
  default = False
  ```
  
  #### Parameters
  Window size - size in nucleotides upstream and downstream of the candidate variant 
  ```
  '-w', '--window'
  default = 7
  ```
  Max Polynucleotide - maximum length in nucleotides of a homopolymer / polynucleotide run. Candidates equal to or above this will be flagged with 'QH_Filter_exceeds_polyn_{n}', where '{n}' is the size of the polynucleotide run.
  ```
  '-p', '--max_polynucleotide'
  default = 5
  ```
  Max ISA - Maximum ISA - maximum number of unrelated strains allowable before a variant is filtered. Candidates with equal to or greater will be flagged with 'QH_fails_ISA_{n}', where '{n}' is the number of unrelated strains the variant is found in. 
  ```
  '-m', '--max_isa'
  default = 3
  ```
  Evaluate Indel - Indels, unlike SNPs have an additional complication where the alternative sequence is more easily described as a continuation of the reference sequence, such that the indel is equivalent to the refernce with no modifier. Candidates with 'ALT' sequences equal to the 'REF' sequence will be marked 'QH_Filter_low_quality_insertion' or 'QH_Filter_low_quality_deletion'.
  ```
  '-no_indel', '--no_evaluate_indel', type=bool, default = False)
  ```
  Indel flanking search region - The determination of the nucleotide where a indel starts is less precise than for SNPs. To determine if an indel has been called in a region in another strain a flanking region around each indel is searched for additional hits. By default this is 7 nucleotides up and down from the called start. 
  ```
  '-f', '--flank', type=int, default = 7)
  ```
  Percent Alignment Match - Indels can be called with nearly identical sequences. For indels less than the flanking size only 1 mismatch is allowed. For indels longer than the flanking size that but less than ten it is (n-1)/n (eg 7/8, 8/9, 9/10). For any indels greater than ten the match is 80% (0.8) by default but can be manually changed.   
  ```
  '-pct', '--pct_align', type=float, default = 0.8)
  ```
 
### Output:

Output will be a modified VCF file inluding all variants present in the original, but with those candidates that failed the polynucelotide or ISA criteria flagged as described above. 

##### Modified VCF files
* Original: 
In the original VCF file four variants are marked with **PASS** and would be evaluated by downstream analysis
```
XIII	680936	.	T	C	1282.06	PASS	AC=2;AF=1.00;AN=2;DP=30;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=34.37;SOR=1.143	GT:AD:DP:GQ:PL	1/1:0,30:30:90:1296,90,0
XIII	680940	.	C	T	1226.06	PASS	AC=2;AF=1.00;AN=2;DP=29;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=28.70;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,28:28:84:1240,84,0
XIII	809198	.	A	G	2052.06	PASS	AC=2;AF=1.00;AN=2;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.78;SOR=0.846	GT:AD:DP:GQ:PL	1/1:0,54:54:99:2066,162,0
XIII	908174	.	G	T	567.64	PASS	AC=1;AF=0.500;AN=2;BaseQRankSum=0.372;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=55.71;MQRankSum=-2.084;QD=11.58;ReadPosRankSum=-3.264;SOR=0.681	GT:AD:DP:GQ:PL	0/1:32,17:49:99:575,0,1292
```
* Output (using default ```remove_filtered = True```):
Heranca finds two of theses variants fail and omits these records
```
XIII	680936	.	T	C	1282.06	PASS	AC=2;AF=1.00;AN=2;DP=30;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=34.37;SOR=1.143	GT:AD:DP:GQ:PL	1/1:0,30:30:90:1296,90,0
XIII	809198	.	A	G	2052.06	PASS	AC=2;AF=1.00;AN=2;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.78;SOR=0.846	GT:AD:DP:GQ:PL	1/1:0,54:54:99:2066,162,0
```  

* Output (using ```remove_filtered = False```):
Setting Heranca to include filtered reads, we find rows 2 (XIII:680940) and 4 (XIII:908174)
```
XIII    680936  .   T   C   1282.06 PASS    AC=2;AF=1.00;AN=2;DP=30;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=34.37;SOR=1.143    GT:AD:DP:GQ:PL  1/1:0,30:30:90:1296,90,0
XIII    680940  .   C   T   1226.06 QH_Filter_exceeds_polyn_5   AC=2;AF=1.00;AN=2;DP=29;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=28.70;SOR=1.179    GT:AD:DP:GQ:PL  1/1:0,28:28:84:1240,84,0
XIII    809198  .   A   G   2052.06 PASS    AC=2;AF=1.00;AN=2;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.78;SOR=0.846    GT:AD:DP:GQ:PL  1/1:0,54:54:99:2066,162,0
XIII    908174  .   G   T   567.64  QH_fails_ISA_5  AC=1;AF=0.500;AN=2;BaseQRankSum=0.372;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=55.71;MQRankSum=-2.084;QD=11.58;ReadPosRankSum=-3.264;SOR=0.681    GT:AD:DP:GQ:PL  0/1:32,17:49:99:575,0,1292
```

##### Summary log files
Heranca will also generate summary log files of the counts of variant filter status categories. A complete description of the 'QH' filters is available [here](https://github.com/pspealman/heranca)

* Example of Ancestral log file:
```
MQ_FILTER	4
PASS	55
QH_exceeds_polyn	23
QH_low_genotype_relative_likelihood	6
SOR_FILTER	4
```
Here we see variants filtered by GATK4 include 4 filtered for MQ, and SOR criteria. 23 are filtered by heranca for exceeding the polyN criteria and 6 fail to meet the relative genotype likelihood criteria. Finally, 55 variants meet all the criteria. Notably, no variants are filtered for violating the ISA criteria as this is a ancestral strain.

* Example of Descendant log file:
```
PASS	27
QH_exceeds_polyn	22
QH_indel_fails_ISA	15
QH_low_genotype_relative_likelihood	1
```
Here we see that no variants were filtered by GATK4. 22 are filtered by heranca for exceeding the polyN criteria and 1 fails to meet the relative genotype likelihood criteria. 15 variants are filtered for violating ISA. Only, 27 variants meet all the criteria. 





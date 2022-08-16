# heranca
A small script for reducing false positive variant calls using lineage information

### Overview
  False positive variant calls can add extra labor to the task of making sense of a genotype. This tool uses lineage information supplied from the user to identify instances where a variant breaks from the infinite site assumption (ISA) and is more likely a false positive than a real variant. By defualt this is if it occurs in three or more unrelated strains.
  ![simple_lineages](https://user-images.githubusercontent.com/32845376/184647554-1ed94eeb-20b7-4978-a88a-99e67f8aa2ac.png)


  Using the figure above a candidate variant that shows up in S3, S4, and S6 is unlikely to be real as these strains lack a common ancestor that also has the same variant. Similarly, if the same variant was found in the ancestor (Anc) then it would be reasonable to assume it was likely to be real. 

   Additionally, heranca, will flag variants that continue a homopolymer run (default of 5 nucelotides of the same base over a window with a default size of 7). 

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
  Output 
  ```
  '-o',"--output_path"
  default = 'heranca_'
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
### Output:

Output will be a modified VCF file inluding all variants present in the original, but with those candidates that failed the polynucelotide or ISA criteria flagged as described above. 

Original: 
```
XIII	680936	.	T	C	1282.06	PASS	AC=2;AF=1.00;AN=2;DP=30;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=34.37;SOR=1.143	GT:AD:DP:GQ:PL	1/1:0,30:30:90:1296,90,0
XIII	680940	.	C	T	1226.06	PASS	AC=2;AF=1.00;AN=2;DP=29;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=28.70;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,28:28:84:1240,84,0
XIII	809198	.	A	G	2052.06	PASS	AC=2;AF=1.00;AN=2;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.78;SOR=0.846	GT:AD:DP:GQ:PL	1/1:0,54:54:99:2066,162,0
XIII	908174	.	G	T	567.64	PASS	AC=1;AF=0.500;AN=2;BaseQRankSum=0.372;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=55.71;MQRankSum=-2.084;QD=11.58;ReadPosRankSum=-3.264;SOR=0.681	GT:AD:DP:GQ:PL	0/1:32,17:49:99:575,0,1292
```
Output:
```
XIII	680936	.	T	C	1282.06	PASS	AC=2;AF=1.00;AN=2;DP=30;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=34.37;SOR=1.143	GT:AD:DP:GQ:PL	1/1:0,30:30:90:1296,90,0
XIII	680940	.	C	T	1226.06	QH_Filter_exceeds_polyn_5	AC=2;AF=1.00;AN=2;DP=29;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=28.70;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,28:28:84:1240,84,0
XIII	809198	.	A	G	2052.06	PASS	AC=2;AF=1.00;AN=2;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=32.78;SOR=0.846	GT:AD:DP:GQ:PL	1/1:0,54:54:99:2066,162,0
XIII	908174	.	G	T	567.64	QH_fails_ISA_5	AC=1;AF=0.500;AN=2;BaseQRankSum=0.372;DP=57;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=55.71;MQRankSum=-2.084;QD=11.58;ReadPosRankSum=-3.264;SOR=0.681	GT:AD:DP:GQ:PL	0/1:32,17:49:99:575,0,1292
```  

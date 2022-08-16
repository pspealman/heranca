# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:01:51 2022

heranca / herança - inheritance

- qc vcfs with both sequence and lineage considerations 

ver 0.1 - beta (hide elapse)
ver 0.2 - beta (printer ambiguous)
ver 0.3 - public (closed suspicion)
    _x_ enable non-"PASS" variants to propagate to the final file
    _x_ enabled lineage counting

@author: pspealman

# Example format for lineage metadata file (tab separated):

Anc	demo/Anc.vcf	Anc
S1	demo/S1.vcf	Anc
S2	demo/S2.vcf	Anc
S3	demo/S3.vcf	Anc
S4	demo/S4.vcf	Anc
S5	demo/S5.vcf	Anc
S6	demo/S6.vcf	Anc, S5

"""

import argparse

parser = argparse.ArgumentParser()

#io
parser.add_argument('-i', '--input_metadata_file', nargs='?', type=str, 
                    default = 'C:/Gresham/Project_Carolino_new/combine_vcf_metadata_indels.txt')
parser.add_argument('-fa', '--fasta_file', nargs='?', type=str, 
                    default = 'C:/Gresham/genomes/NCBI_R64/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa')
parser.add_argument('-o',"--output_path", nargs='?', type=str,
                    default = 'heranca_')
#parameters
parser.add_argument('-w', '--window', nargs='?', type=int, default = 7)
parser.add_argument('-p', '--max_polynucleotide', nargs='?', type=int, default = 5)
parser.add_argument('-m', '--max_isa', nargs='?', type=int, default = 3)

args = parser.parse_args()

genome_sequence_dict = {}
universal_vc = {}
locus_to_vc = {}
ancestor_lookup = {}
strain_variant_catalog = {}
metadata_dict = {}

def load_genome():
    global genome_sequence_dict
    
    fasta_file = open(args.fasta_file)
    
    abrir = False
    seq = ''
    name = ''
    
    for line in fasta_file:
        if line[0] == '>':
            
            if abrir:
                genome_sequence_dict[name] = seq
                abrir = False
            
            if not abrir:    
                print(line)
                abrir = True
                name = line.split('>')[1].split(' ')[0]
                seq = ''           
                
        else:
            seq += line.strip()
          
    fasta_file.close()
    
    if abrir:
        genome_sequence_dict[name] = seq
        abrir = False
        

def first_pass_vcf(strain, filename):
    global strain_variant_catalog
    
    if strain not in strain_variant_catalog:
        strain_variant_catalog[strain] = {}
                
    infile = open(filename)
    
    for line in infile:
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG7CVAFXX_mini02_40
    #I	2146	.	G	A	63.64	QD_filter	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.105;DP=46;ExcessHet=3.0103;FS=1.855;MLEAC=1;MLEAF=0.500;MQ=42.75;MQRankSum=-2.856;QD=1.45;ReadPosRankSum=-1.956;SOR=1.127	GT:AD:DP:GQ:PL	0/1:38,6:44:71:71,0,1199
    #I	25340	.	C	A	1173.06	PASS	AC=2;AF=1.00;AN=2;DP=30;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=40.48;QD=28.73;SOR=1.255	GT:AD:DP:GQ:PL	1/1:0,29:29:87:1187,87,0
        if line[0] != '#':
            filter_value = line.split('\t')[6]
            #if line.split('\t')[6] == 'PASS':
            chromo = line.split('\t')[0]
            #pos = int(line.split('\t')[1])
            start = int(line.split('\t')[1])
            ref = line.split('\t')[3]
            alt = line.split('\t')[4]
            
            variant = ('{chromo},{start},{ref},{alt}').format(
                chromo = chromo, start = start, ref = ref, alt = alt)
            
            strain_variant_catalog[strain][variant] = filter_value

    infile.close()
    
def get_strain_count(strain, variant, metadata_dict):
    global strain_variant_catalog
    strain_ct = set()   
    
    
    strain_ancestor_set = metadata_dict[strain]['ancestor_set']
        
    
    if strain not in strain_ancestor_set: 
        for anc_strain in strain_ancestor_set:
            if variant in strain_variant_catalog[anc_strain]:
                return(0)
        
        for other_strain in strain_variant_catalog:
            if other_strain != strain:
                if other_strain not in strain_ancestor_set: 
                    if variant in strain_variant_catalog[other_strain]:
                        strain_ct.add(other_strain)
                        
    return(len(strain_ct))

def parse_vcf(strain, filename, outfile_uid, metadata_dict):
    global universal_vc, locus_to_vc, genome_sequence_dict, strain_variant_catalog
    
    
                
    infile = open(filename)
    print(filename)
    edited_filename = ('{}_heranca.vcf').format(filename.split('.vcf')[0])
    edited_file = open(edited_filename, 'w')
    
    for line in infile:
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG7CVAFXX_mini02_40
    #I	2146	.	G	A	63.64	QD_filter	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.105;DP=46;ExcessHet=3.0103;FS=1.855;MLEAC=1;MLEAF=0.500;MQ=42.75;MQRankSum=-2.856;QD=1.45;ReadPosRankSum=-1.956;SOR=1.127	GT:AD:DP:GQ:PL	0/1:38,6:44:71:71,0,1199
    #I	25340	.	C	A	1173.06	PASS	AC=2;AF=1.00;AN=2;DP=30;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=40.48;QD=28.73;SOR=1.255	GT:AD:DP:GQ:PL	1/1:0,29:29:87:1187,87,0
        if line[0] == '#':
            edited_file.write(line)
            
        if line[0] != '#':
            original_line = line
            filter_value = line.split('\t')[6].upper()
            

            chromo = line.split('\t')[0]

            start = int(line.split('\t')[1])
            ref = line.split('\t')[3]
            alt = line.split('\t')[4]
            
            variant = ('{chromo},{start},{ref},{alt}').format(
                chromo = chromo, start = start, ref = ref, alt = alt)

            strain_ct = get_strain_count(strain, variant, metadata_dict)
                        
            if strain_ct >= args.max_isa:
                filter_value = ("QH_fails_ISA_{}").format(strain_ct)
                                
            if filter_value == "PASS":
                filter_value = eval_poly(chromo, start)
             
            if filter_value == "PASS":
                if len(ref) != len(alt):
                    filter_value = eval_indel(chromo, start, 0.80, alt, ref)
                else:
                    filter_value = eval_snp(chromo, start, alt)
                    
            print('original_line', original_line)
            print('filter_value', filter_value)
            line_new = original_line.replace("PASS", filter_value)
            edited_file.write(line_new)

    infile.close()
    edited_file.close()
    
        
def count_max_poly(seq):
    seq = seq.lower()
    ct = 0
    base = ''
    if len(seq) > 0:
        for nuc in ['a','t','c','g','n','u']:  
            nuc_ct = seq.count(nuc)
            if nuc_ct > ct:
                ct = nuc_ct
                base = nuc
    return(ct, base)

def eval_poly(chromo, start):
    result = "PASS"
    
    pre_seq = genome_sequence_dict[chromo][start-(args.window)-1:start-1]
    pre_ct, _base = count_max_poly(pre_seq)
    
    post_seq = genome_sequence_dict[chromo][start:start+(args.window)+1]
    post_ct, _base = count_max_poly(post_seq)
    
    if pre_ct >= args.max_polynucleotide or post_ct >= args.max_polynucleotide:
        result = ("QH_Filter_exceeds_polyn_{}").format(max(pre_ct,post_ct))
        
    return(result)
    

def eval_indel(chromo, start, cutoff, alt, ref):
    result = "PASS"
        
    length = len(alt)-len(ref)
    
    #'insertion'
    if length >= 2:
        stop = start + (len(alt)-len(ref)) -1 
        insertion = alt[len(ref):]
        seq = genome_sequence_dict[chromo][(start - len(ref) -1):stop+1]
        
        if insertion in seq:
            result = "QH_Filter_low_quality_insertion"
            return(result)
            
    #'deletion'
    if length <= -2:
        skip = start + (len(ref)-len(alt))
        stop = start + 2*(len(ref)-len(alt))
        deletion = ref[(len(alt)-len(ref)):]
        seq = genome_sequence_dict[chromo][skip:stop]
                
        if deletion in seq:
            result = "QH_Filter_low_quality_deletion"
            return(result)
             
    return(result)

def poly_run(seq, alt, runmode):
    print('seq, alt, runmode', seq, alt, runmode)
    
    result = 'PASS'
    
    for ct in range(1,11):        
        if runmode == 'before':
            field = seq[(-1*ct):]+alt

        if runmode == 'after':
            field = seq[:(ct)-1]+alt

        hit = field.count(alt)
        search_string = ct*alt
                        
        if ct >= 10:
            return(result)
                        
        if (ct >= 8):
            if (hit/ct) >= 0.8:
                result = ('QH_low_complexity_region_{}').format(hit)
                return(result)
            
        if (ct > 4):
            if search_string in field:
                result = ('QH_forms_polyN_{}').format(len(search_string))
                return(result)
         
def eval_snp(chromo, start, alt):    
    result = "PASS"
            
    
    seq = genome_sequence_dict[chromo][start-(args.window)-1:start-1]
    seq = seq.lower()
    alt = alt.lower()
    
    result = poly_run(seq, alt, 'before')
    
    if result == "PASS":
        
        seq = genome_sequence_dict[chromo][start:start+(args.window)+1]
        seq = seq.lower()
        alt = alt.lower()
        
        result = poly_run(seq, alt, 'after')
                
    return(result)
        
def load_metadata():
    global ancestor_lookup, metadata_dict
    
    input_metadata_file = open(args.input_metadata_file)
    outfile_uid = open(args.output_path + '_uid.tab','w')
    
    for line in input_metadata_file:
        if line[0] != '#':
            print(line)
            line = line.strip()
            strain, filename, ancestor = line.split('\t')
            
            ancestor_set = set()
            
            if ',' in ancestor:
                for each in list(ancestor.split(',')):
                    ancestor_set.add(each.strip())
            else:
                ancestor_set.add(ancestor.strip())
                
            ancestor_lookup[strain] = ancestor_set
            
            if strain in metadata_dict:
                print('Error: Duplicate Strain Name ...', strain)
                1/0
                
            metadata_dict[strain] = {'filename': filename,
                                     'ancestor_set':ancestor_set}
            
            first_pass_vcf(strain, filename)
    
    for strain in metadata_dict:
        filename = metadata_dict[strain]['filename']
        
        parse_vcf(strain, filename, outfile_uid, metadata_dict)
        
    outfile_uid.close()
                    
#body    
load_genome() 
load_metadata()

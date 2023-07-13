# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:01:51 2022

heranca / herança - inheritance

- qc vcfs with both sequence and lineage considerations 

ver 0.7 - public (Cooperation Fisherman)
    _x_ added file output rename operation  
    _x_ add 'inherited' label
ver 0.8 - public (dispelvels)
    _x_ added logic filters for output
    

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
import os
import pathlib
import argparse
from Bio import pairwise2


parser = argparse.ArgumentParser()
#io
parser.add_argument('-i', '--input_metadata_file', nargs='?', type=str,
                    default = 'demo/vcf_metadata.txt')
                    #default = 'C:/Gresham/Project_Carolino_new/combine_vcf_metadata_indels.txt')
                    #default = 'C:/Gresham/Project_Carolino_new/combine_vcf_metadata_snps.txt')

parser.add_argument('-fa', '--fasta_file', nargs='?', type=str,
                    default = 'demo/s288c_chr3.fa')
                    #default = 'C:/Gresham/genomes/Ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa')

parser.add_argument('-o',"--output_path",
                    default = 'demo/heranca_')
                    #default = 'C:\\Gresham\\Project_Carolino_new\\vcf\\heranca\\heranca_')
                    #default = 'C:/Gresham/Project_Carolino_new/supplemental/vcf/vcf_files_filtered/')

#parameters
parser.add_argument('-log', '--enable_log', action='store_false')
parser.add_argument('-anno', '--export_annotation', action='store_true')
parser.add_argument('-no_indel', '--no_evaluate_indel', action='store_true')
parser.add_argument('-w', '--window', nargs='?', type=int, default = 7)
parser.add_argument('-p', '--max_polynucleotide', nargs='?', type=int, default = 5)
parser.add_argument('-m', '--max_isa', nargs='?', type=int, default = 3)
parser.add_argument('-pct', '--pct_align', type=float, default = 0.8)
parser.add_argument('-f', '--flank', type=int, default = 7)
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-rename', '--rename', action='store_false')
parser.add_argument('-filter', '--remove_filtered', action='store_false')
parser.add_argument('-inherit', '--label_inherited', action='store_false')

args = parser.parse_args()

genome_sequence_dict = {}
universal_vc = {}
locus_to_vc = {}
ancestor_lookup = {}
strain_variant_catalog = {}
metadata_dict = {}
indel_dict = {}

if args.output_path:
    output_path = args.output_path
else:
    output_path = os.getcwd()

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
                #print(line)
                line=line.strip()
                abrir = True
                name = line.split('>')[1].split(' ')[0]
                seq = ''           
                
        else:
            seq += line.strip()
          
    fasta_file.close()
    
    if abrir:
        genome_sequence_dict[name] = seq
        abrir = False
        
        
def build_indel_set(strain, chromo, start, ref, alt):
    global indel_dict
    
    indow = args.flank
    
    if len(ref) > len(alt):
        istype = 'del'
    else:
        istype = 'ins'
        
    if istype not in indel_dict:
        indel_dict[istype]={}
    
    if chromo not in indel_dict[istype]:
        indel_dict[istype][chromo] = {}
        
    for nt in range(start - indow, start + indow + 1):
        if nt not in indel_dict[istype][chromo]:
            indel_dict[istype][chromo][nt] = {}
            
        if strain not in indel_dict[istype][chromo][nt]:
            indel_dict[istype][chromo][nt][strain] = {}
            
        if ref not in indel_dict[istype][chromo][nt][strain]:
            indel_dict[istype][chromo][nt][strain][ref] = set()
        
        indel_dict[istype][chromo][nt][strain][ref].add(alt)
                            
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
            start = int(line.split('\t')[1])
            ref = line.split('\t')[3]
            alt = line.split('\t')[4]
            
            variant = ('{chromo},{start},{ref},{alt}').format(
                chromo = chromo, start = start, ref = ref, alt = alt)
            
            strain_variant_catalog[strain][variant] = filter_value
            
            if len(ref) != len(alt):
                build_indel_set(strain, chromo, start, ref, alt)

    infile.close()
    
def get_strain_count(strain, variant, metadata_dict):
    global strain_variant_catalog
    strain_ct = set()
    
    strain_ancestor_set = metadata_dict[strain]['ancestor_set']
    
    if strain not in strain_ancestor_set: 
        for anc_strain in strain_ancestor_set:
            if variant in strain_variant_catalog[anc_strain]:
                return(True, 0)
        
        for other_strain in strain_variant_catalog:
            if other_strain != strain:
                if other_strain not in strain_ancestor_set: 
                    if variant in strain_variant_catalog[other_strain]:
                        strain_ct.add(other_strain)
                        
    return(False, len(strain_ct))

def eval_minimum_relative_likelihood(values):
    '''
    minimum relative likelihood > 100: 
        (Genotype_Likelihood * Genotype_ALT_reads) / (Depth - Genotype_ALT_reads + 1)
    '''
    #1:0,7:7:99:262,0
    # Note: this only supports haploid currently
    #   Diploid example:
    #      GT:AD:DP:GQ:PL	0/1:38,6:44:71:71,0,1199
    #   Haploid example:
    #      GT:AD:DP:GQ:PL   1:1,53:54:99:1413,0
    
    if '/' not in values:
        gl = int(values.split(':')[4].split(',')[0])
        ga = int(values.split(':')[1].split(',')[1])
        dp = int(values.split(':')[2])
        
        mrl = (gl * ga) / max((dp - ga + 1),1)
        
        if mrl < 100:
            if args.verbose:
                outline = ('QH_low_genotype_relative_likelihood_{}').format(mrl)
            else:
                outline = ('QH_low_genotype_relative_likelihood')
                
            return(outline)
    
    return('PASS')
    
    

def parse_vcf(strain, filename, outfile_uid, metadata_dict):
    global universal_vc, locus_to_vc, genome_sequence_dict, strain_variant_catalog
                           
    log_dict = {}
    
    infile = open(filename)
    
    file_stem = pathlib.Path(filename).stem
        
    if args.rename:
        file_stem = strain

    edited_filename = ('{}{}_heranca.vcf').format(output_path, file_stem)
    edited_file = open(edited_filename, 'w')

    print('Edited vcf file saved to: ', edited_filename)
    
    if args.export_annotation:
        anno_filename = ('{}{}_heranca.anno.tab').format(output_path, file_stem)
        anno_file = open(anno_filename, 'w')
        print('Annotation summary saved to: ', anno_filename)
            
    for line in infile:
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG7CVAFXX_mini02_40
    #II	389428	.	T	G	1205.04	PASS	AC=1;AF=1.00;AN=1;DP=27;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=60.00;QD=23.23;SOR=0.765	GT:AD:DP:GQ:PL	1:0,27:27:99:1215,0
    #III	143131	.	T	C	1403.04	PASS	AC=1;AF=1.00;AN=1;BaseQRankSum=-1.551;DP=54;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=60.00;MQRankSum=0.000;QD=25.98;ReadPosRankSum=-1.030;SOR=0.238	GT:AD:DP:GQ:PL	1:1,53:54:99:1413,0
        if line[0] == '#':
            edited_file.write(line)
            
        if line[0] != '#':
            original_line = line
            chromo = line.split('\t')[0]

            start = int(line.split('\t')[1])
            ref = line.split('\t')[3]
            alt = line.split('\t')[4]
            
            filter_value = line.split('\t')[6].upper()
            
            if len(ref) != len(alt):
                var_type = 'indel'
            else:
                var_type = 'snp'
            
            variant = ('{chromo},{start},{ref},{alt}').format(
                chromo = chromo, start = start, ref = ref, alt = alt)
            
            if var_type == 'snp':
                inherited, strain_ct = get_strain_count(strain, variant, metadata_dict)
                
                if args.label_inherited:
                    if inherited:
                        filter_value = 'INHERITED'
                
                if not inherited:
                    if strain_ct >= args.max_isa:
                        if args.verbose:
                            filter_value = ("QH_fails_ISA{}").format(strain_ct)
                        else:
                            filter_value = ("QH_fails_ISA")
                
                if filter_value == "PASS":
                    filter_value = eval_snp(chromo, start, alt)
            
            if var_type == 'indel' and (not args.no_evaluate_indel):
                filter_value = eval_indel(strain, chromo, start, 0.80, alt, ref)
            
            if filter_value == "PASS":
                filter_value = eval_poly(chromo, start)
                    
            if filter_value == "PASS":
                if 'GT:AD:DP:GQ:PL' in line:
                    values = line.split('\t')[9].strip()
                    filter_value = eval_minimum_relative_likelihood(values)
                                
            if filter_value not in log_dict:
                log_dict[filter_value] = 0
            log_dict[filter_value] += 1
                                    
            line_new = original_line.replace("PASS", filter_value)
            
            if args.remove_filtered:
                if filter_value == 'PASS':
                    edited_file.write(line_new)
            else:
                edited_file.write(line_new)
            
            if args.export_annotation:
                if 'ANN=' in line.split('\t')[7].upper():
                    #anno_deets = line.split('\t')[7].split('ANN=')[1]
                    anno_set = set(line.split('\t')[7].split('ANN=')[1].split(','))
                    
                    for anno in anno_set:
                        annotation, impact, _genename, gene, feature = anno.split('|')[1:6]
                        #print(annotation, impact, _genename, gene, feature)
                        outline = ('{chromo}\t{start}\t{ref}\t{alt}\t'
                                   '{annotation}\t{impact}\t{gene}\t'
                                   '{feature}\t{filter_value}\n').format(
                                       chromo = chromo, start = start, ref = ref, alt = alt,
                                       annotation = annotation, impact = impact, gene = gene, 
                                       feature = feature, filter_value = filter_value)
                    
                        anno_file.write(outline)
                
    infile.close()
    edited_file.close()
    
    if args.export_annotation:
        anno_file.close()
        
    if args.enable_log:
        enable_filename = ('{}{}_heranca.log').format(output_path, file_stem)
        print('Log file saved to: ', enable_filename)
        enable_file = open(enable_filename, 'w')
        
        values = list(log_dict.keys())
        values.sort()
        
        for value in values:
            outline = ('{}\t{}\n').format(value, log_dict[value])
            enable_file.write(outline)
            print(outline.strip())
            
        enable_file.close()
                    
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
        if args.verbose:
            result = ("QH_exceeds_polyn_{}").format(max(pre_ct,post_ct))
        else:
            result = ("QH_exceeds_polyn")
        
    return(result)

def compare_seq(seq1, seq2):
    
    run_length = max(len(seq1), len(seq2))
    
    if run_length >= 10:
        pct_align_threshold = args.pct_align
    if (run_length <= 10):
        pct_align_threshold = (run_length-1)/run_length
    if (run_length <= args.flank):
        pct_align_threshold = 1
        
    #print('Running pairwise ...', seq1, seq2)
    
    for align in pairwise2.align.globalxx(seq1, seq2):
        pct_align = ((align.score)/run_length)
        
        if pct_align >= pct_align_threshold:
            return(True)
        
    return(False)

def check_indel_anc(strain, ancestor_set, istype, chromo, nt):
    #Check if sequence is in ancestor:
    #in_anc = False
    for anc_strain in ancestor_set:
        #print(istype, chromo, nt, strain)
        ref_set = indel_dict[istype][chromo][nt][strain]
        
        if anc_strain in indel_dict[istype][chromo][nt]:
            anc_ref_set = indel_dict[istype][chromo][nt][anc_strain]
            
            for ref, alt_set in ref_set.items():
                for anc_ref, anc_alt_set in anc_ref_set.items():
                    #If the ref is in an ancestor then ...
                    if compare_seq(ref, anc_ref):
                        #If the alt is in an ancestor then ...
                        for alt in alt_set:
                            for anc_alt in anc_alt_set:
                                if compare_seq(alt, anc_alt):      
                                    return(True)
    
    #print('in_anc', strain, ancestor_set,  istype, chromo, nt)
    return(False)

def eval_indel(strain, chromo, start, cutoff, alt, ref):
    
    global indel_dict, ancestor_lookup
    
    indow = args.flank
    
    ancestor_set = ancestor_lookup[strain]
    
    strain_set = set()
    registered_in_anc = False
               
    if strain not in ancestor_set:
    
        if len(ref) > len(alt):
            istype = 'del'
        else:
            istype = 'ins'
    
        for nt in range(start - indow, start + indow +1):
            if nt in indel_dict[istype][chromo]:
                other_strains_set = indel_dict[istype][chromo][nt]

                #Check if sequence is in ancestor:
                in_anc = check_indel_anc(strain, ancestor_set, istype, chromo, nt)
                                
                #If even matches once in ancestor - never mind the other strains
                if in_anc:
                    registered_in_anc = True
                
                #If it isn't in an ancestor then ...
                if not in_anc and not registered_in_anc:
                    for other_strain in other_strains_set: 
                        if (other_strain != strain):
                            ref_set = indel_dict[istype][chromo][nt][strain]
                            other_ref_set = indel_dict[istype][chromo][nt][other_strain]
                                                  
                            for ref, alt_set in ref_set.items():
                                for other_ref, other_alt_set in other_ref_set.items():
                                    #If the ref is in the other then ...          
                                    if compare_seq(ref, other_ref):
                                        #If the alt is the other then ... 
                                        for alt in alt_set:
                                            for other_alt in other_alt_set:
                                                if compare_seq(alt, other_alt):   
                                                    strain_set.add(other_strain)

    result = "PASS"
    if args.label_inherited:
        if registered_in_anc:
            result = 'INHERITED'
        
    if not registered_in_anc:
        if len(strain_set) >= args.max_isa:
            if args.verbose:
                result = ("QH_indel_fails_ISA_{}").format(len(strain_set))
            else:
                result = ("QH_indel_fails_ISA")

    return(result)

def poly_run(seq, alt, runmode):
    #print('seq, alt, runmode', seq, alt, runmode)
    
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
                if args.verbose:
                    result = ('QH_low_complexity_region_{}').format(hit)
                else:
                    result = ('QH_low_complexity_region')
                    
                return(result)
            
        if (ct > 4):
            if search_string in field:
                if args.verbose:
                    result = ('QH_exceeds_polyN_{}').format(len(search_string))
                else:
                    result = ('QH_exceeds_polyN')
                    
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
    outfile_uid = open(output_path + '_heranca_uid.log','w')
    
    for line in input_metadata_file:
        if line[0] != '#':
            #print(line)
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


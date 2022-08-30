# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 22:19:45 2022

@author: pspea
"""

import argparse

parser = argparse.ArgumentParser()
#io
parser.add_argument('-utr', '--input_utr_file', nargs='?', type=str, 
                    default = 'C:/Gresham/Project_Carolino_new/metadata/sgd_utr_results.tsv')
                    #default = 'C:/Gresham/genomes/McManus/saccharomyces_cerevisiae.gff')
parser.add_argument('-utr_source', '--input_utr_source', nargs='?', type=str, 
                    default = 'yeastmine')
                    #default = 'spealman_naik_2019')

parser.add_argument('-gff', '--input_gff_file', nargs='?', type=str, 
                    default = 'C:/Gresham/Project_Carolino_new/ensembl_50/Saccharomyces_cerevisiae.R64-1-1.50.gtf')
parser.add_argument('-o',"--output_gff_file", nargs='?', type=str,
                    default = 'C:/Gresham/Project_Carolino_new/metadata/Saccharomyces_cerevisiae.R64-1-1.50_cjm-utr-edit.gtf')

args = parser.parse_args()

#parameters


infile = open(args.input_utr_file)

orf_dict = {}
utr_dict = {}

if args.input_utr_source == 'yeastmine':
    for line in infile:
        '''
        S000156980-3prime-utr	YBR265W	chrII	739545	739677	1
        S000156980-5prime-utr	YBR265W	chrII	738559	738581	1
        '''
        orf = line.split('\t')[1]
        if orf not in orf_dict:
            orf_dict[orf] = set()
            
        istype = line.split('\t')[0].rsplit('-',2)[1]
        orf_dict[orf].add(istype)
        
        if orf not in utr_dict:
            utr_dict[orf] = set()
            
        utr_dict[orf].add(int(line.split('\t')[3]))
        utr_dict[orf].add(int(line.split('\t')[4]))    
        
if args.input_utr_source == 'spealman_naik_2019':
    for line in infile:
        '''
        chrI	AWN	five_prime_UTR	136874	136913	.	+	.	PARENT=YAL008W_mRNA
        chrI	AWN	three_prime_UTR	137511	137616	.	+	.	PARENT=YAL008W_mRNA
        '''
        
        if '_prime_UTR' in line:
            line = line.strip()
            orf = line.split('\t')[8].upper()
            
            if 'PARENT=' in orf:
                orf = orf.split('PARENT=')[1]
                
            if '_MRNA' in orf:
                orf = orf.split('_MRNA')[0]
                
            if orf not in orf_dict:
                orf_dict[orf] = set()
                
            istype = line.split('\t')[2]
            orf_dict[orf].add(istype)
            
            if orf not in utr_dict:
                utr_dict[orf] = set()
                
            utr_dict[orf].add(int(line.split('\t')[3]))
            utr_dict[orf].add(int(line.split('\t')[4]))    
        
infile.close()

count_type_dict = {}

for orf in orf_dict:
    ct_type = len(orf_dict[orf])
    
    if ct_type not in count_type_dict:
        count_type_dict[ct_type] = set()
    
    count_type_dict[ct_type].add(orf)
    
    if ct_type >= 3:
        print(orf_dict[orf])

for ct_type in count_type_dict:             
    print(ct_type, len(count_type_dict[ct_type]))
    
    
infile = open(args.input_gff_file)
gff_file = open(args.output_gff_file, 'w')

gene_sets_w_utr = set()

for line in infile:
    if line[0]=='#':
        gff_file.write(line)
        
    if line[0]!='#':
        #print(line)
        '''
        V	sgd	gene	138918	139763	.	-	.	gene_id "YEL009C"; gene_name "GCN4"; gene_source "sgd"; gene_biotype "protein_coding";
        V	sgd	transcript	138918	139763	.	-	.	gene_id "YEL009C"; transcript_id "YEL009C_mRNA"; gene_name "GCN4"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_name "GCN4"; transcript_source "sgd"; transcript_biotype "protein_coding";
        V	sgd	exon	138918	139763	.	-	.	gene_id "YEL009C"; transcript_id "YEL009C_mRNA"; exon_number "1"; gene_name "GCN4"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_name "GCN4"; transcript_source "sgd"; transcript_biotype "protein_coding"; exon_id "YEL009C_mRNA-E1";
        V	sgd	CDS	138921	139763	.	-	0	gene_id "YEL009C"; transcript_id "YEL009C_mRNA"; exon_number "1"; gene_name "GCN4"; gene_source "sgd"; gene_biotype "protein_coding"; transcript_name "GCN4"; transcript_source "sgd"; transcript_biotype "protein_coding"; protein_id "YEL009C";
        '''
        if (line.split('\t')[2] == 'gene')  or (line.split('\t')[2] == 'transcript'):
            gene_id = line.split('\t')[8].split('";')[0].split('"')[1]
            print(gene_id)
            if gene_id in utr_dict and (gene_id not in gene_sets_w_utr):
                gene_sets_w_utr.add(gene_id)
                sign = line.split('\t')[6]
                original_left = int(line.split('\t')[3])
                original_right = int(line.split('\t')[4])
                
                qleft = min(utr_dict[gene_id])
                qright = max(utr_dict[gene_id])
                            
                if qleft < original_left:
                    new_left = qleft
                    make_utr = True
                else:
                    new_left = original_left
                    
                if qright > original_right:
                    new_right = qright
                    make_utr = True
                else:
                    new_right = original_right
                
                chromo = line.split('\t')[0]
                source = line.split('\t')[1]
                feature = line.split('\t')[2]
                dot_1 = line.split('\t')[5]
                sign = line.split('\t')[6]
                dot_2 = line.split('\t')[7]
                details = line.split('\t')[8]
                #
                
                new_line = ('{chromo}\t{source}_modified\t{feature}\t'
                            '{new_left}\t{new_right}\t{dot_1}\t'
                            '{sign}\t{dot_2}\t{details}').format(
                                chromo = chromo, source = source, feature = feature,
                                new_left = new_left, new_right = new_right, dot_1 = dot_1,
                                sign = sign, dot_2 = dot_2, details = details)
                    
                gff_file.write(new_line)
                
                if make_utr:
                    if sign == '+':
                        five_prime_utr_start = new_left
                        five_prime_utr_stop = original_left
                        
                        three_prime_utr_start = new_right
                        three_prime_utr_stop = original_right
                        
                    if sign == '-':
                        three_prime_utr_start = new_left
                        three_prime_utr_stop = original_left
                        
                        five_prime_utr_start = new_right
                        five_prime_utr_stop = original_right
                        
                    five_prime_utr_line = ('{chromo}\t{source}_modified\t{feature}\t'
                                '{left}\t{right}\t{dot_1}\t'
                                '{sign}\t{dot_2}\t{details}').format(
                                    chromo = chromo, source = source, feature = 'five_prime_utr',
                                    left = five_prime_utr_start, right = five_prime_utr_stop, 
                                    dot_1 = dot_1, sign = sign, dot_2 = dot_2, details = details)
                        
                    gff_file.write(five_prime_utr_line)
                    
                    three_prime_utr_line = ('{chromo}\t{source}_modified\t{feature}\t'
                                '{left}\t{right}\t{dot_1}\t'
                                '{sign}\t{dot_2}\t{details}').format(
                                    chromo = chromo, source = source, feature = 'three_prime_utr',
                                    left = three_prime_utr_start, right = three_prime_utr_stop, 
                                    dot_1 = dot_1, sign = sign, dot_2 = dot_2, details = details)
                        
                    gff_file.write(three_prime_utr_line)
                        
        else:
            gff_file.write(line)

infile.close()
gff_file.close()
            


    
    
    
    
    
    
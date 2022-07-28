#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 09:02:17 2022

@author: jonwinkelman
"""
from Bio import SeqIO





path_to_genbank = './data/references/Reference_FERM_BP3421.gbk'
bam_filepath = './RNAseq_bam_files/BAN-1.bam'

        
        
def build_genbank_dict(path_to_genbank):
    """
    builds dict {chromosomeID:{featureID:gbk_feature_object}} from a genbank annotation file
    
    """
    chromosome_dict = {}
    gb_obj = SeqIO.parse(path_to_genbank, "genbank")
    for chromosome in gb_obj:  
        for feature in chromosome.features:
            if feature.type == 'gene':
                if chromosome_dict.get(chromosome.name): 
                    chromosome_dict[chromosome.name].update({feature.qualifiers['locus_tag'][0]:feature})
                else:
                    chromosome_dict[chromosome.name] = {feature.qualifiers['locus_tag'][0]:feature}
    return chromosome_dict
    
    

    

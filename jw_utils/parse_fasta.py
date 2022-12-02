#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 14:08:25 2022

@author: jonwinkelman

"""
from jw_utils import file_utils as fu

def get_seq_dict(path_to_fasta):
    '''
    parse fast file and return as a dictionary
    
    parameters:
        path_to_fasta (str): path to a fasta proteome
        
        return (dict): eqID as key and sequence as value
    '''
    seq = ''
    seq_dict = {}
    prot_id = None
    with open(path_to_fasta, 'r') as f:
        for line in f:
            if line[0] == '>':
                if prot_id:
                    seq_dict[prot_id] = seq
                prot_id = line.split(' ')[0][1:].strip()
                seq = ''
            else:
                seq = seq + line.strip()
        seq_dict[prot_id] = seq
            
    return seq_dict


def get_protein_subset(path_to_fasta, seq_ids):
    '''
    return a dict with seqID as key and sequence as value
    
    parameters:
        path_to_fasta (str): path to a fasta proteome
    '''
    
    seq_dict_subset = {}
    seq_dict = get_seq_dict(path_to_fasta)
    for i in seq_ids:
        seq_dict_subset[i] = seq_dict.get(i)
    return seq_dict_subset


def write_to_fasta(seq_dict, path):
    """Write a dictionary to a fasta file"""  
    with open(path, 'w') as f:
        for name, seq in seq_dict.items():
            f.write(f'>{name}\n{seq}\n')
            
            
            
def concat_mult_fastas(dir_with_fastas, write_file=True, path=None, return_dict=False):
    """Combine multiple fasta files into one dict and write to file.
    
        Arguments:
        dir_with_fastas (str):
        write_file (bool): if True, write file to input path argument
        path (str): path to the file to be written
    """
    new_dict = {}
    fps = fu.get_filepaths_in_dir(dir_with_fastas,end_to_exclude='.DS_Store')
    for path in fps:
        d = get_seq_dict(path)
        for key in d.keys():
            new_dict[key] = d[key]
    if write_file:
        write_to_fasta(new_dict, path)
    if return_dict:
        return new_dict
            
            
            
            
def fasta_to_phyl(fasta_aln_path, output_path):
    """
    turn fasta format into simple phylip format
    Note: other info in fasta will be lost in this phylip format
    """

    fasta_dict = get_seq_dict(fasta_aln_path)
    num_algnments = len(list(fasta_dict.values()))
    aln_length = len(list(fasta_dict.values())[0])
    with open(output_path, 'w') as f:
        f.write(str(num_algnments)+' ' + str(aln_length) + '\n')
        for protID, seq in fasta_dict.items():
            f.write(f'{protID} {seq}\n')
            
            
        
def get_seq_region(path_to_fasta, feature_name, start_coord, end_coord):
    """Return the sequence within given fasta element from the given coordinates.

    Arguments:
        feature name (str): protein, gene, contig ID, genome etc... (the thing after the '<')
        start_coord (int): beginning of sequence (DNA or protein)
        end_coord (int): end of sequence, DNA or protein
        *start and end coordinates are 1-based. 1 = first element, 
        1,3 will return the first-third elements of the sequences
    """
    fasta_dict = get_seq_dict(path_to_fasta)
    feature = fasta_dict.get(feature_name)
    if not feature:
        raise Exception('Feature name not in fasta dictionary')
    else:
        if start_coord <1:
            raise Exception('This is 1-based. Start coordinate must be >=1')
        if end_coord > len(feature):
            raise Exception(f'end coordinate is > length {len(feature)} of {feature_name}')
        start_coord = start_coord-1 # to make 1-based
        return feature[start_coord:end_coord] 
        
            
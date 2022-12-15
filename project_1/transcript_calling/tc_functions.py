import pandas as pd
import pysam
import bisect
from jw_utils import parse_gff as pgf
from datetime import date

def make_coverage_df(path_to_transcript_annots):
    """return coverage df
    
    parameters:
    """
    transcrpt_obj_dict = pgf.make_seq_object_dict(path_to_transcript_annots, feature_type='transcript')
    coverage = []
    name = []
    for trnscrpt_obj in transcrpt_obj_dict.values():
        coverage.append(float(trnscrpt_obj.coverage))
        name.append(trnscrpt_obj.ID)
    df = pd.DataFrame(coverage, index=name)
    df.columns = ['coverage']  
    return df


def get_sorted_df_dict(path_to_gff):
    contig_names = pgf.get_contig_names(path_to_gff)
    df = pgf.make_simple_annot_df(path_to_gff,start_end=True, contig=True)
    pl_gbo = df.loc[df.loc[:,'strand'] =='+',:].groupby('contig_name')
    mi_gbo = df.loc[df.loc[:,'strand'] =='-',:].groupby('contig_name')
    minus_df_dict = {contig:mi_gbo.get_group(contig).sort_values('start') for contig in contig_names}
    plus_df_dict = {contig:pl_gbo.get_group(contig).sort_values('start') for contig in contig_names}
    return minus_df_dict, plus_df_dict


        
def is_coordinate_within_gene(strand, coordinate, contig, strand_df_dict):
    "determine if the given gene coordinate on a contig is within a gene"
    if strand == '+':
        index = bisect.bisect_right(strand_df_dict[contig]['start'], coordinate)-1
        if index<0:
            return False
        elif coordinate <= strand_df_dict[contig].iloc[index, 4]:
            return True
        else: return False
    elif strand == '-':
        index = bisect.bisect_right(strand_df_dict[contig]['start'], coordinate)-1
        if index<0:
            return False
        elif coordinate <= strand_df_dict[contig].iloc[index, 4]:
            return True
        else: return False



# This calls the function directly above to get all transcription start sites that
# that do not begin within a gene on the same strand
def transcrpts_begin_intergene(path_to_strintie_GTF, path_to_gff, 
                            return_type='str', annot_type='gff',
                           return_within_gene =False):
    """return stringite transcripts IDs that do NOT start within a gene ON THE SAME STRAND
    
    parameters:
    return_type (str): type of item in dictionary, can be 1) 'str' returns seq obj ID or 2) 'obj' returns entire sequence object
    """
    #df_dict_minus = make_feature_startStop_df('-', path_to_gff)
    #df_dict_plus = make_feature_startStop_df('+', path_to_gff)
    if annot_type == 'gff':
        stringtie_seqobj_dict = pgf.make_seq_object_dict(path_to_strintie_GTF, feature_type='transcript')
    elif annot_type == 'gtf':  
        stringtie_seqobj_dict = make_GTF_seq_obj(path_to_strintie_GTF, feature_type='transcript')
    
    minus_df_dict, plus_df_dict = get_sorted_df_dict(path_to_gff)
    intergene_start_trnspts = {'+':[], '-':[]}
    withingene_start_trnspts = {'+':[], '-':[]}
    for trnscrpt_obj in stringtie_seqobj_dict.values():
        if trnscrpt_obj.strand=='-':
            starts_in_gene = is_coordinate_within_gene('-', trnscrpt_obj.end, trnscrpt_obj.chromosome, minus_df_dict)
            if not starts_in_gene:
                if return_type=='obj':
                    intergene_start_trnspts['-'].append(trnscrpt_obj)
                else:
                    intergene_start_trnspts['-'].append(trnscrpt_obj.ID)
            else:
                if return_type=='obj':
                    withingene_start_trnspts['-'].append(trnscrpt_obj)
                else:
                    withingene_start_trnspts['-'].append(trnscrpt_obj.ID)

        elif trnscrpt_obj.strand=='+':
            starts_in_gene = is_coordinate_within_gene('+', trnscrpt_obj.start, trnscrpt_obj.chromosome, plus_df_dict)
            if not starts_in_gene:
                if return_type=='obj':
                    intergene_start_trnspts['+'].append(trnscrpt_obj)
                else:
                    intergene_start_trnspts['+'].append(trnscrpt_obj.ID)
            else:
                if return_type=='obj':
                    withingene_start_trnspts['+'].append(trnscrpt_obj)
                else:
                    withingene_start_trnspts['+'].append(trnscrpt_obj.ID)

    if return_within_gene:
        return withingene_start_trnspts
    else:         
        return intergene_start_trnspts




class seq_attributes:
    'make a sequence object from a partially processed line in a stringtie GTF file'
    def __init__(self, line, anot_lst, attributes_dict, *args, **kwargs):
        self.line = line
        self.ID = attributes_dict['gene_id']
        self.transcript_id = attributes_dict.get('transcript_id')
        self.coverage = attributes_dict.get('cov')
        self.FPKM = attributes_dict.get('FPKM')
        self.TPM = attributes_dict.get('TPM')

        self.chromosome = anot_lst[0]
        self.source = anot_lst[1]
        self.feature_type = anot_lst[2]
        self.start = int(anot_lst[3])
        self.end = int(anot_lst[4])
        self.score = anot_lst[5]
        self.strand = anot_lst[6]   
        self.phase = anot_lst[7]



def make_GTF_seq_obj(path_to_GTF, feature_type='transcript'):
    "make a dict of seq objects from a stringtie GTF"
    
    with open(path_to_GTF, 'r') as f:
        seq_dict = {} 
        contig_dict = {}
        for i,line in enumerate(f):
            line = line.replace('"','')
            line = line.replace('\n','')[:-1]
            if line[0] != '#' and line[0] != ' ':
                anot_lst = line.split('\t')
                if len(anot_lst) <6:
                    raise Exception(f'line {i+1} in {path_to_GTF} does not contain all annotation fields')
                if anot_lst[2] == feature_type:
                    attributes = anot_lst[-1].replace('; ', ';')
                    if attributes.find('gene_id') !=-1:
                        temp = [key_val.split(' ') for key_val in attributes.split(';')]
                        attributes_dict = {key_val[0].strip():key_val[1].strip() for key_val in temp}
                        ID = attributes_dict['gene_id']
                        seq_dict[ID] = seq_attributes(line, anot_lst, attributes_dict)
    return seq_dict



def _make_wig_from_pilup(sam_align_obj, wig_path, contig):
    "make pilup wig file from bam alignment, genome position is 1-based"
    df = pd.DataFrame()
    pilup = sam_align_obj.pileup(contig=contig)
    file_lines = []
    header_line = f'variableStep chrom={contig}'
    file_lines.append(header_line)
    for col in pilup:
        line = f'{col.reference_pos+1} {col.nsegments}'
        file_lines.append(line)
    with open(wig_path, 'w') as f:
        for line in file_lines:
            f.write(line + '\n')
            
            
def make_wig_pilups(bam_filepath, contig_names):
    sam_align_obj = pysam.AlignmentFile(bam_filepath, 'rb')
    filename = bam_filepath.split('/')[-1].replace('.bam', '')
    for contig in contig_names:
        wig_path = f'./{contig}_pilup_{filename}.wig'
        _make_wig_from_pilup(sam_align_obj, wig_path, contig=contig)



def extend_to_stopCodon(path_to_strintie_GTF, path_to_gff):
    """
    For a each transcript object, if stop of the object exist within gene on same strand, then
    extend the stop of that object to the stop codon of that gene
    
    """
    annot_df = pgf.make_simple_annot_df(path_to_gff, start_end=True, contig=True).sort_values('end')
    plus_annot_grp = annot_df.loc[annot_df['strand']=='+',:].groupby('contig_name')
    minus_annot_grp = annot_df.loc[annot_df['strand']=='-',:].groupby('contig_name')
    contig_names = pgf.get_contig_names(path_to_gff)
    #determin if trn obj is within gene
    contigs = {contig:{} for contig in contig_names}
    df_contig_dict_minus = {contig:minus_annot_grp.get_group(contig) for contig in contig_names}
    df_contig_dict_plus = {contig:plus_annot_grp.get_group(contig) for contig in contig_names}
    trnascript_obj_dict = pgf.make_seq_object_dict(path_to_strintie_GTF, feature_type='transcript')
    for trans_obj in trnascript_obj_dict.values():
        if trans_obj.strand == '-':
            eot =  trans_obj.start  #eot = end of transcript
            contig=trans_obj.chromosome
            #if not is_coordinate_within_gene('-',eot,contig,)
            df = df_contig_dict_minus[trans_obj.chromosome]
            indr = bisect.bisect_right(list(df.loc[:,'start']), eot)-1 
            if indr>=0: 
                nearest_gene = df.iloc[indr,:].name # index? of series is .name
                if eot < df.loc[nearest_gene, 'end']: #if trnscrpt is within gene on same strand
                    new_trnscript_3pr_end = df.loc[nearest_gene, 'start']
                    trans_obj.start = new_trnscript_3pr_end
            contigs[trans_obj.chromosome][trans_obj.ID] = trans_obj
        if trans_obj.strand == '+':
            eot =  trans_obj.end  #eot = end of transcript
            df = df_contig_dict_plus[trans_obj.chromosome]
            indr = bisect.bisect_right(list(df.loc[:,'start']), eot)-1
            if indr>=0: 
                nearest_gene = df.iloc[indr,:].name # index? of series is .name
                if eot < df.loc[nearest_gene, 'end']:
                    new_trnscript_3pr_end = df.loc[nearest_gene, 'end']
                    trans_obj.end = new_trnscript_3pr_end
            contigs[trans_obj.chromosome][trans_obj.ID] = trans_obj
    return  contigs


def filter_transcripts(path_to_stringtie_gtf, path_to_gff, 
                       coverage_thresh,annot_type):
    """return nested dict {contig:{id:obj}} of seq objects that passed filters
    
    filters: coverage filter and intergene filter"""
    #keep transcripts that passed earlier filters: density and i
    new_trnscript_coords = extend_to_stopCodon(path_to_stringtie_gtf, path_to_gff)
    intergene_transcripts = transcrpts_begin_intergene(path_to_stringtie_gtf, 
                                                       path_to_gff, annot_type=annot_type)
    coverage_df = make_coverage_df(path_to_stringtie_gtf)
    filt = coverage_df.loc[:,'coverage']> coverage_thresh
    coverage_df_filt = coverage_df.loc[filt,:]
    new_trnscript_coords_filt = {contig:{} for contig in new_trnscript_coords}
    inter = intergene_transcripts['+'] + intergene_transcripts['-']
    for contig in new_trnscript_coords:
        for trnscrpt_id in new_trnscript_coords[contig]:
            if trnscrpt_id in inter:
                if trnscrpt_id in coverage_df_filt.index:
                    new_trnscript_coords_filt[contig][trnscrpt_id] = new_trnscript_coords[contig][trnscrpt_id]

    return new_trnscript_coords_filt


def merge_pos_neg_sttie_objs(path_pos_gtf, path_neg_gtf, path_to_gff):
    "renames and merges sequence annotation objects into single nested dict"

    stringtie_seqobj_dict_negs = make_GTF_seq_obj(path_neg_gtf, feature_type='transcript')
    stringtie_seqobj_dict_pos = make_GTF_seq_obj(path_pos_gtf, feature_type='transcript')
    merged_nested_obj_dict = {contig:{} for contig in pgf.get_contig_names(path_to_gff)}
    i=0
    for obj in stringtie_seqobj_dict_pos.values():
        i+=1
        obj.ID = f'STRG.{i}'
        obj.transcript_id = f'STRG.{i}'
        merged_nested_obj_dict[obj.chromosome][obj.ID]=obj
    for obj in stringtie_seqobj_dict_negs.values():
        i+=1
        obj.ID = f'STRG.{i}'
        obj.transcript_id = f'STRG.{i}'
        merged_nested_obj_dict[obj.chromosome][obj.ID]=obj
        
    return merged_nested_obj_dict


def write_gff_from_seqobj(trans_obj_dict, new_gtf_filepath):
    """Write a new gff file from dict seq objects.
    
    trans_obj_dict (dict): {contig:{id:seq object}}
    """
    with open(new_gtf_filepath, 'w') as f:
        f.write(f'#stringtie output filtered and edited with filter_stringtie_results.ipynb on {date.today()}\n')
    for contig in trans_obj_dict.keys():
        lines = {}
        for transcpt_obj in trans_obj_dict[contig].values():
            ft=transcpt_obj.feature_type
            s = transcpt_obj.start
            e = transcpt_obj.end
            sc = '.'
            st = transcpt_obj.strand
            p = transcpt_obj.phase
            ID = transcpt_obj.ID
            TI = transcpt_obj.transcript_id
            cv = transcpt_obj.coverage
            FPKM = transcpt_obj.FPKM
            TPM = transcpt_obj.TPM
            ln=f'{contig}\tTrestleBio\t{ft}\t{s}\t{e}\t{sc}\t{st}\t{p}\t'
            attr = f'ID={ID};transcript_id={TI};cov={cv};FPKM={FPKM};TPM={TPM}'
            lines[s] = ln+attr+'\n'
        df = pd.DataFrame(lines.keys(), lines.values())
        df = df.sort_values(df.columns[0])
        with open(new_gtf_filepath, 'a') as f:
            for line in df.index:
                f.write(line)


from __future__ import division

import csv
import sys
import copy
import operator
import numpy
import math
import scipy
import pandas

termini_percentage_const = 10
th_asa_1 = 0.1
th_asa_2 = 0.07
th_iupred_bin = 0.5
th_profbval_bin_1 = -0.1
th_profbval_bin_2 = 0.5

cached_ID = ''
cached_AA_fractions = None
cached_AA_indices = None
cached_AA_ranks = None



# def util_get_amino_acid_indices():
#     df= pandas.read_csv("amino_acid_data.csv",index_col= "ID")
#     data = df.transpose().to_dict()
#     return data

def util_get_amino_acid_indices():
  return  {'Flexibility': {'A': 0.35700000000000004,
  'R': 0.529,
  'N': 0.46299999999999997,
  'D': 0.511,
  'C': 0.34600000000000003,
  'Q': 0.493,
  'E': 0.49700000000000005,
  'G': 0.544,
  'H': 0.32299999999999995,
  'I': 0.462,
  'L': 0.365,
  'K': 0.466,
  'M': 0.295,
  'F': 0.314,
  'P': 0.509,
  'S': 0.507,
  'T': 0.444,
  'W': 0.305,
  'Y': 0.42,
  'V': 0.386},
 'B-Value': {'A': 0.9840000000000001,
  'R': 1.008,
  'N': 1.048,
  'D': 1.068,
  'C': 0.9059999999999999,
  'Q': 1.037,
  'E': 1.094,
  'G': 1.031,
  'H': 0.95,
  'I': 0.927,
  'L': 0.935,
  'K': 1.102,
  'M': 0.9520000000000001,
  'F': 0.915,
  'P': 1.0490000000000002,
  'S': 1.046,
  'T': 0.997,
  'W': 0.904,
  'Y': 0.929,
  'V': 0.9309999999999999},
 'Top-IDP': {'A': 0.06,
  'R': 0.18,
  'N': 0.006999999999999999,
  'D': 0.192,
  'C': 0.02,
  'Q': 0.318,
  'E': 0.736,
  'G': 0.166,
  'H': 0.303,
  'I': -0.486,
  'L': -0.326,
  'K': 0.586,
  'M': -0.397,
  'F': -0.6970000000000001,
  'P': 0.987,
  'S': 0.341,
  'T': 0.059000000000000004,
  'W': -0.884,
  'Y': -0.51,
  'V': -0.121},
 'FoldUnfold': {'A': 19.89,
  'R': 21.03,
  'N': 18.49,
  'D': 17.41,
  'C': 23.52,
  'Q': 19.23,
  'E': 17.46,
  'G': 17.11,
  'H': 21.72,
  'I': 25.71,
  'L': 25.36,
  'K': 18.19,
  'M': 24.82,
  'F': 27.18,
  'P': 17.43,
  'S': 17.67,
  'T': 19.81,
  'W': 28.48,
  'Y': 25.93,
  'V': 23.93},
 'DisProt': {'A': 0.042,
  'R': 0.21100000000000002,
  'N': -0.106,
  'D': 0.127,
  'C': -0.546,
  'Q': 0.381,
  'E': 0.469,
  'G': 0.095,
  'H': -0.127,
  'I': -393.0,
  'L': -0.26,
  'K': 0.37,
  'M': 0.19699999999999998,
  'F': -0.381,
  'P': 0.419,
  'S': 0.201,
  'T': -0.11599999999999999,
  'W': -0.465,
  'Y': -0.42700000000000005,
  'V': -0.302}}

def util_get_amino_acid_ranks():
    global cached_AA_indices
    if cached_AA_indices == None:
        AA_data = util_get_amino_acid_indices()
        cached_AA_indices = AA_data
    else:
        AA_data = cached_AA_indices
    AA_ranks = {}
    for idx in AA_data:
        AAs  = numpy.array(list(AA_data[idx].keys()))
        AA_values  = numpy.array(list(AA_data[idx].values()))
        sorted_postions  = numpy.argsort(AA_values)
        AA_ranks[idx] = AAs[sorted_postions]
    return AA_ranks






def util_mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)



def util_find_AAs_with_specified_ranks_for_an_AA_index(index_name,rank_range):
    global cached_AA_ranks
    rank_start, rank_end = rank_range
    if cached_AA_ranks==None:
        ranks = util_get_amino_acid_ranks()
        cached_AA_ranks = ranks
    else:
        ranks = cached_AA_ranks
    sorted_byidx=ranks[index_name]
    AA_set = sorted_byidx[rank_start:rank_end]
    return AA_set




def template_sliding_window(chains, chain_level_prediction_name, op_chains, op_sw, sw_size, termini_percent, th_asa):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl, sw_size)
    empty_window_threshold = 3
    non_empty_window_exist = False
    #the window size used should be no less than the size of the smallest chain
    ch_values = []
    if op_sw == max:
        default_window_value = 0
    elif op_sw == min:
        default_window_value = 10E6
    for ch in chains:
        if len(ch['sequence']) * 80/100 <= window_size:
            return None
        #ch value is the value to be calculated for current chain which then will be aggregated between them
        ch_value = 0
        temp = 0
        num = 1
        operation_range_beginning = int(math.floor((termini_percent/100)*len(ch['sequence'])))
        operation_range_end = len(ch['sequence'])-operation_range_beginning
        aa = ch[chain_level_prediction_name]
        asaquick = ch["asaquick"]
        #calculating the value in the first window from the begining of the working range and counting the number of the elements that are used (on surface) based on their asaquick
        for i in range(operation_range_beginning, operation_range_beginning+window_size):
            if i>= len(asaquick):
                a = 7
            if asaquick[i] > th_asa:
                temp += aa[i]
                num+=1
        #if the window is not empty then we use the obtained value
        if num >= empty_window_threshold:
            non_empty_window_exist =True
            ch_value = temp/num

        else: #otherwise we use default value
            ch_value = default_window_value

        for i in range(operation_range_beginning, operation_range_end - window_size):
            if asaquick[i] > th_asa:
                temp -= aa[i]
                num -=1
            if asaquick[i + window_size] > th_asa:
                temp += aa[i + window_size]
                num+= 1
            if num >= empty_window_threshold:
                non_empty_window_exist = True
                ch_value = op_sw(temp/num, ch_value)
            # else:
            # do nothing with min and max values (ignore this window)
        ch_values.append(ch_value)
    as_value = op_chains(ch_values)
    if not non_empty_window_exist:
        return None
    else:
        return as_value


def template_max_min_sliding_window(chains, chain_level_prediction_name, op_chains, sw_size, termini_percent, th_asa):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl, sw_size)
    ch_values = []
    empty_window_threshold = 3
    #the window size used should be no less than the size of the smallest chain
    default_min = 10E6
    default_max = 0
    non_empty_window_exist = False
    for ch in chains:
        if len(ch['sequence']) * 80/100 <= window_size:
            return None
        ch_value = 0
        temp = 0
        num = 0
        ch_max = default_max
        ch_min = default_min
        aa = ch[chain_level_prediction_name]
        asaquick = ch["asaquick"]
        operation_range_beginning = int(math.floor((termini_percent/100)*len(ch['sequence'])))
        operation_range_end = len(ch['sequence'])-operation_range_beginning
        for i in range(operation_range_beginning, operation_range_beginning+window_size):
            if asaquick[i] > th_asa:
                temp += aa[i]
                num+=1
        if num >empty_window_threshold:
            non_empty_window_exist = True
            ch_max = max(temp/num,ch_max)
            ch_min = min(temp/num,ch_min)
        else:
            ch_max = default_max
            ch_min = default_min
        for i in range(operation_range_beginning, operation_range_end - window_size):
            if asaquick[i] > th_asa:
                temp -= aa[i]
                num -= 1
            if asaquick[i + window_size] > th_asa:
                temp += aa[i + window_size]
                num += 1
            if num > empty_window_threshold:
                non_empty_window_exist = True
                ch_max = max(temp/num, ch_max)
                ch_min = min(temp/num, ch_min)
            # else:
                #do nothing with min and max values (ignore this window)
        ch_values.append(ch_max-ch_min)
    as_value = op_chains(ch_values)
    if not non_empty_window_exist:
        return None
    else:
        return as_value



def util_amino_acid_fractions(chains):
    global cached_AA_fractions
    global cached_ID
    ID = chains[0]["id"]
    if ID == cached_ID:
        return cached_AA_fractions
    cached_AA_fractions = {'A' : [] ,'R' : [] ,'N' : [],'D' : [] ,'C' : [],'Q' : [],'E' : [],'G' : [],'H': [],'I':[],'L':[],'K':[],'M':[],'F':[],'P':[],'S':[],'T':[],'W':[],'Y':[],'V':[]}
    i = 0
    values  = {}
    for ch in chains:
        seq = ch["sequence"]
        chain_length = len(seq)
        num = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 'M': 0,'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
        #
        # for k in num.keys():
        #     num[k].append(0)
        for residue in seq:
            if residue in num.keys():
                num[residue] += 1
            else:
                print( "Observed an out of range amino acid in " , ch["PDB ID"] ,"  ",ch["Chain ID"], "  : ",residue)
        i+=1
        for k in cached_AA_fractions.keys():
            cached_AA_fractions[k].append(float(num[k])/chain_length)
    cached_ID = ID
    return cached_AA_fractions





def template_amino_acid_set_fraction(chains, AAset, op_chain):
    frac = util_amino_acid_fractions(chains)
    sum = 0
    for t in AAset:
        sum += op_chain(frac[t])
        # amino_acid_fractions_max(chains)
    return sum


def frac_maxc_sidechainclass_aliphatic(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'I', 'L', 'V'], max)

def frac_maxc_sidechainclass_basic(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'H', 'K'], max)


def frac_maxc_sidechainclass_amide(chains):
    return template_amino_acid_set_fraction(chains, ['N', 'Q'], max)


def frac_maxc_sidechainclass_acid(chains):
    return template_amino_acid_set_fraction(chains, ['E', 'D'], max)

def frac_maxc_sidechainclass_sulfurcontaining(chains):
    return template_amino_acid_set_fraction(chains, ['C', 'M'], max)

def frac_maxc_sidechainclass_aromatic(chains):
    return template_amino_acid_set_fraction(chains, ['H', 'F', 'W', 'Y'], max)


def frac_maxc_sidechainpolarity_nonpolar(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'C', 'G', 'I', 'M', 'F', 'P', 'W', 'V'], max)

def frac_maxc_sidechainpolarity_basicpolar(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'H', 'K'], max)


def frac_maxc_sidechainpolarity_polar(chains):
    return template_amino_acid_set_fraction(chains, ['N', 'Q', 'S', 'T', 'Y'], max)

def frac_maxc_sidechainpolarity_acidicpolar(chains):
    return template_amino_acid_set_fraction(chains, ['D', 'E'], max)

def frac_maxc_sidechaincharge_neutral(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'H'], max)



def frac_maxc_sidechaincharge_positive(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'K'], max)


def frac_maxc_sidechaincharge_negative(chains):
    return template_amino_acid_set_fraction(chains, ['D', 'E'], max)

def frac_maxc_heavyHFRYW(chains):
    return template_amino_acid_set_fraction(chains, ['H', 'F', 'R', 'Y', 'W'], max)
def frac_maxc_lightGASPV(chains):
    return template_amino_acid_set_fraction(chains, ['G', 'A', 'S', 'P', 'V'], max)

def frac_maxc_lightGASPV(chains):
    return template_amino_acid_set_fraction(chains, ['G', 'A', 'S', 'P', 'V'], max)

def frac_maxc_hydrophob(chains):
    return template_amino_acid_set_fraction(chains, ['F', 'Y', 'W', 'M', 'I', 'L', 'V', 'C', 'G', 'A', 'C', 'T', 'H', 'K'], max)

def frac_maxc_tinyAGCS(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'C', 'S'], max)

def frac_maxc_smallAGCSVTDNP(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'C', 'S', 'V', 'T', 'D', 'N', 'P'], max)







###############################################



def frac_minc_sidechainclass_aliphatic(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'I', 'L', 'V'], min)

def frac_minc_sidechainclass_basic(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'H', 'K'], min)


def frac_minc_sidechainclass_amide(chains):
    return template_amino_acid_set_fraction(chains, ['N', 'Q'], min)


def frac_minc_sidechainclass_acid(chains):
    return template_amino_acid_set_fraction(chains, ['E', 'D'], min)

def frac_minc_sidechainclass_sulfurcontaining(chains):
    return template_amino_acid_set_fraction(chains, ['C', 'M'], min)

def frac_minc_sidechainclass_aromatic(chains):
    return template_amino_acid_set_fraction(chains, ['H', 'F', 'W', 'Y'], min)


def frac_minc_sidechainpolarity_nonpolar(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'C', 'G', 'I', 'M', 'F', 'P', 'W', 'V'], min)

def frac_minc_sidechainpolarity_basicpolar(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'H', 'K'], min)


def frac_minc_sidechainpolarity_polar(chains):
    return template_amino_acid_set_fraction(chains, ['N', 'Q', 'S', 'T', 'Y'], min)

def frac_minc_sidechainpolarity_acidicpolar(chains):
    return template_amino_acid_set_fraction(chains, ['D', 'E'], min)

def frac_minc_sidechaincharge_neutral(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'H'], min)



def frac_minc_sidechaincharge_positive(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'K'], min)


def frac_minc_sidechaincharge_negative(chains):
    return template_amino_acid_set_fraction(chains, ['D', 'E'], min)

def frac_minc_heavyHFRYW(chains):
    return template_amino_acid_set_fraction(chains, ['H', 'F', 'R', 'Y', 'W'], min)
def frac_minc_lightGASPV(chains):
    return template_amino_acid_set_fraction(chains, ['G', 'A', 'S', 'P', 'V'], min)

def frac_minc_lightGASPV(chains):
    return template_amino_acid_set_fraction(chains, ['G', 'A', 'S', 'P', 'V'], min)

def frac_minc_hydrophob(chains):
    return template_amino_acid_set_fraction(chains, ['F', 'Y', 'W', 'M', 'I', 'L', 'V', 'C', 'G', 'A', 'C', 'T', 'H', 'K'], min)

def frac_minc_tinyAGCS(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'C', 'S'], min)

def frac_minc_smallAGCSVTDNP(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'C', 'S', 'V', 'T', 'D', 'N', 'P'], min)



####################################################################







def frac_avgc_sidechainclass_aliphatic(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'I', 'L', 'V'], util_mean)

def frac_avgc_sidechainclass_basic(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'H', 'K'], util_mean)


def frac_avgc_sidechainclass_amide(chains):
    return template_amino_acid_set_fraction(chains, ['N', 'Q'], util_mean)


def frac_avgc_sidechainclass_acid(chains):
    return template_amino_acid_set_fraction(chains, ['E', 'D'], util_mean)

def frac_avgc_sidechainclass_sulfurcontaining(chains):
    return template_amino_acid_set_fraction(chains, ['C', 'M'], util_mean)

def frac_avgc_sidechainclass_aromatic(chains):
    return template_amino_acid_set_fraction(chains, ['H', 'F', 'W', 'Y'], util_mean)


def frac_avgc_sidechainpolarity_nonpolar(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'C', 'G', 'I', 'M', 'F', 'P', 'W', 'V'], util_mean)

def frac_avgc_sidechainpolarity_basicpolar(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'H', 'K'], util_mean)


def frac_avgc_sidechainpolarity_polar(chains):
    return template_amino_acid_set_fraction(chains, ['N', 'Q', 'S', 'T', 'Y'], util_mean)

def frac_avgc_sidechainpolarity_acidicpolar(chains):
    return template_amino_acid_set_fraction(chains, ['D', 'E'], util_mean)

def frac_avgc_sidechaincharge_neutral(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'H'], util_mean)



def frac_avgc_sidechaincharge_positive(chains):
    return template_amino_acid_set_fraction(chains, ['R', 'K'], util_mean)


def frac_avgc_sidechaincharge_negative(chains):
    return template_amino_acid_set_fraction(chains, ['D', 'E'], util_mean)

def frac_avgc_heavyHFRYW(chains):
    return template_amino_acid_set_fraction(chains, ['H', 'F', 'R', 'Y', 'W'], util_mean)
def frac_avgc_lightGASPV(chains):
    return template_amino_acid_set_fraction(chains, ['G', 'A', 'S', 'P', 'V'], util_mean)

def frac_avgc_lightGASPV(chains):
    return template_amino_acid_set_fraction(chains, ['G', 'A', 'S', 'P', 'V'], util_mean)

def frac_avgc_hydrophob(chains):
    return template_amino_acid_set_fraction(chains, ['F', 'Y', 'W', 'M', 'I', 'L', 'V', 'C', 'G', 'A', 'C', 'T', 'H', 'K'], util_mean)

def frac_avgc_tinyAGCS(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'C', 'S'], util_mean)

def frac_avgc_smallAGCSVTDNP(chains):
    return template_amino_acid_set_fraction(chains, ['A', 'G', 'C', 'S', 'V', 'T', 'D', 'N', 'P'], util_mean)
















def template_max_sliding_window(chains, chain_level_prediction_name, op_chains, sw_size, termini_percent, th_asa):
    return template_sliding_window(chains, chain_level_prediction_name, op_chains, max, sw_size, termini_percent, th_asa)




def template_max_consecutive_1s_in_bin(chains, chain_level_prediction_name, op_chains, binerization_threshold, termini_percent):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        operation_range_beginning = int(math.floor((termini_percent/100)*len(ch['sequence'])))
        operation_range_end = len(ch['sequence'])-operation_range_beginning
        aa_values_f = ch[chain_level_prediction_name]
        aa_values = []
        i = 0
        for f in aa_values_f:
            if f>binerization_threshold:
                aa_values.append(1)
            else:
                aa_values.append(0)
        count = 0
        max_consecutives = 0
        for i in range(operation_range_beginning,operation_range_end):
            if aa_values[i] ==1:
                count+=1
                max_consecutives = max(count, max_consecutives)
            else:
                count = 0
        ch_values.append(max_consecutives)
    as_value = op_chains(ch_values)
    return as_value


def template_max_consecutive_1s_in_bin_per_length(chains, chain_level_prediction_name, op_chains, binerization_threshold, termini_percent):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        operation_range_beginning = int(math.floor((termini_percent/100)*len(ch['sequence'])))
        operation_range_end = len(ch['sequence'])-operation_range_beginning
        aa_values_f = ch[chain_level_prediction_name]
        aa_values  = []
        i = 0
        for f in aa_values_f:
            if f>binerization_threshold:
                aa_values.append(1)
            else:
                aa_values.append(0)
        count = 0
        max_consecutives = 0

        for i in range(operation_range_beginning,operation_range_end):
            if aa_values[i] ==1:
                count+=1
                max_consecutives = max(count, max_consecutives)
            else:
                count = 0
        ch_values.append(float(max_consecutives)/(operation_range_end-operation_range_beginning))
    as_value = op_chains(ch_values)
    return as_value




def template_frac_of_1s_in_binary(chains, chain_level_prediction_name, op_chains, binerization_threshold, termini_percent, th_asa):
    as_max = 0
    ch_values = []
    for ch in chains:
        ch_max = 0
        temp = 0
        operation_range_beginning = int(math.floor((termini_percent/100)*len(ch['sequence'])))
        operation_range_end = len(ch['sequence'])-operation_range_beginning
        aa_values = ch[chain_level_prediction_name]
        aa_values_bin =[]
        i = 0
        asaquick = ch["asaquick"]
        for i in range(operation_range_beginning,operation_range_end):
            f = aa_values[i]
            if asaquick[i]>th_asa:
                if f>binerization_threshold:
                    aa_values_bin.append(1)
                else:
                    aa_values_bin.append(0)
        frac = float(sum(aa_values_bin))/len(aa_values_bin)
        ch_values.append(frac)
    as_value= op_chains(ch_values)
    return as_value








def template_avg(chains, chain_level_prediction_name, op_chains, termini_percent, th_asa):
    ch_values =[]
    for ch in chains:
        operation_range_beginning = int(math.floor((termini_percent/100)*len(ch['sequence'])))
        operation_range_end = len(ch['sequence'])-operation_range_beginning
        ch_max = 0
        temp = 0
        aa= ch[chain_level_prediction_name]
        filtered_aa = []
        asaquick = ch['asaquick']

        for i in range(operation_range_beginning,operation_range_end):
            if asaquick[i]>th_asa:
                filtered_aa.append(aa[i])
        ch_values.append(util_mean(filtered_aa))
    as_max = op_chains(ch_values)
    return as_max

























def no_unique_chains(chains):
    seqs = []
    for ch in chains:
        seqs.append(ch["sequence"])
    res = len(set(seqs))
    return res

def sum_chain_length(chains):
    sum = 0
    for ch in chains:
        sum += len(ch["sequence"])
    return sum




def average_chain_length(chains):
    sum = 0
    num = 0
    for ch in chains:
        sum += len(ch["sequence"])
        num +=1
    avg = sum / num
    return avg

def shortest_chain_length(chains):
    min = sys.maxsize
    for ch in chains:
        if (len(ch["sequence"])<min):
            min = len(ch["sequence"])
    return min

def longest_chain_length(chains):
    max = 0
    for ch in chains:
        if (len(ch["sequence"])>max):
            max = len(ch["sequence"])
    return max

def num_chains(chains):
    return len(chains)

def stdev_chains(chains):
    ls = []
    for ch in chains:
        ls.append(len(ch["sequence"]))
    sd = float(numpy.std(numpy.array(ls)))
    return sd






def max_avg_in_top10_percent_profbval(chains):
    perc = 30
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval = ch['profbval']
        i = 0
        sorted_profbval = sorted(profbval)
        boundary = numpy.floor((1-perc/100)*len(profbval))
        sum =0
        count =0
        for i in range(int(boundary),len(profbval)):
            sum += sorted_profbval[i]
            count += 1
        avg = sum/count
    return avg




def min_avg_in_bottom10_precent_profbval(chains):
    perc = 30
    for ch in chains:
        ch_max = 0
        temp = 0
        profbval = ch['profbval']
        sorted_profbval = sorted(profbval)
        boundary = numpy.floor((perc/100)*len(profbval))
        sum =0
        count =0
        for i in range(0,int(boundary)):
            sum+=sorted_profbval[i]
            count+=1
        avg = sum/count
    return avg





#######################################################




def template_fraction_of_bottom3AA(chains,index_name,op_chain):
    AA_set = util_find_AAs_with_specified_ranks_for_an_AA_index(index_name,(0,3))
    return template_amino_acid_set_fraction(chains, AA_set, op_chain)



def template_fraction_of_top3AA(chains,index_name,op_chain):
    AA_set = util_find_AAs_with_specified_ranks_for_an_AA_index(index_name,(17,20))
    return template_amino_acid_set_fraction(chains, AA_set, op_chain)






def frac_maxc_top3Flexbilityidx(chains):
    return template_fraction_of_top3AA(chains,'Flexibility',max)


def frac_maxc_top3bvalueidx(chains):
    return template_fraction_of_top3AA(chains,'B-Value',max)


def frac_maxc_top3topidp(chains):
    return template_fraction_of_top3AA(chains, 'Top-IDP', max)

def frac_maxc_top3foldunfold(chains):
    return template_fraction_of_top3AA(chains, 'FoldUnfold', max)


def frac_maxc_top3disprot(chains):
    return template_fraction_of_top3AA(chains, 'DisProt', max)


####################################################


def frac_maxc_bottom3Flexbilityidx(chains):
    return template_fraction_of_bottom3AA(chains,'Flexibility',max)

def frac_maxc_bottom3bvalueidx(chains):
    return template_fraction_of_bottom3AA(chains,'B-Value',max)

def frac_maxc_bottom3topidp(chains):
    return template_fraction_of_bottom3AA(chains, 'Top-IDP', max)

def frac_maxc_bottom3foldunfold(chains):
    return template_fraction_of_bottom3AA(chains, 'FoldUnfold', max)

def frac_maxc_bottom3disprot(chains):
    return template_fraction_of_bottom3AA(chains, 'DisProt', max)

#####################################################




def frac_minc_top3Flexbilityidx(chains):
    return template_fraction_of_top3AA(chains,'Flexibility',min)


def frac_minc_top3bvalueidx(chains):
    return template_fraction_of_top3AA(chains,'B-Value',min)


def frac_minc_top3topidp(chains):
    return template_fraction_of_top3AA(chains, 'Top-IDP', min)

def frac_minc_top3foldunfold(chains):
    return template_fraction_of_top3AA(chains, 'FoldUnfold', min)


def frac_minc_top3disprot(chains):
    return template_fraction_of_top3AA(chains, 'DisProt', min)


####################################################


def frac_minc_bottom3Flexbilityidx(chains):
    return template_fraction_of_bottom3AA(chains,'Flexibility',min)

def frac_minc_bottom3bvalueidx(chains):
    return template_fraction_of_bottom3AA(chains,'B-Value',min)

def frac_minc_bottom3topidp(chains):
    return template_fraction_of_bottom3AA(chains, 'Top-IDP', min)

def frac_minc_bottom3foldunfold(chains):
    return template_fraction_of_bottom3AA(chains, 'FoldUnfold', min)

def frac_minc_bottom3disprot(chains):
    return template_fraction_of_bottom3AA(chains, 'DisProt', min)


####################################################


def frac_avgc_top3Flexbilityidx(chains):
    return template_fraction_of_top3AA(chains,'Flexibility', util_mean)


def frac_avgc_top3bvalueidx(chains):
    return template_fraction_of_top3AA(chains,'B-Value', util_mean)


def frac_avgc_top3topidp(chains):
    return template_fraction_of_top3AA(chains, 'Top-IDP', util_mean)

def frac_avgc_top3foldunfold(chains):
    return template_fraction_of_top3AA(chains, 'FoldUnfold', util_mean)


def frac_avgc_top3disprot(chains):
    return template_fraction_of_top3AA(chains, 'DisProt', util_mean)


####################################################


def frac_avgc_bottom3Flexbilityidx(chains):
    return template_fraction_of_bottom3AA(chains,'Flexibility', util_mean)

def frac_avgc_bottom3bvalueidx(chains):
    return template_fraction_of_bottom3AA(chains,'B-Value', util_mean)

def frac_avgc_bottom3topidp(chains):
    return template_fraction_of_bottom3AA(chains, 'Top-IDP', util_mean)

def frac_avgc_bottom3foldunfold(chains):
    return template_fraction_of_bottom3AA(chains, 'FoldUnfold', util_mean)

def frac_avgc_bottom3disprot(chains):
    return template_fraction_of_bottom3AA(chains, 'DisProt', util_mean)


####################################################




def frac_maxc_AA_A(chains):
    return template_amino_acid_set_fraction(chains, ['A'], max)
def frac_maxc_AA_R(chains):
    return template_amino_acid_set_fraction(chains, ['R'], max)
def frac_maxc_AA_N(chains):
    return template_amino_acid_set_fraction(chains, ['N'], max)
def frac_maxc_AA_D(chains):
    return template_amino_acid_set_fraction(chains, ['D'], max)
def frac_maxc_AA_C(chains):
    return template_amino_acid_set_fraction(chains, ['C'], max)
def frac_maxc_AA_Q(chains):
    return template_amino_acid_set_fraction(chains, ['Q'], max)
def frac_maxc_AA_E(chains):
    return template_amino_acid_set_fraction(chains, ['E'], max)
def frac_maxc_AA_G(chains):
    return template_amino_acid_set_fraction(chains, ['G'], max)
def frac_maxc_AA_H(chains):
    return template_amino_acid_set_fraction(chains, ['H'], max)
def frac_maxc_AA_I(chains):
    return template_amino_acid_set_fraction(chains, ['I'], max)
def frac_maxc_AA_L(chains):
    return template_amino_acid_set_fraction(chains, ['L'], max)
def frac_maxc_AA_K(chains):
    return template_amino_acid_set_fraction(chains, ['K'], max)
def frac_maxc_AA_M(chains):
    return template_amino_acid_set_fraction(chains, ['M'], max)
def frac_maxc_AA_F(chains):
    return template_amino_acid_set_fraction(chains, ['F'], max)
def frac_maxc_AA_P(chains):
    return template_amino_acid_set_fraction(chains, ['P'], max)
def frac_maxc_AA_S(chains):
    return template_amino_acid_set_fraction(chains, ['S'], max)
def frac_maxc_AA_T(chains):
    return template_amino_acid_set_fraction(chains, ['T'], max)
def frac_maxc_AA_W(chains):
    return template_amino_acid_set_fraction(chains, ['W'], max)
def frac_maxc_AA_Y(chains):
    return template_amino_acid_set_fraction(chains, ['Y'], max)
def frac_maxc_AA_V(chains):
    return template_amino_acid_set_fraction(chains, ['V'], max)





def frac_minc_AA_A(chains):
    return template_amino_acid_set_fraction(chains, ['A'], min)
def frac_minc_AA_R(chains):
    return template_amino_acid_set_fraction(chains, ['R'], min)
def frac_minc_AA_N(chains):
    return template_amino_acid_set_fraction(chains, ['N'], min)
def frac_minc_AA_D(chains):
    return template_amino_acid_set_fraction(chains, ['D'], min)
def frac_minc_AA_C(chains):
    return template_amino_acid_set_fraction(chains, ['C'], min)
def frac_minc_AA_Q(chains):
    return template_amino_acid_set_fraction(chains, ['Q'], min)
def frac_minc_AA_E(chains):
    return template_amino_acid_set_fraction(chains, ['E'], min)
def frac_minc_AA_G(chains):
    return template_amino_acid_set_fraction(chains, ['G'], min)
def frac_minc_AA_H(chains):
    return template_amino_acid_set_fraction(chains, ['H'], min)
def frac_minc_AA_I(chains):
    return template_amino_acid_set_fraction(chains, ['I'], min)
def frac_minc_AA_L(chains):
    return template_amino_acid_set_fraction(chains, ['L'], min)
def frac_minc_AA_K(chains):
    return template_amino_acid_set_fraction(chains, ['K'], min)
def frac_minc_AA_M(chains):
    return template_amino_acid_set_fraction(chains, ['M'], min)
def frac_minc_AA_F(chains):
    return template_amino_acid_set_fraction(chains, ['F'], min)
def frac_minc_AA_P(chains):
    return template_amino_acid_set_fraction(chains, ['P'], min)
def frac_minc_AA_S(chains):
    return template_amino_acid_set_fraction(chains, ['S'], min)
def frac_minc_AA_T(chains):
    return template_amino_acid_set_fraction(chains, ['T'], min)
def frac_minc_AA_W(chains):
    return template_amino_acid_set_fraction(chains, ['W'], min)
def frac_minc_AA_Y(chains):
    return template_amino_acid_set_fraction(chains, ['Y'], min)
def frac_minc_AA_V(chains):
    return template_amino_acid_set_fraction(chains, ['V'], min)




def frac_avgc_AA_A(chains):
    return template_amino_acid_set_fraction(chains, ['A'], util_mean)
def frac_avgc_AA_R(chains):
    return template_amino_acid_set_fraction(chains, ['R'], util_mean)
def frac_avgc_AA_N(chains):
    return template_amino_acid_set_fraction(chains, ['N'], util_mean)
def frac_avgc_AA_D(chains):
    return template_amino_acid_set_fraction(chains, ['D'], util_mean)
def frac_avgc_AA_C(chains):
    return template_amino_acid_set_fraction(chains, ['C'], util_mean)
def frac_avgc_AA_Q(chains):
    return template_amino_acid_set_fraction(chains, ['Q'], util_mean)
def frac_avgc_AA_E(chains):
    return template_amino_acid_set_fraction(chains, ['E'], util_mean)
def frac_avgc_AA_G(chains):
    return template_amino_acid_set_fraction(chains, ['G'], util_mean)
def frac_avgc_AA_H(chains):
    return template_amino_acid_set_fraction(chains, ['H'], util_mean)
def frac_avgc_AA_I(chains):
    return template_amino_acid_set_fraction(chains, ['I'], util_mean)
def frac_avgc_AA_L(chains):
    return template_amino_acid_set_fraction(chains, ['L'], util_mean)
def frac_avgc_AA_K(chains):
    return template_amino_acid_set_fraction(chains, ['K'], util_mean)
def frac_avgc_AA_M(chains):
    return template_amino_acid_set_fraction(chains, ['M'], util_mean)
def frac_avgc_AA_F(chains):
    return template_amino_acid_set_fraction(chains, ['F'], util_mean)
def frac_avgc_AA_P(chains):
    return template_amino_acid_set_fraction(chains, ['P'], util_mean)
def frac_avgc_AA_S(chains):
    return template_amino_acid_set_fraction(chains, ['S'], util_mean)
def frac_avgc_AA_T(chains):
    return template_amino_acid_set_fraction(chains, ['T'], util_mean)
def frac_avgc_AA_W(chains):
    return template_amino_acid_set_fraction(chains, ['W'], util_mean)
def frac_avgc_AA_Y(chains):
    return template_amino_acid_set_fraction(chains, ['Y'], util_mean)
def frac_avgc_AA_V(chains):
    return template_amino_acid_set_fraction(chains, ['V'], util_mean)


################################################################################################

def template_max_number_of_occurance_in_sliding_window(chains,AA_set,window_size,op_chain):
    min_cl = shortest_chain_length(chains)
    window_size = min(min_cl,window_size)
    ch_values =[]
    for ch in chains:
        temp = 0
        seq = ch['sequence']
        for i in range(0, window_size):
            if seq[i] in AA_set:
                temp += 1
        ch_max=temp
        for i in range(0,len(seq)-window_size):
            if seq[i] in AA_set:
                temp -= 1
            if seq[i+window_size] in AA_set:
                temp += 1
            if temp > ch_max:
                ch_max = temp
        ch_values.append(ch_max)
    return op_chain(ch_values)




def idx_maxc_maxswnum20_top3_flex(chains):
    AA_set = util_find_AAs_with_specified_ranks_for_an_AA_index("Flexibility",(17,20))
    return template_max_number_of_occurance_in_sliding_window(chains, AA_set, 20, max)

def idx_maxc_maxswnum20_top3_bvalue(chains):
    AA_set = util_find_AAs_with_specified_ranks_for_an_AA_index('B-Value',(17,20))
    return template_max_number_of_occurance_in_sliding_window(chains, AA_set, 20, max)



def idx_maxc_maxswnum20_top3_topidp(chains):
    AA_set = util_find_AAs_with_specified_ranks_for_an_AA_index('Top-IDP',(17,20))
    return template_max_number_of_occurance_in_sliding_window(chains, AA_set, 20, max)


def idx_maxc_maxswnum20_top3_foldunfold(chains):
    AA_set = util_find_AAs_with_specified_ranks_for_an_AA_index('FoldUnfold',(17,20))
    return template_max_number_of_occurance_in_sliding_window(chains, AA_set, 20, max)



def idx_maxc_maxswnum20_top3_disprot(chains):
    AA_set = util_find_AAs_with_specified_ranks_for_an_AA_index('DisProt',(17,20))
    return template_max_number_of_occurance_in_sliding_window(chains, AA_set, 20, max)



###########################################################################################################





def	profbval_maxc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 10, termini_percentage_const, 0)
def	profbval_maxc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 20, termini_percentage_const, 0)
def	profbval_maxc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 25, termini_percentage_const, 0)
def	profbval_maxc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 30, termini_percentage_const, 0)
def	profbval_maxc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 35, termini_percentage_const, 0)
def	profbval_maxc_maxcons1sbinth1		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'profbval', max, th_profbval_bin_1, termini_percentage_const)
def	profbval_maxc_maxcons1sbinth1pl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'profbval', max, th_profbval_bin_1, termini_percentage_const)
def	profbval_maxc_maxcons1sbinth2		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'profbval', max, th_profbval_bin_2, termini_percentage_const)
def	profbval_maxc_maxcons1sbinth2pl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'profbval', max, th_profbval_bin_2, termini_percentage_const)


def	profbval_minc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 10, termini_percentage_const, 0)
def	profbval_minc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 20, termini_percentage_const, 0)
def	profbval_minc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 25, termini_percentage_const, 0)
def	profbval_minc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 30, termini_percentage_const, 0)
def	profbval_minc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 35, termini_percentage_const, 0)
def	profbval_minc_maxcons1sbinth1		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'profbval', min, th_profbval_bin_1, termini_percentage_const)
def	profbval_minc_maxcons1sbinth1pl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'profbval', min, th_profbval_bin_1, termini_percentage_const)
def	profbval_minc_maxcons1sbinth2		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'profbval', min, th_profbval_bin_2, termini_percentage_const)
def	profbval_minc_maxcons1sbinth2pl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'profbval', min, th_profbval_bin_2, termini_percentage_const)

def	profbval_avgc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 10, termini_percentage_const, 0)
def	profbval_avgc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 20, termini_percentage_const, 0)
def	profbval_avgc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 25, termini_percentage_const, 0)
def	profbval_avgc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 30, termini_percentage_const, 0)
def	profbval_avgc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 35, termini_percentage_const, 0)

def	profbval_avgc_maxcons1sbinth1		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'profbval', util_mean, th_profbval_bin_1, termini_percentage_const)
def	profbval_avgc_maxcons1sbinth1pl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'profbval', util_mean, th_profbval_bin_1, termini_percentage_const)

def	profbval_avgc_maxcons1sbinth2		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'profbval', util_mean, th_profbval_bin_2, termini_percentage_const)
def	profbval_avgc_maxcons1sbinth2pl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'profbval', util_mean, th_profbval_bin_2, termini_percentage_const)



def	profbval_maxc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 10, termini_percentage_const, th_asa_1)
def	profbval_maxc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 20, termini_percentage_const, th_asa_1)
def	profbval_maxc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 25, termini_percentage_const, th_asa_1)
def	profbval_maxc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 30, termini_percentage_const, th_asa_1)
def	profbval_maxc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 35, termini_percentage_const, th_asa_1)


def	profbval_minc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 10, termini_percentage_const, th_asa_1)
def	profbval_minc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 20, termini_percentage_const, th_asa_1)
def	profbval_minc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 25, termini_percentage_const, th_asa_1)
def	profbval_minc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 30, termini_percentage_const, th_asa_1)
def	profbval_minc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 35, termini_percentage_const, th_asa_1)


def	profbval_avgc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 10, termini_percentage_const, th_asa_1)
def	profbval_avgc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 20, termini_percentage_const, th_asa_1)
def	profbval_avgc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 25, termini_percentage_const, th_asa_1)
def	profbval_avgc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 30, termini_percentage_const, th_asa_1)
def	profbval_avgc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 35, termini_percentage_const, th_asa_1)

def	profbval_maxc_avgtop10p		(chains):
    return max_avg_in_top10_percent_profbval(chains)
def	profbval_minc_avgbottom10p		(chains):
    return min_avg_in_bottom10_precent_profbval(chains)



###########################################





def	profbval_maxc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 10, termini_percentage_const, th_asa_2)
def	profbval_maxc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 20, termini_percentage_const, th_asa_2)
def	profbval_maxc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 25, termini_percentage_const, th_asa_2)
def	profbval_maxc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 30, termini_percentage_const, th_asa_2)
def	profbval_maxc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', max, 35, termini_percentage_const, th_asa_2)


def	profbval_minc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 10, termini_percentage_const, th_asa_2)
def	profbval_minc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 20, termini_percentage_const, th_asa_2)
def	profbval_minc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 25, termini_percentage_const, th_asa_2)
def	profbval_minc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 30, termini_percentage_const, th_asa_2)
def	profbval_minc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', min, 35, termini_percentage_const, th_asa_2)


def	profbval_avgc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 10, termini_percentage_const, th_asa_2)
def	profbval_avgc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 20, termini_percentage_const, th_asa_2)
def	profbval_avgc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 25, termini_percentage_const, th_asa_2)
def	profbval_avgc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 30, termini_percentage_const, th_asa_2)
def	profbval_avgc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'profbval', util_mean, 35, termini_percentage_const, th_asa_2)













###################################################################################

def	iul_maxc_avg    (chains):
    return template_avg(chains,'iupred_long',max,termini_percentage_const,0)
def	iul_maxc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 10, termini_percentage_const, 0)
def	iul_maxc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 20, termini_percentage_const, 0)
def	iul_maxc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 25, termini_percentage_const, 0)
def	iul_maxc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 30, termini_percentage_const, 0)
def	iul_maxc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 35, termini_percentage_const, 0)
def	iul_maxc_frac1sbin		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', max, th_iupred_bin, termini_percentage_const, 0)
def	iul_maxc_maxcons1sbin		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'iupred_long', max, th_iupred_bin, termini_percentage_const)
def	iul_maxc_maxcons1sbinpl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'iupred_long', max, th_iupred_bin, termini_percentage_const)
def	iul_minc_avg		(chains):
    return template_avg(chains,'iupred_long',min,termini_percentage_const,0)
def	iul_minc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 10, termini_percentage_const, 0)
def	iul_minc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 20, termini_percentage_const, 0)
def	iul_minc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 25, termini_percentage_const, 0)
def	iul_minc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 30, termini_percentage_const, 0)
def	iul_minc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 35, termini_percentage_const, 0)
def	iul_minc_frac1sbin		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', min, th_iupred_bin, termini_percentage_const, 0)
def	iul_minc_maxcons1sbin		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'iupred_long', min, th_iupred_bin, termini_percentage_const)
def	iul_minc_maxcons1sbinpl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'iupred_long', min, th_iupred_bin, termini_percentage_const)
def	iul_avgc_avg		(chains):
    return template_avg(chains,'iupred_long', util_mean, termini_percentage_const, 0)
def	iul_avgc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 10, termini_percentage_const, 0)
def	iul_avgc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 20, termini_percentage_const, 0)
def	iul_avgc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 25, termini_percentage_const, 0)
def	iul_avgc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 30, termini_percentage_const, 0)
def	iul_avgc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 35, termini_percentage_const, 0)
def	iul_avgc_frac1sbin		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', util_mean, th_iupred_bin, termini_percentage_const, 0)
def	iul_avgc_maxcons1sbin		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'iupred_long', util_mean, th_iupred_bin, termini_percentage_const)
def	iul_avgc_maxcons1sbinpl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'iupred_long', util_mean, th_iupred_bin, termini_percentage_const)
def	ius_maxc_avg		(chains):
    return template_avg(chains,'iupred_short',max,termini_percentage_const,0)
def	ius_maxc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 10, termini_percentage_const, 0)
def	ius_maxc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 20, termini_percentage_const, 0)
def	ius_maxc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 25, termini_percentage_const, 0)
def	ius_maxc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 30, termini_percentage_const, 0)
def	ius_maxc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 35, termini_percentage_const, 0)
def	ius_maxc_frac1sbin		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', max, th_iupred_bin, termini_percentage_const, 0)
def	ius_maxc_maxcons1sbin		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'iupred_short', max, th_iupred_bin, termini_percentage_const)
def	ius_maxc_maxcons1sbinpl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'iupred_short', max, th_iupred_bin, termini_percentage_const)
def	ius_minc_avg		(chains):
    return template_avg(chains,'iupred_short',min,termini_percentage_const,0)
def	ius_minc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 10, termini_percentage_const, 0)
def	ius_minc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 20, termini_percentage_const, 0)
def	ius_minc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 25, termini_percentage_const, 0)
def	ius_minc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 30, termini_percentage_const, 0)
def	ius_minc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 35, termini_percentage_const, 0)
def	ius_minc_frac1sbin		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', min, th_iupred_bin, termini_percentage_const, 0)
def	ius_minc_maxcons1sbin		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'iupred_short', min, th_iupred_bin, termini_percentage_const)
def	ius_minc_maxcons1sbinpl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'iupred_short', min, th_iupred_bin, termini_percentage_const)
def	ius_avgc_avg		(chains):
    return template_avg(chains,'iupred_short', util_mean, termini_percentage_const, 0)
def	ius_avgc_maxswsum10		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 10, termini_percentage_const, 0)
def	ius_avgc_maxswsum20		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 20, termini_percentage_const, 0)
def	ius_avgc_maxswsum25		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 25, termini_percentage_const, 0)
def	ius_avgc_maxswsum30		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 30, termini_percentage_const, 0)
def	ius_avgc_maxswsum35		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 35, termini_percentage_const, 0)
def	ius_avgc_frac1sbin		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', util_mean, th_iupred_bin, termini_percentage_const, 0)
def	ius_avgc_maxcons1sbin		(chains):
    return template_max_consecutive_1s_in_bin(chains, 'iupred_short', util_mean, th_iupred_bin, termini_percentage_const)
def	ius_avgc_maxcons1sbinpl		(chains):
    return template_max_consecutive_1s_in_bin_per_length(chains, 'iupred_short', util_mean, th_iupred_bin, termini_percentage_const)
def	iul_maxc_avg_asa_th1		(chains):
    return template_avg(chains,'iupred_long',max,termini_percentage_const,th_asa_1)
def	iul_maxc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 10, termini_percentage_const, th_asa_1)
def	iul_maxc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 20, termini_percentage_const, th_asa_1)
def	iul_maxc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 25, termini_percentage_const, th_asa_1)
def	iul_maxc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 30, termini_percentage_const, th_asa_1)
def	iul_maxc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 35, termini_percentage_const, th_asa_1)
def	iul_maxc_frac1sbin_asa_th1		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', max, th_iupred_bin, termini_percentage_const, th_asa_1)
def	iul_minc_avg_asa_th1		(chains):
    return template_avg(chains,'iupred_long',min,termini_percentage_const,th_asa_1)
def	iul_minc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 10, termini_percentage_const, th_asa_1)
def	iul_minc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 20, termini_percentage_const, th_asa_1)
def	iul_minc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 25, termini_percentage_const, th_asa_1)
def	iul_minc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 30, termini_percentage_const, th_asa_1)
def	iul_minc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 35, termini_percentage_const, th_asa_1)
def	iul_minc_frac1sbin_asa_th1		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', min, th_iupred_bin, termini_percentage_const, th_asa_1)
def	iul_avgc_avg_asa_th1		(chains):
    return template_avg(chains,'iupred_long', util_mean, termini_percentage_const, th_asa_1)
def	iul_avgc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 10, termini_percentage_const, th_asa_1)
def	iul_avgc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 20, termini_percentage_const, th_asa_1)
def	iul_avgc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 25, termini_percentage_const, th_asa_1)
def	iul_avgc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 30, termini_percentage_const, th_asa_1)
def	iul_avgc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 35, termini_percentage_const, th_asa_1)
def	iul_avgc_frac1sbin_asa_th1		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', util_mean, th_iupred_bin, termini_percentage_const, th_asa_1)

def	ius_maxc_avg_asa_th1		(chains):
    return template_avg(chains,'iupred_short',max,termini_percentage_const,th_asa_1)
def	ius_maxc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 10, termini_percentage_const, th_asa_1)
def	ius_maxc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 20, termini_percentage_const, th_asa_1)
def	ius_maxc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 25, termini_percentage_const, th_asa_1)
def	ius_maxc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 30, termini_percentage_const, th_asa_1)
def	ius_maxc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 35, termini_percentage_const, th_asa_1)
def	ius_maxc_frac1sbin_asa_th1		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', max, th_iupred_bin, termini_percentage_const, th_asa_1)
def	ius_minc_avg_asa_th1		(chains):
    return template_avg(chains,'iupred_short',min,termini_percentage_const,th_asa_1)
def	ius_minc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 10, termini_percentage_const, th_asa_1)
def	ius_minc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 20, termini_percentage_const, th_asa_1)
def	ius_minc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 25, termini_percentage_const, th_asa_1)
def	ius_minc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 30, termini_percentage_const, th_asa_1)
def	ius_minc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 35, termini_percentage_const, th_asa_1)
def	ius_minc_frac1sbin_asa_th1		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', min, th_iupred_bin, termini_percentage_const, th_asa_1)

def	ius_avgc_avg_asa_th1		(chains):
    return template_avg(chains,'iupred_short', util_mean, termini_percentage_const, th_asa_1)
def	ius_avgc_maxswsum10_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 10, termini_percentage_const, th_asa_1)
def	ius_avgc_maxswsum20_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 20, termini_percentage_const, th_asa_1)
def	ius_avgc_maxswsum25_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 25, termini_percentage_const, th_asa_1)
def	ius_avgc_maxswsum30_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 30, termini_percentage_const, th_asa_1)
def	ius_avgc_maxswsum35_asa_th1		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 35, termini_percentage_const, th_asa_1)
def	ius_avgc_frac1sbin_asa_th1		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', util_mean, th_iupred_bin, termini_percentage_const, th_asa_1)



#################################################





def	iul_maxc_avg_asa_th2		(chains):
    return template_avg(chains,'iupred_long',max,termini_percentage_const,th_asa_2)
def	iul_maxc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 10, termini_percentage_const, th_asa_2)
def	iul_maxc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 20, termini_percentage_const, th_asa_2)
def	iul_maxc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 25, termini_percentage_const, th_asa_2)
def	iul_maxc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 30, termini_percentage_const, th_asa_2)
def	iul_maxc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', max, 35, termini_percentage_const, th_asa_2)
def	iul_maxc_frac1sbin_asa_th2		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', max, th_iupred_bin, termini_percentage_const, th_asa_2)
def	iul_minc_avg_asa_th2		(chains):
    return template_avg(chains,'iupred_long',min,termini_percentage_const,th_asa_2)
def	iul_minc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 10, termini_percentage_const, th_asa_2)
def	iul_minc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 20, termini_percentage_const, th_asa_2)
def	iul_minc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 25, termini_percentage_const, th_asa_2)
def	iul_minc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 30, termini_percentage_const, th_asa_2)
def	iul_minc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', min, 35, termini_percentage_const, th_asa_2)
def	iul_minc_frac1sbin_asa_th2		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', min, th_iupred_bin, termini_percentage_const, th_asa_2)
def	iul_avgc_avg_asa_th2		(chains):
    return template_avg(chains,'iupred_long', util_mean, termini_percentage_const, th_asa_2)
def	iul_avgc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 10, termini_percentage_const, th_asa_2)
def	iul_avgc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 20, termini_percentage_const, th_asa_2)
def	iul_avgc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 25, termini_percentage_const, th_asa_2)
def	iul_avgc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 30, termini_percentage_const, th_asa_2)
def	iul_avgc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_long', util_mean, 35, termini_percentage_const, th_asa_2)
def	iul_avgc_frac1sbin_asa_th2		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_long', util_mean, th_iupred_bin, termini_percentage_const, th_asa_2)

def	ius_maxc_avg_asa_th2		(chains):
    return template_avg(chains,'iupred_short',max,termini_percentage_const,th_asa_2)
def	ius_maxc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 10, termini_percentage_const, th_asa_2)
def	ius_maxc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 20, termini_percentage_const, th_asa_2)
def	ius_maxc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 25, termini_percentage_const, th_asa_2)
def	ius_maxc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 30, termini_percentage_const, th_asa_2)
def	ius_maxc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', max, 35, termini_percentage_const, th_asa_2)
def	ius_maxc_frac1sbin_asa_th2		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', max, th_iupred_bin, termini_percentage_const, th_asa_2)
def	ius_minc_avg_asa_th2		(chains):
    return template_avg(chains,'iupred_short',min,termini_percentage_const,th_asa_2)
def	ius_minc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 10, termini_percentage_const, th_asa_2)
def	ius_minc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 20, termini_percentage_const, th_asa_2)
def	ius_minc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 25, termini_percentage_const, th_asa_2)
def	ius_minc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 30, termini_percentage_const, th_asa_2)
def	ius_minc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', min, 35, termini_percentage_const, th_asa_2)
def	ius_minc_frac1sbin_asa_th2		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', min, th_iupred_bin, termini_percentage_const, th_asa_2)

def	ius_avgc_avg_asa_th2		(chains):
    return template_avg(chains,'iupred_short', util_mean, termini_percentage_const, th_asa_2)
def	ius_avgc_maxswsum10_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 10, termini_percentage_const, th_asa_2)
def	ius_avgc_maxswsum20_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 20, termini_percentage_const, th_asa_2)
def	ius_avgc_maxswsum25_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 25, termini_percentage_const, th_asa_2)
def	ius_avgc_maxswsum30_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 30, termini_percentage_const, th_asa_2)
def	ius_avgc_maxswsum35_asa_th2		(chains):
    return template_max_sliding_window(chains, 'iupred_short', util_mean, 35, termini_percentage_const, th_asa_2)
def	ius_avgc_frac1sbin_asa_th2		(chains):
    return template_frac_of_1s_in_binary(chains, 'iupred_short', util_mean, th_iupred_bin, termini_percentage_const, th_asa_2)








#################################################






def profbval_minc_minswsum10(chains):
    return template_sliding_window(chains, 'profbval', min, min, 10, termini_percentage_const, 0)
def profbval_minc_minswsum20(chains):
    return template_sliding_window(chains, 'profbval', min, min, 20, termini_percentage_const, 0)
def profbval_minc_minswsum25(chains):
    return template_sliding_window(chains, 'profbval', min, min, 25, termini_percentage_const, 0)
def profbval_minc_minswsum30(chains):
    return template_sliding_window(chains, 'profbval', min, min, 30, termini_percentage_const, 0)
def profbval_minc_minswsum35(chains):
    return template_sliding_window(chains, 'profbval', min, min, 35, termini_percentage_const, 0)
def profbval_maxc_maxminswsum10(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 10, termini_percentage_const, 0)
def profbval_maxc_maxminswsum20(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 20, termini_percentage_const, 0)
def profbval_maxc_maxminswsum25(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 25, termini_percentage_const, 0)
def profbval_maxc_maxminswsum30(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 30, termini_percentage_const, 0)
def profbval_maxc_maxminswsum35(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 35, termini_percentage_const, 0)




def profbval_minc_minswsum10_asa_th1(chains):
    return template_sliding_window(chains, 'profbval', min, min, 10, termini_percentage_const, th_asa_1)
def profbval_minc_minswsum20_asa_th1(chains):
    return template_sliding_window(chains, 'profbval', min, min, 20, termini_percentage_const, th_asa_1)
def profbval_minc_minswsum25_asa_th1(chains):
    return template_sliding_window(chains, 'profbval', min, min, 25, termini_percentage_const, th_asa_1)
def profbval_minc_minswsum30_asa_th1(chains):
    return template_sliding_window(chains, 'profbval', min, min, 30, termini_percentage_const, th_asa_1)
def profbval_minc_minswsum35_asa_th1(chains):
    return template_sliding_window(chains, 'profbval', min, min, 35, termini_percentage_const, th_asa_1)    
def profbval_maxc_maxminswsum10_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 10, termini_percentage_const, th_asa_1)
def profbval_maxc_maxminswsum20_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 20, termini_percentage_const, th_asa_1)
def profbval_maxc_maxminswsum25_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 25, termini_percentage_const, th_asa_1)
def profbval_maxc_maxminswsum30_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 30, termini_percentage_const, th_asa_1)
def profbval_maxc_maxminswsum35_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'profbval', max, 35, termini_percentage_const, th_asa_1)




def iul_maxc_maxminswsum10(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 10, termini_percentage_const, 0)
def iul_maxc_maxminswsum20(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 20, termini_percentage_const, 0)
def iul_maxc_maxminswsum25(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 25, termini_percentage_const, 0)
def iul_maxc_maxminswsum30(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 30, termini_percentage_const, 0)
def iul_maxc_maxminswsum35(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 35, termini_percentage_const, 0) 
def iul_minc_minswsum10(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 10, termini_percentage_const, 0)
def iul_minc_minswsum20(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 20, termini_percentage_const, 0)
def iul_minc_minswsum25(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 25, termini_percentage_const, 0)
def iul_minc_minswsum30(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 30, termini_percentage_const, 0)
def iul_minc_minswsum35(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 35, termini_percentage_const, 0)



def ius_maxc_maxminswsum10(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 10, termini_percentage_const, 0)
def ius_maxc_maxminswsum20(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 20, termini_percentage_const, 0)
def ius_maxc_maxminswsum25(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 25, termini_percentage_const, 0)
def ius_maxc_maxminswsum30(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 30, termini_percentage_const, 0)
def ius_maxc_maxminswsum35(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 35, termini_percentage_const, 0) 
def ius_minc_minswsum10(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 10, termini_percentage_const, 0)
def ius_minc_minswsum20(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 20, termini_percentage_const, 0)
def ius_minc_minswsum25(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 25, termini_percentage_const, 0)
def ius_minc_minswsum30(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 30, termini_percentage_const, 0)
def ius_minc_minswsum35(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 35, termini_percentage_const, 0)






def iul_maxc_maxminswsum10_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 10, termini_percentage_const, th_asa_1)
def iul_maxc_maxminswsum20_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 20, termini_percentage_const, th_asa_1)
def iul_maxc_maxminswsum25_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 25, termini_percentage_const, th_asa_1)
def iul_maxc_maxminswsum30_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 30, termini_percentage_const, th_asa_1)
def iul_maxc_maxminswsum35_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_long', min, min, 35, termini_percentage_const, th_asa_1) 
def iul_minc_minswsum10_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 10, termini_percentage_const, th_asa_1)
def iul_minc_minswsum20_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 20, termini_percentage_const, th_asa_1)
def iul_minc_minswsum25_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 25, termini_percentage_const, th_asa_1)
def iul_minc_minswsum30_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 30, termini_percentage_const, th_asa_1)
def iul_minc_minswsum35_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_long', max, 35, termini_percentage_const, th_asa_1)



def ius_maxc_maxminswsum10_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 10, termini_percentage_const, th_asa_1)
def ius_maxc_maxminswsum20_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 20, termini_percentage_const, th_asa_1)
def ius_maxc_maxminswsum25_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 25, termini_percentage_const, th_asa_1)
def ius_maxc_maxminswsum30_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 30, termini_percentage_const, th_asa_1)
def ius_maxc_maxminswsum35_asa_th1(chains):
    return template_sliding_window(chains, 'iupred_short', min, min, 35, termini_percentage_const, th_asa_1) 
def ius_minc_minswsum10_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 10, termini_percentage_const, th_asa_1)
def ius_minc_minswsum20_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 20, termini_percentage_const, th_asa_1)
def ius_minc_minswsum25_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 25, termini_percentage_const, th_asa_1)
def ius_minc_minswsum30_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 30, termini_percentage_const, th_asa_1)
def ius_minc_minswsum35_asa_th1(chains):
    return template_max_min_sliding_window(chains, 'iupred_short', max, 35, termini_percentage_const, th_asa_1)








features_to_compute = [
    frac_maxc_sidechainclass_aliphatic,
    frac_maxc_sidechainclass_basic,
    frac_maxc_sidechainclass_amide,
    frac_maxc_sidechainclass_acid,
    frac_maxc_sidechainclass_sulfurcontaining,
    frac_maxc_sidechainclass_aromatic,
    frac_maxc_sidechainpolarity_nonpolar,
    frac_maxc_sidechainpolarity_basicpolar,
    frac_maxc_sidechainpolarity_polar,
    frac_maxc_sidechainpolarity_acidicpolar,
    frac_maxc_sidechaincharge_neutral,
    frac_maxc_sidechaincharge_positive,
    frac_maxc_sidechaincharge_negative,
    frac_maxc_heavyHFRYW,
    frac_maxc_lightGASPV,
    frac_maxc_lightGASPV,
    frac_maxc_hydrophob,
    frac_maxc_tinyAGCS,
    frac_maxc_smallAGCSVTDNP,
    frac_minc_sidechainclass_aliphatic,
    frac_minc_sidechainclass_basic,
    frac_minc_sidechainclass_amide,
    frac_minc_sidechainclass_acid,
    frac_minc_sidechainclass_sulfurcontaining,
    frac_minc_sidechainclass_aromatic,
    frac_minc_sidechainpolarity_nonpolar,
    frac_minc_sidechainpolarity_basicpolar,
    frac_minc_sidechainpolarity_polar,
    frac_minc_sidechainpolarity_acidicpolar,
    frac_minc_sidechaincharge_neutral,
    frac_minc_sidechaincharge_positive,
    frac_minc_sidechaincharge_negative,
    frac_minc_heavyHFRYW,
    frac_minc_lightGASPV,
    frac_minc_lightGASPV,
    frac_minc_hydrophob,
    frac_minc_tinyAGCS,
    frac_minc_smallAGCSVTDNP,
    frac_avgc_sidechainclass_aliphatic,
    frac_avgc_sidechainclass_basic,
    frac_avgc_sidechainclass_amide,
    frac_avgc_sidechainclass_acid,
    frac_avgc_sidechainclass_sulfurcontaining,
    frac_avgc_sidechainclass_aromatic,
    frac_avgc_sidechainpolarity_nonpolar,
    frac_avgc_sidechainpolarity_basicpolar,
    frac_avgc_sidechainpolarity_polar,
    frac_avgc_sidechainpolarity_acidicpolar,
    frac_avgc_sidechaincharge_neutral,
    frac_avgc_sidechaincharge_positive,
    frac_avgc_sidechaincharge_negative,
    frac_avgc_heavyHFRYW,
    frac_avgc_lightGASPV,
    frac_avgc_lightGASPV,
    frac_avgc_hydrophob,
    frac_avgc_tinyAGCS,
    frac_avgc_smallAGCSVTDNP,
    no_unique_chains,
    sum_chain_length,
    average_chain_length,
    shortest_chain_length,
    longest_chain_length,
    num_chains,
    stdev_chains,
    # max_avg_in_top10_percent_profbval,
    # min_avg_in_bottom10_precent_profbval,
    frac_maxc_top3Flexbilityidx,
    frac_maxc_top3bvalueidx,
    frac_maxc_top3topidp,
    frac_maxc_top3foldunfold,
    frac_maxc_top3disprot,
    frac_maxc_bottom3Flexbilityidx,
    frac_maxc_bottom3bvalueidx,
    frac_maxc_bottom3topidp,
    frac_maxc_bottom3foldunfold,
    frac_maxc_bottom3disprot,
    frac_minc_top3Flexbilityidx,
    frac_minc_top3bvalueidx,
    frac_minc_top3topidp,
    frac_minc_top3foldunfold,
    frac_minc_top3disprot,
    frac_minc_bottom3Flexbilityidx,
    frac_minc_bottom3bvalueidx,
    frac_minc_bottom3topidp,
    frac_minc_bottom3foldunfold,
    frac_minc_bottom3disprot,
    frac_avgc_top3Flexbilityidx,
    frac_avgc_top3bvalueidx,
    frac_avgc_top3topidp,
    frac_avgc_top3foldunfold,
    frac_avgc_top3disprot,
    frac_avgc_bottom3Flexbilityidx,
    frac_avgc_bottom3bvalueidx,
    frac_avgc_bottom3topidp,
    frac_avgc_bottom3foldunfold,
    frac_avgc_bottom3disprot,
    frac_maxc_AA_A,
    frac_maxc_AA_R,
    frac_maxc_AA_N,
    frac_maxc_AA_D,
    frac_maxc_AA_C,
    frac_maxc_AA_Q,
    frac_maxc_AA_E,
    frac_maxc_AA_G,
    frac_maxc_AA_H,
    frac_maxc_AA_I,
    frac_maxc_AA_L,
    frac_maxc_AA_K,
    frac_maxc_AA_M,
    frac_maxc_AA_F,
    frac_maxc_AA_P,
    frac_maxc_AA_S,
    frac_maxc_AA_T,
    frac_maxc_AA_W,
    frac_maxc_AA_Y,
    frac_maxc_AA_V,
    frac_minc_AA_A,
    frac_minc_AA_R,
    frac_minc_AA_N,
    frac_minc_AA_D,
    frac_minc_AA_C,
    frac_minc_AA_Q,
    frac_minc_AA_E,
    frac_minc_AA_G,
    frac_minc_AA_H,
    frac_minc_AA_I,
    frac_minc_AA_L,
    frac_minc_AA_K,
    frac_minc_AA_M,
    frac_minc_AA_F,
    frac_minc_AA_P,
    frac_minc_AA_S,
    frac_minc_AA_T,
    frac_minc_AA_W,
    frac_minc_AA_Y,
    frac_minc_AA_V,
    frac_avgc_AA_A,
    frac_avgc_AA_R,
    frac_avgc_AA_N,
    frac_avgc_AA_D,
    frac_avgc_AA_C,
    frac_avgc_AA_Q,
    frac_avgc_AA_E,
    frac_avgc_AA_G,
    frac_avgc_AA_H,
    frac_avgc_AA_I,
    frac_avgc_AA_L,
    frac_avgc_AA_K,
    frac_avgc_AA_M,
    frac_avgc_AA_F,
    frac_avgc_AA_P,
    frac_avgc_AA_S,
    frac_avgc_AA_T,
    frac_avgc_AA_W,
    frac_avgc_AA_Y,
    frac_avgc_AA_V,
    idx_maxc_maxswnum20_top3_flex,
    idx_maxc_maxswnum20_top3_bvalue,
    idx_maxc_maxswnum20_top3_topidp,
    idx_maxc_maxswnum20_top3_foldunfold,
    idx_maxc_maxswnum20_top3_disprot,
    # profbval_maxc_maxswsum10,
    # profbval_maxc_maxswsum20,
    # profbval_maxc_maxswsum25,
    # profbval_maxc_maxswsum30,
    # profbval_maxc_maxswsum35,
    # profbval_maxc_maxcons1sbinth1,
    # profbval_maxc_maxcons1sbinth1pl,
    # profbval_maxc_maxcons1sbinth2,
    # profbval_maxc_maxcons1sbinth2pl,
    # profbval_minc_maxswsum10,
    # profbval_minc_maxswsum20,
    # profbval_minc_maxswsum25,
    # profbval_minc_maxswsum30,
    # profbval_minc_maxswsum35,
    # profbval_minc_maxcons1sbinth1,
    # profbval_minc_maxcons1sbinth1pl,
    # profbval_minc_maxcons1sbinth2,
    # profbval_minc_maxcons1sbinth2pl,
    # profbval_avgc_maxswsum10,
    # profbval_avgc_maxswsum20,
    # profbval_avgc_maxswsum25,
    # profbval_avgc_maxswsum30,
    # profbval_avgc_maxswsum35,
    # profbval_avgc_maxcons1sbinth1,
    # profbval_avgc_maxcons1sbinth1pl,
    # profbval_avgc_maxcons1sbinth2,
    # profbval_avgc_maxcons1sbinth2pl,
    # profbval_maxc_maxswsum10_asa_th1,
    # profbval_maxc_maxswsum20_asa_th1,
    # profbval_maxc_maxswsum25_asa_th1,
    # profbval_maxc_maxswsum30_asa_th1,
    # profbval_maxc_maxswsum35_asa_th1,
    # profbval_minc_maxswsum10_asa_th1,
    # profbval_minc_maxswsum20_asa_th1,
    # profbval_minc_maxswsum25_asa_th1,
    # profbval_minc_maxswsum30_asa_th1,
    # profbval_minc_maxswsum35_asa_th1,
    # profbval_avgc_maxswsum10_asa_th1,
    # profbval_avgc_maxswsum20_asa_th1,
    # profbval_avgc_maxswsum25_asa_th1,
    # profbval_avgc_maxswsum30_asa_th1,
    # profbval_avgc_maxswsum35_asa_th1,
    # profbval_maxc_avgtop10p,
    # profbval_minc_avgbottom10p,
    # profbval_maxc_maxswsum10_asa_th2,
    # profbval_maxc_maxswsum20_asa_th2,
    # profbval_maxc_maxswsum25_asa_th2,
    # profbval_maxc_maxswsum30_asa_th2,
    # profbval_maxc_maxswsum35_asa_th2,
    # profbval_minc_maxswsum10_asa_th2,
    # profbval_minc_maxswsum20_asa_th2,
    # profbval_minc_maxswsum25_asa_th2,
    # profbval_minc_maxswsum30_asa_th2,
    # profbval_minc_maxswsum35_asa_th2,
    # profbval_avgc_maxswsum10_asa_th2,
    # profbval_avgc_maxswsum20_asa_th2,
    # profbval_avgc_maxswsum25_asa_th2,
    # profbval_avgc_maxswsum30_asa_th2,
    # profbval_avgc_maxswsum35_asa_th2,
    iul_maxc_avg,
    iul_maxc_maxswsum10,
    iul_maxc_maxswsum20,
    iul_maxc_maxswsum25,
    iul_maxc_maxswsum30,
    iul_maxc_maxswsum35,
    iul_maxc_frac1sbin,
    iul_maxc_maxcons1sbin,
    iul_maxc_maxcons1sbinpl,
    iul_minc_avg,
    iul_minc_maxswsum10,
    iul_minc_maxswsum20,
    iul_minc_maxswsum25,
    iul_minc_maxswsum30,
    iul_minc_maxswsum35,
    iul_minc_frac1sbin,
    iul_minc_maxcons1sbin,
    iul_minc_maxcons1sbinpl,
    iul_avgc_avg,
    iul_avgc_maxswsum10,
    iul_avgc_maxswsum20,
    iul_avgc_maxswsum25,
    iul_avgc_maxswsum30,
    iul_avgc_maxswsum35,
    iul_avgc_frac1sbin,
    iul_avgc_maxcons1sbin,
    iul_avgc_maxcons1sbinpl,
    ius_maxc_avg,
    ius_maxc_maxswsum10,
    ius_maxc_maxswsum20,
    ius_maxc_maxswsum25,
    ius_maxc_maxswsum30,
    ius_maxc_maxswsum35,
    ius_maxc_frac1sbin,
    ius_maxc_maxcons1sbin,
    ius_maxc_maxcons1sbinpl,
    ius_minc_avg,
    ius_minc_maxswsum10,
    ius_minc_maxswsum20,
    ius_minc_maxswsum25,
    ius_minc_maxswsum30,
    ius_minc_maxswsum35,
    ius_minc_frac1sbin,
    ius_minc_maxcons1sbin,
    ius_minc_maxcons1sbinpl,
    ius_avgc_avg,
    ius_avgc_maxswsum10,
    ius_avgc_maxswsum20,
    ius_avgc_maxswsum25,
    ius_avgc_maxswsum30,
    ius_avgc_maxswsum35,
    ius_avgc_frac1sbin,
    ius_avgc_maxcons1sbin,
    ius_avgc_maxcons1sbinpl,
    iul_maxc_avg_asa_th1,
    iul_maxc_maxswsum10_asa_th1,
    iul_maxc_maxswsum20_asa_th1,
    iul_maxc_maxswsum25_asa_th1,
    iul_maxc_maxswsum30_asa_th1,
    iul_maxc_maxswsum35_asa_th1,
    iul_maxc_frac1sbin_asa_th1,
    iul_minc_avg_asa_th1,
    iul_minc_maxswsum10_asa_th1,
    iul_minc_maxswsum20_asa_th1,
    iul_minc_maxswsum25_asa_th1,
    iul_minc_maxswsum30_asa_th1,
    iul_minc_maxswsum35_asa_th1,
    iul_minc_frac1sbin_asa_th1,
    iul_avgc_avg_asa_th1,
    iul_avgc_maxswsum10_asa_th1,
    iul_avgc_maxswsum20_asa_th1,
    iul_avgc_maxswsum25_asa_th1,
    iul_avgc_maxswsum30_asa_th1,
    iul_avgc_maxswsum35_asa_th1,
    iul_avgc_frac1sbin_asa_th1,
    ius_maxc_avg_asa_th1,
    ius_maxc_maxswsum10_asa_th1,
    ius_maxc_maxswsum20_asa_th1,
    ius_maxc_maxswsum25_asa_th1,
    ius_maxc_maxswsum30_asa_th1,
    ius_maxc_maxswsum35_asa_th1,
    ius_maxc_frac1sbin_asa_th1,
    ius_minc_avg_asa_th1,
    ius_minc_maxswsum10_asa_th1,
    ius_minc_maxswsum20_asa_th1,
    ius_minc_maxswsum25_asa_th1,
    ius_minc_maxswsum30_asa_th1,
    ius_minc_maxswsum35_asa_th1,
    ius_minc_frac1sbin_asa_th1,
    ius_avgc_avg_asa_th1,
    ius_avgc_maxswsum10_asa_th1,
    ius_avgc_maxswsum20_asa_th1,
    ius_avgc_maxswsum25_asa_th1,
    ius_avgc_maxswsum30_asa_th1,
    ius_avgc_maxswsum35_asa_th1,
    ius_avgc_frac1sbin_asa_th1,
    iul_maxc_avg_asa_th2,
    iul_maxc_maxswsum10_asa_th2,
    iul_maxc_maxswsum20_asa_th2,
    iul_maxc_maxswsum25_asa_th2,
    iul_maxc_maxswsum30_asa_th2,
    iul_maxc_maxswsum35_asa_th2,
    iul_maxc_frac1sbin_asa_th2,
    iul_minc_avg_asa_th2,
    iul_minc_maxswsum10_asa_th2,
    iul_minc_maxswsum20_asa_th2,
    iul_minc_maxswsum25_asa_th2,
    iul_minc_maxswsum30_asa_th2,
    iul_minc_maxswsum35_asa_th2,
    iul_minc_frac1sbin_asa_th2,
    iul_avgc_avg_asa_th2,
    iul_avgc_maxswsum10_asa_th2,
    iul_avgc_maxswsum20_asa_th2,
    iul_avgc_maxswsum25_asa_th2,
    iul_avgc_maxswsum30_asa_th2,
    iul_avgc_maxswsum35_asa_th2,
    iul_avgc_frac1sbin_asa_th2,
    ius_maxc_avg_asa_th2,
    ius_maxc_maxswsum10_asa_th2,
    ius_maxc_maxswsum20_asa_th2,
    ius_maxc_maxswsum25_asa_th2,
    ius_maxc_maxswsum30_asa_th2,
    ius_maxc_maxswsum35_asa_th2,
    ius_maxc_frac1sbin_asa_th2,
    ius_minc_avg_asa_th2,
    ius_minc_maxswsum10_asa_th2,
    ius_minc_maxswsum20_asa_th2,
    ius_minc_maxswsum25_asa_th2,
    ius_minc_maxswsum30_asa_th2,
    ius_minc_maxswsum35_asa_th2,
    ius_minc_frac1sbin_asa_th2,
    ius_avgc_avg_asa_th2,
    ius_avgc_maxswsum10_asa_th2,
    ius_avgc_maxswsum20_asa_th2,
    ius_avgc_maxswsum25_asa_th2,
    ius_avgc_maxswsum30_asa_th2,
    ius_avgc_maxswsum35_asa_th2,
    ius_avgc_frac1sbin_asa_th2,
    # profbval_minc_minswsum10,
    # profbval_minc_minswsum20,
    # profbval_minc_minswsum25,
    # profbval_minc_minswsum30,
    # profbval_minc_minswsum35,
    # profbval_maxc_maxminswsum10,
    # profbval_maxc_maxminswsum20,
    # profbval_maxc_maxminswsum25,
    # profbval_maxc_maxminswsum30,
    # profbval_maxc_maxminswsum35,
    # profbval_minc_minswsum10_asa_th1,
    # profbval_minc_minswsum20_asa_th1,
    # profbval_minc_minswsum25_asa_th1,
    # profbval_minc_minswsum30_asa_th1,
    # profbval_minc_minswsum35_asa_th1,
    # profbval_maxc_maxminswsum10_asa_th1,
    # profbval_maxc_maxminswsum20_asa_th1,
    # profbval_maxc_maxminswsum25_asa_th1,
    # profbval_maxc_maxminswsum30_asa_th1,
    # profbval_maxc_maxminswsum35_asa_th1,
    iul_maxc_maxminswsum10,
    iul_maxc_maxminswsum20,
    iul_maxc_maxminswsum25,
    iul_maxc_maxminswsum30,
    iul_maxc_maxminswsum35,
    iul_minc_minswsum10,
    iul_minc_minswsum20,
    iul_minc_minswsum25,
    iul_minc_minswsum30,
    iul_minc_minswsum35,
    ius_maxc_maxminswsum10,
    ius_maxc_maxminswsum20,
    ius_maxc_maxminswsum25,
    ius_maxc_maxminswsum30,
    ius_maxc_maxminswsum35,
    ius_minc_minswsum10,
    ius_minc_minswsum20,
    ius_minc_minswsum25,
    ius_minc_minswsum30,
    ius_minc_minswsum35,
    iul_maxc_maxminswsum10_asa_th1,
    iul_maxc_maxminswsum20_asa_th1,
    iul_maxc_maxminswsum25_asa_th1,
    iul_maxc_maxminswsum30_asa_th1,
    iul_maxc_maxminswsum35_asa_th1,
    iul_minc_minswsum10_asa_th1,
    iul_minc_minswsum20_asa_th1,
    iul_minc_minswsum25_asa_th1,
    iul_minc_minswsum30_asa_th1,
    iul_minc_minswsum35_asa_th1,
    ius_maxc_maxminswsum10_asa_th1,
    ius_maxc_maxminswsum20_asa_th1,
    ius_maxc_maxminswsum25_asa_th1,
    ius_maxc_maxminswsum30_asa_th1,
    ius_maxc_maxminswsum35_asa_th1,
    ius_minc_minswsum10_asa_th1,
    ius_minc_minswsum20_asa_th1,
    ius_minc_minswsum25_asa_th1,
    ius_minc_minswsum30_asa_th1,
    ius_minc_minswsum35_asa_th1
]


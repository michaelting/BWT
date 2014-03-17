#!/usr/bin/env python

"""
#=======================================================#
# Burrows-Wheeler Transform Test Code                   #
# Michael Ting                                          #
# 10 October 2013                                       #
#=======================================================#
"""

import sys
sys.path.append("./tests/")
from bwt import *
from nose.tools import ok_, eq_, raises, nottest

def test_process_FASTA():
    """
    test_process_FASTA
    
    Ensures that FASTA files are read correctly
    """
    testfile1 = open("./tests/test1.fasta","r")
    seqgen = process_FASTA(testfile1)
    
    for name, seq in seqgen:
        eq_(name, ">test1", "FASTA parser name test 1 failed!")
        eq_(seq, "ACGT", "FASTA parser seq test 1 failed!")
    
    testfile1.close()

def test_prep_DC3():
    """
    test_prep_DC3
    
    Ensures that sequences are correctly converted to numerical lists of ascii values
    """
    s = "ABCabc"
    t = prep_DC3(s)
    correct =[65, 66, 67, 97, 98, 99]

    eq_(t, correct, "Prep DC3 test failed!")
 
def test_DC3():
    """
    test_DC3
    """

    seq = "ACGTACGTACGT$"
    numlst = prep_DC3(seq)
    suffixarray = DC3(numlst)
    correct = {0: 12, 1: 8, 2: 4, 3: 0, 4: 9, 5: 5, 6: 1, 7: 10, 8: 6, 9: 2, 10: 11, 11: 7, 12: 3}
    
    eq_(suffixarray, correct, "DC3 suffix array output not correct!")
    
def test_DC3_sample():
    """
    test_DC3_sample
    """
    
    numlst = [0,1,2,3,4,5,6,7]
    B0, B1, B2, C = DC3_sample(numlst)
    correctB0 = [0,3,6]
    correctB1 = [1,4,7]
    correctB2 = [2,5]
    correctC = [1,4,7,2,5]
    
    eq_(B0, correctB0, "B0 values not correct!")
    eq_(B1, correctB1, "B1 values not correct!")
    eq_(B2, correctB2, "B2 values not correct!")
    eq_(C, correctC, "C values not correct!")
    
def test_DC3_num_to_suffix():
    """
    test_DC3_num_to_suffix
    """
    numR = [1,1,2,3,4,5,6]
    samind = [1,4,2,5]
    suffixes = DC3_num_to_suffix(numR, samind)
    correct = [([1, 2, 3],1), ([4, 5, 6],4), ([2, 3, 4],2), ([5, 6, 0],5)]
    
    eq_(suffixes, correct, "DC3 numbers to numeral suffix not correct!")
    
def test_DC3_samsort():
    """
    test_DC3_samsort
    """
    """
    seq = "ACGTA"
    samind = [1,4,2]
    # samples = ['CGT','A$$','GTA']
    sortedsam = ['A$$','CGT','GTA']
    ordered, result = DC3_samsort(seq, samind)
    
    eq_(result, sortedsam, "Radix-Sorted B1UB2 samples from indices not correct!")
    
    seq2 = "ACGTACGTA"
    samind2 = [1,4,7,2,5,8]
    # samples = ['CGT,'ACG','TA$','GTA','CGT','A$$']
    sortedsam2 = ['A$$','ACG','CGT','CGT','GTA','TA$']
    ordered = result2 = DC3_samsort(seq2, samind2)
    
    print ""
    print result2
    print sortedsam2
    
    eq_(result2, sortedsam2, "Radix-Sorted B1UB2 samples 2 from indices not correct!")
    """
    
    s = [([1, 2, 3],1), ([4, 5, 6],4), ([2, 3, 4],2), ([5, 6, 0],5)]
    ordresult1, indrmapresult1 = DC3_samsort(s)
    correctord1 = [([1, 2, 3], 1), ([2, 3, 4], 2), ([4, 5, 6], 4), ([5, 6, 0], 5)]
    correctmap1 = {1: 1, 2: 2, 4: 3, 5: 4}
    
    eq_(ordresult1, correctord1, "Radix-sorted B1UB2 samples not correct!")
    eq_(indrmapresult1, correctmap1, "Radix-sorted mappings not correct!")
    
def test_DC3_radixsort():
    """
    test_DC3_radixsort
    """
    
    """
    samlst = ['AAA','ACA','C$$','bbA','CCC','A$$','aCC','cgG','cGG','bba']
    result = DC3_radixsort(samlst, 0)
    correct = ['A$$','AAA','ACA','C$$','CCC','aCC','bbA','bba','cGG','cgG']
    """
    
    #s = [[1,2,3],[1,2,4],[5,2,1],[3,3,3],[2,4,1],[3,2,4],[1,2,3]]
    s = [([1, 2, 3],1), ([4, 5, 6],4), ([2, 3, 4],2), ([5, 6, 0],5)]
    result = DC3_radixsort(s, 0)
    correct = [([1, 2, 3],1), ([2, 3, 4], 2), ([4, 5, 6],4), ([5, 6, 0],5)]

    eq_(result, correct, "Radix-Sorted sampled list not correct!")
    

def test_DC3_check_radix():
    """
    test_DC3_check_radix
    """
    
    """
    sortedsamlst = ['A$$','AAA','ABC','CCC']
    flag, newstr, checked = DC3_check_radix(sortedsamlst)
    eq_(flag, False, "DC3 check radix boolean return incorrect!")
    eq_(newstr, [1,2,3,4], "DC3 check radix new string return incorrect!")
    
    sortedsamlst.append('AAA')
    flag2, newstr2, checked2 = DC3_check_radix(sortedsamlst)
    eq_(flag2, True, "DC3 check radix boolean 2 return incorrect!")
    eq_(newstr2, 'ABCDB', "DC3 check radix new string 2 return incorrect!")
    """
    
    """
    s1 = [[1, 2, 3], [1, 2, 4], [2, 4, 1], [3, 2, 4], [3, 3, 3], [5, 2, 1]]
    flag1, newR1, checked1 = DC3_check_radix(s1)
    eq_(flag1, False, "DC2 check radix boolean return 1 incorrect!")
    eq_(newR1, [1,2,3,4,5,6], "DC3 check radix new R return 1 incorrect!")
    
    s2 = [[1, 2, 3], [1, 2, 3], [1, 2, 4], [2, 4, 1], [3, 2, 4], [3, 3, 3], [5, 2, 1]]
    flag2, newR2, checked2 = DC3_check_radix(s2)
    eq_(flag2, True, "DC3 check radix boolean return 2 incorrect!")
    eq_(newR2, [1, 1, 2, 3, 4, 5, 6], "DC3 check radix new R return 2 incorrect!")
    """
    
    sortedtuplst1 = [([1, 2, 3],1), ([2, 3, 4], 2), ([4, 5, 6],4), ([5, 6, 0],5)]
    flag1, indexrankmap1, checked1 = DC3_check_radix(sortedtuplst1)
    correctindexrankmap1 = {1: 1, 2: 2, 4: 3, 5: 4}
    eq_(flag1, False, "DC3 check radix boolean return 1 incorrect!")
    for key in correctindexrankmap1.keys():
        eq_(indexrankmap1[key], correctindexrankmap1[key], "DC3 check radix indexrankmap return 1 incorrect!")
        
    sortedtuplst2 = [([1, 2, 3],1), ([1, 2, 3],7), ([2, 3, 4], 2), ([4, 5, 6],4), ([5, 6, 0],5)]
    flag2, indexrankmap2, checked2 = DC3_check_radix(sortedtuplst2)
    correctindexrankmap2 = {1: 1, 2: 2, 4: 3, 5: 4, 7: 1}
    eq_(flag2, True, "DC3 check radix boolean return 1 incorrect!")
    for key in correctindexrankmap2.keys():
        eq_(indexrankmap2[key], correctindexrankmap2[key], "DC3 check radix indexrankmap return 1 incorrect!")


def test_DC3_nonsamsort():
    """
    test_DC3_nonsamsort
    """
    
    indexrankmap = {1:1,
                    2:4,
                    4:2,
                    5:6,
                    7:5,
                    8:3,
                    10:7,
                    11:8,
                    13:0,
                    14:0}

    nonsortedsamtups = [([121,97,98], 0), ([98, 97, 100], 3), ([97, 98, 98], 6), ([97, 100, 111], 9) ]
    numlst = [121, 97, 98, 98, 97, 100, 97, 98, 98, 97, 100, 111]
    result = DC3_nonsamsort(nonsortedsamtups, numlst, indexrankmap)
    correct = [((97, 5), 6), ((97, 7), 9), ((98, 2), 3), ((121, 1), 0)]
    
    eq_(result, correct, "DC3 Non-sample sorting step incorrect!")
        
def test_DC3_nonsamradixsort():
    """
    test_DC3_nonsamradixsort
    """
    
    nonsamtups = [((1,7),9),((1,5),6),((0,0),12),((2,2),3),((3,1),0)]
    iternum = 0
    result = DC3_nonsamradixsort(nonsamtups, iternum)
    correct = [((0, 0), 12), ((1, 5), 6), ((1, 7), 9), ((2, 2), 3), ((3, 1), 0)]
    
    eq_(result, correct, "DC3 Non-sample radix sort incorrect!")

@nottest
def test_DC3_mergesamples():
    """
    test_DC3_mergesamples
    """    
    
def test_encode_BWT():
    """
    test_encode_BWT
    """
    s = "ACGT$"
    correctL = 'T$ACG'
    BWT = encode_BWT(s)
    
    eq_(BWT, correctL, "BWT encoding not correct!")
    
def test_decode_BWT():
    """
    test_decode_BWT
    """

def test_get_M():
    """
    test_get_M
    """
    s = "hello$"
    BWT = encode_BWT(s)
    M = get_M(BWT)
    correct = {'h': 2, 'e': 1, '$': 0, 'o': 5, 'l': 3}
    
    eq_(M, correct, "M table retrival incorrect!")

def test_get_occurrences():
    """
    test_get_occurrences
    """
    
    s = "hello$"
    BWT = encode_BWT(s)
    occ = get_occurrences(BWT)
    correct = {0: Counter({'o': 1}),
               1: Counter({'h': 1, 'o': 1}),
               2: Counter({'h': 1, '$': 1, 'o': 1}),
               3: Counter({'h': 1, 'e': 1, '$': 1, 'o': 1}),
               4: Counter({'h': 1, 'e': 1, '$': 1, 'o': 1, 'l': 1}),
               5: Counter({'l': 2, 'h': 1, 'e': 1, '$': 1, 'o': 1})}
    eq_(occ, correct, "Occurrence table calculation incorrect!")

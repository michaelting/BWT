#!/usr/bin/env python

"""
#=======================================================#
# Burrows-Wheeler Transform                             #
# Michael Ting                                          #
# 10 October 2013                                       #
#=======================================================#
"""

from argparse import ArgumentParser
from collections import Counter
import pprint

def process_FASTA(infile):
    """
    Generator that reads a FASTA file and outputs line by line
    """
    name, seq = None, []
    for line in infile:
        if line.startswith(">"):            # sequence header found
            if name:                        # next sequence in file found
                yield (name, ''.join(seq))  # output the sequence
            name, seq = line.strip(), []    # create a new sequence item
        else:
            seq.append(line.strip())        # continuing portion of sequence appended
    if name: 
        yield (name, ''.join(seq))	

def prep_DC3(charseq):
    """
    Prepare for DC3 to run on numerical lists since our input is a sequence
    of characters
    Input:
        charseq - our sequence of interest from our FASTA file
    Output:
        numlst  - the numerical representation of our characters
    """
    # stored as [0,1,]
    numlst = []
    # convert chars to ascii numerical values
    for i in range(len(charseq)):
        char = charseq[i]
        num = int(ord(char))    # ensure value stored as an int
        numlst.append(num)
    return numlst

def DC3(numlst):
    """
    DC3/Skew Algorithm
    Computes the suffix array S for the string seq
    Input:
        numlst - the sequence of interest represented by a numerical list
    Output:
        suffixarray  - the suffix array of the input sequence numlst
    """
    B0, B1, B2, C = DC3_sample(numlst) # {C} = {B1}U{B2}
    
    # extract suffixes and indices of suffix start positions
    samsuftups = DC3_num_to_suffix(numlst, C)
    nonsamsuftups = DC3_num_to_suffix(numlst, B0)
    
    # sort sample suffixes
    sortedsamples, indexrankmap = DC3_samsort(samsuftups)
    
    # pad rank(S_n+1) and rank(S_n+2) with 0's
    n = len(numlst)
    indexrankmap[n]     = 0
    indexrankmap[n+1]   = 0
    
    # pad the ends
    numlst.extend([0,0])
    
    # sort the nonsample suffixes using the indexrankmap
    sortednonsamples = DC3_nonsamsort(nonsamsuftups, numlst, indexrankmap)
    
    # merge the sets of sorted suffixes into one big sorted list
    allsortedsuffixes, suffixarray = DC3_mergesamples(numlst, sortednonsamples, sortedsamples, indexrankmap)
    
    return suffixarray

    
def DC3_sample(numlst):
    """
    Sort the indices and return the sample indices for B0, B1, B2, and B1UB2
    Input:
        numlst - the sequence of interest in numerical list representation
    Output:
        B0  - the set of indices i for which i = 0mod3
        B1  - the set of indices i for which i = 1mod3
        B2  - the set of indices i for which i = 2mod3
        C   - the set of indices that is the union of {B1} and {B2}
    """
    B0, B1, B2 = [], [], []
    # sort our indices modulo 3 into our samples
    for i in range(len(numlst)):
        if i % 3 == 1:
            B1.append(i)
        elif i % 3 == 2:
            B2.append(i)
        elif i % 3 == 0:
            B0.append(i)
        else:
            raise Exception("Invalid index found: %d" % i)
    # join our difference cover sample positions
    C = B1[:]       # copy list by value, not by reference
    C.extend(B2)    # join B1 and B2 by value
    return B0, B1, B2, C
    
def DC3_num_to_suffix(numR, samind):
    """
    Convert numerical sequence representation to suffixes for recursions
    Input:
        numR    - our new string R represented numerically as a list of ints
        samind  - sample indices corresponding to R1'UR2'
    Output:
        suffixes     - a list of tuples, with each tuple appearing:
                    ([65,2,1],1), ([34,0,0],4)
                    So the list looks like
                    [([65,2,1],1), ([34,0,0],4), ...]
        Tuples contain (suffix, index) where index corresponds to the position
        of the beginning of the suffix in the original number list.
    """
    # suffixes looks like [[1,3,2],[4,5,0],...]
    suffixes = []
    for index in samind:
        sufnumlst = numR[index:index+3]    # extract suffix of length 2 starting at index
        while len(sufnumlst) < 3:   # pad short suffixes with 0's
            sufnumlst.append(0)
        # each tuple stores the suffix and the index it came from in the original numlist
        sufnumtuple = (sufnumlst, index)
        suffixes.append(sufnumtuple)
    return suffixes
    
def DC3_samsort(suftuplst):
    """
    Sort the samples from C = B1UB2
    Input:  
        suftuplst   - a list of tuples of (suffix, origindex)
    Output:
        ordered     - list of ordered suffixes
        indexrankmap    - a mapping of original string indices i to the suffix rank, rank(S_i)
    """

    ITER = 0
    
    # lexicographically sorted suffix tuples (suffix, index)
    ordered = DC3_radixsort(suftuplst, ITER)
    
    # determine whether we need to recurse on R' (flag)
    # get a mapping of original suffix index to its rank {index:rank}
    # checked contains exactly one copy of a suffix (elim repeats)
    flag, indexrankmap, checked = DC3_check_radix(ordered)

    # create new string R' using mappings from radix sort on original suffixes
    newR = []
    for suftup in suftuplst:
        suf = suftup[0]     # suffix
        sufind = suftup[1]  # index
        newR.append(indexrankmap[sufind])
    
    # need to recurse on newR since a duplicate suffix is found in the suffixes
    while flag:
        # DC3 should return a list of sorted suffixes
        suffixarray = DC3(newR)  #suffixarray = rank:index
        
        invsuffixarray = dict((v,k) for k,v in suffixarray.items()) # index: rank
        # sort by ranks
        indices = sorted(invsuffixarray, cmp=lambda x,y:cmp(x,y))
        # get ordering of suffixes using suffix array information
        # relabel our suffixes using our ranks, which are now unique from recursion
        arraysortedsamples = []
        for index in indices:
            rank = invsuffixarray[index]
            oldtup = suftuplst[index]
            newtup = ([rank], oldtup[1])  # rank becomes the new suffix name
            arraysortedsamples.append(newtup)
            
        ordered = DC3_radixsort(arraysortedsamples, ITER)
        flag, indexrankmap, checked = DC3_check_radix(ordered)
    
    return ordered, indexrankmap
    
def DC3_radixsort(suftuplst, iternum):
    """
    Sort the suffixes in linear time using radix sort
    Input:
        samlst  - a list of samples (triplets of indices)
        iternum - the iteration number, valid from 0 to 2, base case at 3
    Output:
        ordered - a list of tuples of numerically ordered samples [ ([1,2,3],1), ([2,4,5],5),...]
    """    
    # stop at iternum = 3 since we are in mod 3 and suffixes are max length 3
    if iternum > 2:
        return suftuplst
    
    # buckets have form {1:[1,2,3],[1,7,2],...]
    #                    3:[3,2,5],[3,18,20],...}
    buckets = {}
    # process all sample suffixes using radix sort
    for suftuple in suftuplst:
        # sort suffixes into buckets by character at index iternum {0, 1, 2}
        sufnumlst = suftuple[0]
        label = sufnumlst[iternum]
        # create bucket if it does not yet exist
        if label not in buckets:
            buckets[label] = []
        # buckets will look like [ ([1,2,3],1), ([1,2,4],3), ...]]
        buckets[label].append(suftuple)
    # sort bucket keys numerically increasing, 1,2,3,...
    keysortedbuckets = sorted(buckets, cmp=lambda x,y:cmp(x,y)) # sort by numerical value
    # order the suffixes
    ordered = []
    for key in keysortedbuckets:
        # get list of tuples in the bucket
        samples = buckets[key]
        # recurse on buckets with multiple tuples
        if len(buckets[key]) > 1:   # more than one tuple in a bucket
            samples = DC3_radixsort(samples, iternum+1) # use next indexed num in suffix
        ordered.extend(samples)
    return ordered

def DC3_check_radix(sortedtuplst):
    """
    Determines whether or not we need to recurse with R'
    Input:
        sortedtuplst    - result from running DC3_radixsort
    Output:
        flag    - True if we need to recurse with DC3 on R' (duplicates found)
                - False if all characters are different, no need to recurse
sortedsamranks  - a list of ranks mapping on the same indices to sortedsamlst
        indexrankmap    - a dictionary mapping original suffix index to its rank
        checked - a dictionary of samnumlist:rank, {(1,2,3):1,(2,3,4):2,...}
    """
    checked = {}
    indexrankmap = {}
    rank = 1
    flag = False
    # rank samples from the radix-sorted list of tuples
    for samnumtuple in sortedtuplst:
        samnumlst = samnumtuple[0]  # extract numlst from tuple
        sufindex = samnumtuple[1]   # index of suffix from original string
        samnumlst = tuple(samnumlst)    # need immutable keys       
        # sample is unique
        if samnumlst not in checked:
            checked[samnumlst] = rank
            #sortedsamranks.append(rank)
            indexrankmap[sufindex] = rank
            rank += 1
        # sample is a repeat
        else:
            flag = True
            #sortedsamranks.append(checked[samnumlst])
            indexrankmap[sufindex] = checked[samnumlst] # rank tie, use from table
    return flag, indexrankmap, checked

def DC3_nonsamsort(nonsortedsamtups, numlst, indexrankmap):
    """
    Sort the non-sample suffixes using information from the sorted sample suffixes
    Input:
        nonsortedsamtups    - a list of tuples (suffix, origindex) from B0
        numlst      - the original list of numbers representing our sequence
        indexrankmap    - a mapping from indices i in our original numlst to suffix rank, rank(S_i)
    Output:
        A non-sample indexrankmap dictionary
        OR a sorted list of tuples of non-sample suffixes and their indices
    """
    
    nonsampairtuplst = []
    for samtup in nonsortedsamtups:
        samlst = samtup[0]
        samind = samtup[1]
        # create the pair (t_i, rank(S_i+1))
        t_i = samlst[0] # the first num in our suffix        
        trank = indexrankmap[samind+1]  # rank(S_i+1)
        nonsampair = (t_i, trank)
        nonsampairtup = (nonsampair, samind)    # keep track of original index
        nonsampairtuplst.append(nonsampairtup)
    # radix sort the non-sample suffix tuples
    sortednonsampairtuplst = DC3_nonsamradixsort(nonsampairtuplst, 0)    
    return sortednonsampairtuplst
    
def DC3_nonsamradixsort(nonsampairtuplst, iternum):
    """
    Sort the non-sample suffix tuples (t_i, trank) where trank = rank(S_i+1)
    Input:
        nonsampairtuplst   - a list of non-sample suffix tuples [((t_i, rank(S_i+1)), origindex), ...]
        iternum     - the iteration number, which is in {0, 1} since we only have 2 values for radix sort
    Output:
        ordered - the radix-sorted list of non-sample suffix tuples with their original indices
    """

    # stop at iternum = 2 since we only have 2 values to compare
    if iternum > 1:
        return nonsampairtuplst

    # buckets have form {1:[ (1,1), (1,4), ...]
    #                    3:[ (3,2), (3,20), ...}
    buckets = {}
    # process all sample suffixes using radix sort
    for nonsampairtup in nonsampairtuplst:
        nonsampair = nonsampairtup[0]   # extract the (t_i, rank(S_i+1)) tuple
        # need to sort by trank since chars t are equivalent
        if iternum > 0:
            label = nonsampair[1]   # sorting by trank
            if label not in buckets:
                buckets[label] = []
            buckets[label].append(nonsampairtup)
        # first sort by chars t
        else:
            # sort suffixes into buckets
            label = nonsampair[0]   # sorting into buckets by t_i
            # create bucket if it does not yet exist
            if label not in buckets:
                buckets[label] = []
            # buckets will look like [ (1,1), (2,3), ...]]
            buckets[label].append(nonsampairtup)
    # sort bucket keys numerically increasing, 1,2,3,...
    keysortedbuckets = sorted(buckets, cmp=lambda x,y:cmp(x,y)) # sort by numerical value
    # order the suffixes
    ordered = []
    for key in keysortedbuckets:
        # get list of tuples in the bucket
        samples = buckets[key]
        # recurse on buckets with multiple tuples
        if len(buckets[key]) > 1:   # more than one tuple in a bucket
            samples = DC3_nonsamradixsort(samples, iternum+1) # use next indexed num in suffix
        ordered.extend(samples)
    return ordered

def DC3_mergesamples(numlst, sortedB0, sortedC, indexrankmap):
    """
    Merge the sample sets C and B0
    
    B0 form:    [((0, 0), 12), ((1, 5), 6), ((1, 7), 9), ((2, 2), 3), ((3, 1), 0)]
    C=B1UB2 form: [([1, 2, 3], 1), ([2, 3, 4], 2), ([4, 5, 6], 4), ([5, 6, 0], 5)]
    
    S_i corresponds to C=B1UB2
    S_j corresponds to B0
    
    """
    TUPIND = 0
    SUFFIX = 0
    INDEX = 1    
    
    allsorted = []
    suffixarray = {}  # rank:index
    # loop through every suffix
    k = 0
    while k < len(numlst):
        # no more pairs in either list
        if (not sortedC) and (not sortedB0):
            break
        # no more pairs in C=B1UB2
        if not sortedC:
            suffixarray[k] = sortedB0[TUPIND][INDEX]
            allsorted.append(sortedB0.pop(0)) # remove first item from S_j
            k += 1
            continue
        if not sortedB0:
            suffixarray[k] = sortedC[TUPIND][INDEX]
            allsorted.append(sortedC.pop(0))    # remove first item from S_i
            k += 1
            continue
            
        # extract information from sorted lists
        S_i = sortedC[TUPIND]   # ((0, 0), 12)      (suffix, index)
        S_i_suffix = S_i[SUFFIX]
        S_i_index = S_i[INDEX]
        S_j = sortedB0[TUPIND]    # ([1, 2, 3], 1)    (suffix, index)
        S_j_suffix = S_j[SUFFIX]
        S_j_index = S_j[INDEX]
        
        if S_i_index % 3 == 1:
            # handle B1
            # iff ((t_i, rank(S_i+1)) <= (t_j, rank(S_j+1))), then S_i <= S_j
            # compare:
            ipair = (numlst[S_i_index], indexrankmap[S_i_index+1])
            jpair = (numlst[S_j_index], indexrankmap[S_j_index+1])
        elif S_i_index % 3 == 2:
            # handle B2
            # iff (t_i, t_i+1, rank(S_i+2)) <= (t_j, t_j+1, rank(S_j+2)), then S_i <= S_j
            ipair = (numlst[S_i_index], numlst[S_i_index+1], indexrankmap[S_i_index+2])
            jpair = (numlst[S_j_index], numlst[S_j_index+1], indexrankmap[S_j_index+2])
        else:
            raise Exception("Invalid suffix index found in C=B1UB2")

        comparison = compare(ipair, jpair)
        # ipair > jpair
        if comparison > 0:
            suffixarray[k] = S_j_index
            allsorted.append(sortedB0.pop(0)) # remove first item from S_j
        # ipair < jpair
        elif comparison < 0:
            suffixarray[k] = S_i_index
            allsorted.append(sortedC.pop(0))  # remove first item from S_i
        elif comparison == 0:
            raise Exception("Pairs found to be equal!")
        
        k += 1
    
    return allsorted, suffixarray
            
def compare(x, y):
    """
    Returns:
    1 for x > y
    0 for x = y
    -1 for x < y
    """
    # check lengths are the same
    if len(x) != len(y):
        raise Exception("Lengths of merging pairs to compare not equal!")
        
    result = 0
    # check ordering
    for i in range(len(x)):
        if x[i] < y[i]:
            result = -1
            break
        elif x[i] > y[i]:
            result = 1
            break
        
    return result    

def encode_BWT(seq):
    """
    Constructs the Burrows-Wheeler Transform of the input sequence
    Input:
        seq - the input sequence (assumes $ at end of seq)
    Output:
        BWT - the Burrows-Wheeler Transform of seq (which is L in pi_sorted)
    """
    
    suffixarray = DC3(prep_DC3(seq))
    L = []
    for rank in range(len(seq)):
        L.append(seq[suffixarray[rank]-1])
    BWT = ''.join(L)
    return BWT

def decode_BWT(BWT):
    """
    Performs the Inverse Burrows-Wheeler Transform of a BWT.
    This entails constructing a table mapping L to F, and doing
    a forward decoding starting from the lexicographically lowest character in L.
    Input:
        BWT - A string containing the Burrows-Wheeler Transform of a sequence S
    Output:
        iBWT - The original string
    """
    M = get_M(BWT)
    occ = get_occurrences(BWT)

    finalstr = []
    j = 0
    for i in range(len(BWT)):
        finalstr.append(BWT[j])
        j = bd_index(BWT, M, occ, j)
    
    finalstr = ''.join(finalstr)
    finalstr = finalstr[::-1]
    finalstr = finalstr[1:] + finalstr[:1]
    
    return finalstr
        

def get_M(BWT):
    """
    Computes a table for our character starting positions in F in linear time
    Input:
        BWT
    Output:
        M
    """
    F = sorted(BWT)
    
    # find the starting row positions of each character in F
    M = {}
    for index in range(len(F)):
        char = F[index]
        if char not in M:
            M[char] = index
            
    return M
    
def get_occurrences(BWT):
    """
    Returns an occurrence function scoring table in linear time.
    Input:
        BWT - The Burrows-Wheeler Transform of a sequence, or L
    Output:
        occurrences - a table of occurrence function values
            structures as {L-index:{inclusivecharcount}}
    """
    # occurrences is structured as {L-index:{inclusivecharcount}}
    occurrences = {}
    charcounter = Counter()
    # create table for occurrence function
    for charindex in range(len(BWT)):
        char = BWT[charindex]
        charcounter[char] += 1
        occurrences[charindex] = charcounter.copy()
    return occurrences
    
def bd_index(BWT, M, occ, j):
    """
    Compute the backwards decoding index of j
    Input:
        BWT
        M
        occ
        j
    Output:
        bd_index
    """
    
    bdi = occ[j][BWT[j]] + M[BWT[j]] - 1    # account for off-by-one errors
    
    return bdi
    
def split80(string):
    """
    Splits a string into length 80 or less chunks
    """    
    
    chunks = []
    strcopy = string
    while strcopy:
        if len(strcopy) > 80:
            cut = strcopy[:80]
            strcopy = strcopy[80:]
            chunks.append(cut)
        else:
            chunks.append(strcopy)
            strcopy = ""
    return chunks
    
# Executes program when run from command line
def main():
    
    # Read command-line arguments
    parser = ArgumentParser(description="Choose parameters for BWW or iBWT")
    
    parser.add_argument("infile", metavar="in", help="input FASTA file")
    parser.add_argument("outfile", metavar="out", help="output FASTA file")
    
    parser.add_argument("-bwt", "--bwt", action="store_true", help="Compute the Burrows-Wheeler Transform")
    parser.add_argument("-ibwt", "--invbwt", action="store_true", help="Compute the Inverse Burrows-Wheeler Transform")
    
    """
    parser.add_argument("-e","--exclude", dest="exclude", metavar="exfile",
                        help="overhang exclusion list, see README for details")
    parser.add_argument("-p","--percent", nargs=2, dest="topxtuple", 
                        metavar=("percent", "topxfile"), 
                        help="top x percent combos to retain, name of filtered output file")
    """
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="be more verbose")

    args = parser.parse_args()
    
    infile = args.infile
    outfile = args.outfile
    
    template = open(infile)
    newfile = open(outfile, 'w')
    
    # extract sequences from the FASTA file
    for name, seq in process_FASTA(template):
    
        if args.bwt:
            # compute the BWT
            encodedBWT = encode_BWT(seq)
            # output the BWT in the file
            newfile.write(name+"\n")
            strtowrite = encodedBWT
            for line in split80(strtowrite):
                newfile.write(line+"\n")
        elif args.invbwt:
            # compute the iBWT
            decodedBWT = decode_BWT(seq)
            # output the BWT in the file
            newfile.write(name+"\n")
            strtowrite = decodedBWT
            for line in split80(strtowrite):
                newfile.write(line+"\n")
        else:
            raise Exception("Please enter an argument: -bwt or -ibwt")
    template.close()
    newfile.close()
    
if __name__ == "__main__":
    main()

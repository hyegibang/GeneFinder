# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    return complement[nucleotide]


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_complement = ''

    for nucleotide in dna:
        reverse_complement += get_complement(nucleotide)
    return reverse_complement[::-1]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """

    store = ''
    rest_ORF = ''
    stop_codon = ['TAG', 'TAA','TGA']

    for i in range(0,len(dna),3):
        rest_ORF = dna[i:i+3]
        if rest_ORF in stop_codon:
            break
        store = store + rest_ORF

    return store

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    i = 0
    frame2 = []
    start_codon = 'ATG'

    while i < len(dna):
        if start_codon == dna[i:i+3]:
            frame = rest_of_ORF(dna[i:])
            frame2.append(frame) #putting collected data into list
            if frame is None:
                break
            i += len(frame)

        else:
            i += 3

    return frame2

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    """ Difference between append and extend is that extend adds 'objects' to the list
        whereas append adds list within the list -> making the list nested"""

    possible = []
    for i in range(3): # 3 possible frames
        possible.extend(find_all_ORFs_oneframe(dna[i:])) # looks for start codon
    return possible


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    fivetothree = find_all_ORFs(dna)
    threetofive = find_all_ORFs(get_reverse_complement(dna))
    bothstrand = fivetothree + threetofive

    return bothstrand
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORF = find_all_ORFs_both_strands(dna)
    try:
        longest = max(ORF, key=len)
    except:
        return ""
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest = []
    result = []
    for x in range(num_trials):
        result = longest_ORF(shuffle_string(dna))
        longest.append(result)
    longest.sort(key = len)
    return longest[-1]

    passs


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    length_dna = len(dna)
    remainder = length_dna % 3
    if remainder is not 0:
        dna = dna[:-1*remainder]
    aa = ''
    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        new = aa_table[codon]
        aa = aa + new
    return aa


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = len(longest_ORF_noncoding(dna, 1500))
    print('threshold', threshold)
    sequence = []
    for orf in find_all_ORFs_both_strands(dna):
        if len(orf) > threshold:
            sequence.append(coding_strand_to_AA(orf))
    return sequence

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))

from Bio import Restriction
from Bio.Restriction import *
from Bio.Seq import Seq
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
import itertools
from collections import defaultdict
from Bio.Data import CodonTable
import fnmatch
import re
from docx import Document
from docx.enum.text import *
from docx.shared import RGBColor


def reformat_codon_table(table=None):
    """
    Given a codon table, returns a dictionary mapping each codon to its other \
    possible synonymous codons. Each key is at considered a synonymous codon of \
    itself, so it is included in the list.

    Parameters
    ----------
    table : dict
        Dictionary mapping codon to amino acid. If None, then use
        Bio.Data.CodonTable.standard_dna_table.forward_table

    Returns
    -------
    dict
        Dictionary mapping each codon to the set of alternative codons (that does
        not include itself).

    Example
    -------

    >>> d = {'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M'}
    >>> expected = {
    ... 'ATT': ['ATT', 'ATC', 'ATA'],
    ... 'ATC': ['ATT', 'ATC', 'ATA'],
    ... 'ATA': ['ATT', 'ATC', 'ATA'],
    ... 'ATG': ['ATG'],
    ... 'TAA': ['TAA'],
    ... 'TGA': ['TGA'],
    ... 'TAG': ['TAG']}
    >>> reformat_codon_table(d) == expected
    True
    """
    if table is None:
        table = CodonTable.standard_dna_table.forward_table

    # Convert dict of form codon: AA to dict of form AA: [codons]
    synonymous_per_aa = defaultdict(list)
    for k, v in table.items():
        synonymous_per_aa[v].append(k)

    # Convert dict of form AA: [codons] to dict of form codon: [alt codons]
    synonymous = defaultdict(list)
    for v in synonymous_per_aa.values():
        for i in v:
            synonymous[i] = list(v)
    synonymous_dict = dict(synonymous)
    # add stop codons
    synonymous_dict['TAA'] = ['TAA']
    synonymous_dict['TGA'] = ['TGA']
    synonymous_dict['TAG'] = ['TAG']
    return synonymous_dict

def synonymous_seqs(seq, table=None):
    """
    Generator of all possible synonymous versions of the provided sequence \
    using the provided table.

    Parameters
    ----------
    seq : str or Bio.Seq
        Sequence containing codons from which synonymous versions will be
        generated. Must be divisible by 3.

    table : dict
        Optional codon: amino acid dictionary passed to reformat_codon_table;
        see that function for details.

    Yields
    -------
    Generator
        Generator of synonymous sequences (each one a string of same length as
        original).

    Example
    -------

    >>> list(synonymous_seqs('ATCCAATGG'))
    ['ATTCAATGG', 'ATTCAGTGG', 'ATCCAATGG', 'ATCCAGTGG', 'ATACAATGG', 'ATACAGTGG']
    """
    if isinstance(seq, Seq):
        seq = str(seq)

    if len(seq) % 3 != 0:
        raise ValueError(
            'Sequence is length {0}, not divisible by 3'.format(len(seq)))

    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]

    d = reformat_codon_table(table)
    expanded = [d[c] for c in codons]

    for i in itertools.product(*expanded):
        yield ''.join(i)


def get_endogenous_sequence_and_stats(
    full_seq,
    gRNA_start,
    gRNA_end,
    PAM_start,
    AA_change_site,
    which_PAM
):
    """
    Given a full sequence and the locations of the gRNA, PAM, and AA change \
    site, returns a dictionary of the isolated 273 bp sequence that extends 155 \
    nt upstream and 120 nt downstream (or as close to that as the sequence \
    extends) of the Cas9 cutsite and the adjusted gRNA, PAM, and AA change site \
    values that correspond to the isolated endogenous sequence.

    Parameters
    ----------
    full_seq: string
        The full DNA locus sequence, must extend at least 85 nt upstream and
        50 nt downstream of the Cas9 cutsite.

     gRNA_start: int
        the index location in the sequence of the beginning of the
        region the CRISPR/Cas9 guide RNA will target.

    gRNA_end: int
        the index location in the sequence of the end of the region
        the CRISPR/Cas9 guide RNA will target.

    PAM_start: int
        The index location in the sequence of the beginning of the
        PAM site. For now, the PAM region has to be either 'CCN' or 'NGG'.

    AA_change_site: int
        the index location in the sequence of the beginning of
        the codon that is to be nonsynonymously changed in the CRISPR/Cas9
        process

    which_pam: string
        The sequence of the PAM site (for now 'CCN', 'NGG', 'TTTV', and 'BAAA'
        are supported)

    Returns
    -------
    dict
        Dictionary of the isolated endogenous sequence and the adjusted gRNA, PAM,
        and AA_change_site locations.
    """
      
    # create dictionary with PAM startsites as keys and the locations of the
    # CRISPR cutsites relative to the PAM startsites as values
    pam_dict = {
        'CCN': 6,
        'NGG': -3,
        'TTTV': 23,
        'BAAA': -23
    }

    full_seq = full_seq.replace(" ", "").upper()
    Cas9_cutsite = PAM_start + pam_dict[which_PAM]
    seq_stats = {
        'gRNA_start': gRNA_start,
        'gRNA_end': gRNA_end,
        'PAM_start': PAM_start,
        'AA_change_site': AA_change_site
    }

    #cut it down to 15 nt on either side of cutsite + 70 nt upstream and 35
    #nt downstream + 70 nt spacers on either side.
    new_zero = 0
    new_end = len(full_seq)
    if which_PAM == 'CCN' or which_PAM == 'TTTV':
        if Cas9_cutsite > 155:
            new_zero = Cas9_cutsite - 155
        if new_end - Cas9_cutsite > 120:
            new_end = Cas9_cutsite + 120

    if which_PAM == 'NGG' or which_PAM == 'BAAA':
        if Cas9_cutsite > 120:
            new_zero = Cas9_cutsite - 120
        if new_end - Cas9_cutsite > 155:
            new_end = Cas9_cutsite + 155
    focus_seq = full_seq[new_zero: new_end]

    #adjust all of the seq_stats to correspond focus_seq rather than full_seq
    seq_stats.update((k,v - new_zero) for (k,v) in seq_stats.items())
    #check the frame of the new sequence based the AA change site, which
    #has to be at the beginning of a codon.
    frame = seq_stats['AA_change_site'] % 3
    #adjust the sequence to be in frame and divisible by 3
    focus_seq = focus_seq[frame: 273 + frame]
    #adjust the important locations to match the in frame focus_seq
    seq_stats.update((k, v - frame) for (k, v) in seq_stats.items())
    #add focus_seq to the dictionary before returning it
    seq_stats['sequence'] = focus_seq
    seq_stats['which_PAM'] = which_PAM
    return seq_stats


def generate_exon_seqs(start, end, AA_change_seq):
    """
    Given a sequence, and start and end coordinates for the cutsite region \
    within that sequence, returns a generator of all possible synonymous \
    cutsite regions.

    Parameters
    ----------
    start: int
        The starting location of the cutsite region within the sequence

    end: int
        The ending location of the cutsite region within the sequence

    AA_change_sequence: string
        The full candidate sequence with the desired amino acid change

    Yields
    -------
    generator
        A generator of sequences synonymous to the region within AA_change seq
        specified by the start and end coordinates.
    """
    start_frame = start % 3
    end_frame = end % 3
    in_frame_start = start - start_frame
    in_frame_end = end
    if end_frame != 0:
        in_frame_end = in_frame_end + 3 - end_frame
    seqs = synonymous_seqs(AA_change_seq[in_frame_start:in_frame_end])
    for seq in seqs:
        if seq[:start_frame] == AA_change_seq[in_frame_start: start]:
            if seq[len(seq) - 3 + end_frame:] == AA_change_seq[end: in_frame_end]:
                yield seq[start_frame: len(seq) -3 + end_frame]
            elif end_frame == 0:
                yield seq[start_frame:]

def get_restriction_sites(sequence, Renzymes=None):
    """
    Given a sequence, returns a dictionary of restriction enzymes mapping \
    Restriction enzymes that would cut the sequence to locations in the \
    sequence that they would cut.

    Parameters
    ----------
    sequence: string
        A DNA nucleotide sequence. Must have all capital letters.

    Renzymes: list
        Optional enzyme list: a list of restriction enzymes to search for in the sequence.
        The default list is the Bio.Restriction package's 'AllEnzymes'
        module.

    Returns
    -------
    dict
        Dictionary of restriction enzymes with recognition sites in sequence.
        The keys are RE names and the values are their cutsite locations in
        the sequence. NOTE: the cutsites use the convention that the first
        nucleotide in the sequence is 1 not 0.

    Example
    -------

    >>> result =  get_restriction_sites('ATGCGA')
    >>> expected = {RlaI: [2], PsuGI: [2], BmeDI: [3, 5]}
    >>> assert result == expected
    """
    seq = Seq(sequence)
    if Renzymes is None:
        rb = AllEnzymes
    else:
        rb = RestrictionBatch(list(Renzymes))
    ana = Analysis(rb, seq)
    return ana.with_sites()

def get_Rsite_dict(sequence, Renzymes):
    """
    Given a sequence and either a list or a dictionary of restriction enzymes, \
    returns a dictionary mapping the restriction enzymes to coordinates of \
    their recognition sites in the sequence.

    Parameters
    ----------
    sequence: string
        A DNA nucleotide sequence. Must be all capital letters.

    Renzymes: list or dict
        Either a list of restriction enzyme or a dictionary with restriction
        enzymes as keys. The restriction enzymes must be in the form of
        Bio.Restriction restriction enzyme objects.

    Returns
    -------
    dict
        Dictionary of enzymes with recognition sites in 'sequence', with the
        enzyme names as keys and tuples of the start and end sites of their
        recognition sites in the sequence as values.
        NOTE: the endsites correspond to the first nucleotide in the sequence after
        the recognition site (e.g. recognition site  = sequence[start:end] NOT
        sequence[start:(end + 1)])

    Example
    -------
    >>> result = get_Rsite_dict('ATGCAC', {HpyCH4V: [4], CviRI: [4], RlaI: [2, 3], BmeDI: [3], Gba708II: [1]})
    >>> expected = {
    ... HpyCH4V: [(1, 5)],
    ... CviRI: [(1, 5)],
    ... RlaI: [(1, 4), (2, 5)],
    ... BmeDI: [(2, 3), (3, 4), (5, 6)],
    ... Gba708II: [(0, 5)]}
    >>> assert result == expected, result
    """

    #Use the Renzyme dict or list to make a new dictionary with the recognition site coordinates
    Rcoords = {}
    for Renzyme in Renzymes:
        Rcoords[Renzyme] = [tuple([n.start(), n.start() + len(Renzyme.site)]) for n in Renzyme.compsite.finditer(sequence)]
    return Rcoords

def map_recognition_sites(sequence, Renzyme, num):
    """
    Given a sequence, a restriction enzyme that would cut the sequence, and \
    a nonzero integer ('num'), finds the locations of all the recognition sites \
    for that enzyme in the sequence and maps their locations. 

    Parameters
    ----------
    sequence: string
        A DNA nucleotide sequence. All of the letters must be capitalized.

    Renzyme: Bio.Restriction enzyme object
        A restriction enzyme in the form of a Bio.Restriction enzyme object.

    Returns
    -------
    list
        List of 0's the length of sequence, with the 0's replaced with 'num' at
        any point where the corresponding location in the sequence is part of
        a recognition site in sequence.


    Example
    -------
    #TaqI recognition site: 'TCGA'
    >>> map_recognition_sites('AGCTCGATG', TaqI, 1)
    [0, 0, 0, 1, 1, 1, 1, 0, 0]

    """
    Rsite_dict = get_Rsite_dict(sequence, [Renzyme])
    num_list = [0 for x in range(0, len(sequence))]
    for start, stop in Rsite_dict[Renzyme]:
        for i in range(start, stop):
            num_list[i] += num
    return num_list


def get_substitutions_Renzyme_depends_on(Renzyme, endogenous_seq, candidate_seq, sub_dep_sites):
    """
    Given a restriction enzyme, an endogenous sequence, and a \
    candidate sequence, returns a list of tuples where each tuple corresponds to \
    a recognition site for the enzyme in the sequence and contains the locations of \
    synonymous substitutions in candidate_seq that that recognition site \
    depends on.

    Parameters
    ----------
    Renzyme: Bio.Restriction enzyme object
        A restriction enzyme in the form of a Bio.Restriction enzyme object.

    endogenous_seq: string
        A DNA nucleotide sequence. All letters must be capitalized.

    candidate_seq: string
        A candidate sequence based on the endogenous sequence. In other words, a sequence that's similar to the endogenous sequence, but it has one nonsynonymous substitution and a series of synonymous substitutions in the cutsite region.

    sub_dep_sites: list
        A list of tuples indicating start and end sites for substitution dependent recognition sites for Renzyme in candidate seq.

    Returns
    ------
    list
        List of tuples where each tuple corresponds to a novel recognition site for
        Renzyme (AKA a recognition site for Renzyme that is in the candidate
        sequence but not the endogenous sequence due to substitutions to the
        canidate sequence). The numbers in each tuple correspond to index locations
        in the candidate sequence of substutions that are required to form its
        corresponding recognition site.

    Example
    -------
    >>> get_substitutions_Renzyme_depends_on(Hpy188I, 'TTCTGCAGCCCTGGATGGGCTCTGCTGCTA', 'TTCTTCAGAGTTGGATGGGTTCTGCTGCTA',[(4,9)])
    [(4, 8)]
    """
    nucleotide_code = {'N': ['G', 'A', 'T', 'C'],
                       'W': ['A', 'T'],
                       'S': ['C', 'G'],
                       'M': ['A', 'C'],
                       'K': ['G', 'T'],
                       'R': ['A', 'G'],
                       'Y': ['C', 'T'],
                       'B': ['C', 'G', 'T'],
                       'D': ['A', 'G', 'T'],
                       'H': ['A', 'C', 'T'],
                       'V': ['A', 'C', 'G']
                      }

    # E.g., 
    # nuc_delta = [0000000100010000]
    # where 1 indicates difference from endogenous
    nuc_delta = [0 for x in range(0, len(candidate_seq))]
    for i in range(0, len(nuc_delta)):
        if endogenous_seq[i] != candidate_seq[i]:
            nuc_delta[i] += 1


    # Each RE can have multiple sites; we'll be checking them all here and
    # putting the info for each site in `all_sub_sites`
    all_sub_sites= []
    sites = sub_dep_sites
    for start, stop in sites:
        sub_sites = []
        for i in range(start, stop):
            if nuc_delta[i]==1:

                # If an RE has an ambiguous base, need to handle that here.
                # For example, if a substitution was C to T, but the RE has
                # a "Y" in it's site, then it shouldn't be considered
                # a substitution from the point of view of the RE
                recognition_bp = Renzyme.site[i - start]
                if recognition_bp in nucleotide_code:
                    if endogenous_seq[i] not in nucleotide_code[recognition_bp]:
                        sub_sites.append(i)
                else:
                    sub_sites.append(i)

        all_sub_sites.append(tuple(sub_sites))
    return all_sub_sites

def get_codon_frequency_dict(codon_usage_tsv):
    """
    Takes in a codon usage frequency tsv and converts it into a dictionary

    Parameters
    ----------
    codon_usage_tsv: a tsv file of codons and their corresponding frequencies in the format [codon]>[frequency] with one codon per line

    Returns
    -------
    dict
        A dictionary of codons and the frequency at which they are used in

    """
    codon_frequency_dict = {}
    for line in open(codon_usage_tsv):
        if not line.startswith('#'):
            linelist = line.strip().split('\t')
            codon = linelist[0].replace('U', 'T')
            codon_frequency_dict[codon] = float(linelist[1])
    return codon_frequency_dict

def get_filter_function_dict(filter_list, candidate):
    """
    Given a list of filters, return a dictionary with the filter names in the \
    filter list as keys and their corresponding filter functions as values.
    
    Parameters
    ----------
    filter_list: list
        A list of filters to get corresponding functions for
    candidate: candidate object
        The candidate object to pull the methods from
    
    Returns
    -------
    dict
        A dictionary with the filter names in the
        filter list as keys and their corresponding filter functions as values.
    """

    full_filter_dict = {
        'no_sandwiched_cutsites': candidate.no_sandwiched_cutsites,
        'no_endogenous_cutsites_within_Nnt_of_novel_cutsite': candidate.no_endogenous_cutsites_within_Nnt_of_novel_cutsite,
        'novel_site_depends_on_at_least_two_substitutions': candidate.novel_site_depends_on_at_least_two_substitutions,
        'num_recognition_sites_equals_num_cutsites': candidate.num_recognition_sites_equals_num_cutsites,
        'no_unneeded_substitutions': candidate.no_unneeded_substitutions,
        'at_least_one_has_endogenous_site': candidate.at_least_one_has_endogenous_site,
        'both_have_endogenous_sites': candidate.both_have_endogenous_sites,
        'compatible_polarity': candidate.compatible_polarity,
        'both_Renzymes_depend_on_at_least_one_unique_substitution': candidate.both_Renzymes_depend_on_at_least_one_unique_substitution
    }
    return({filter_name: filter_function for (filter_name, filter_function)
            in full_filter_dict.items() if filter_name in filter_list}) 

def get_codon_frequency_differences(seqA, seqB, codon_dict):
    """
    Given two sequences and a dictionary of usage frequences for all codons, \
    returns the total difference in codon frequency between the two \
    sequences.

    Parameters
    ----------
    seqA: string
        An in frame DNA nucleotide sequence

    seqB: string
        An in frame DNA nucleotide sequence of the same length as seqA

    codon_dict: dict
        A dictionary of codons and the frequency at which they are used

    Returns
    -------
    float
        A float describing the total difference in usage frequency across codons
        between the two input sequences.
    """

    total_frequency_diff = 0
    # iterate over the starting coordinates for each codon in the sequences and
    # establish the codon frequency differences based on codon_dict
    for i in range(0, (int(len(seqA)/3)) - 1):
        seqA_codon = seqA[i * 3: (i * 3) + 3]
        seqB_codon = seqB[i * 3: (i * 3) + 3]
        frequency_diff = abs(codon_dict[seqA_codon] - codon_dict[seqB_codon])
        total_frequency_diff += frequency_diff
    return total_frequency_diff

def add_sequence_to_doc(document, sequence, gRNA_start, gRNA_end, PAM_start, AA_change_site, Renzyme_list, endogenous, which_pam, paired = True):
    """A helper function for make_doc(). The make_doc function calls \
    it to one of the candidate sequences to the document with the Renzyme \
    recognition sites, synonymous substitutions, PAM and gRNA and \
    AA_change_site coordinates marked according to the key.

    Parameters
    ----------
    document: docx.document
        The document that the make_doc() function is creating when it calls
        this function.

    sequence: string
        A candidate sequence to be evalutated as a possible repair oligo for
        the endogenous sequence

    gRNA_start: int
        The location of the beginning of the gRNA target region in the
        sequence.

    gRNA_end: int
        The location of the end of the gRNA target region in the sequence.

    PAM_start: int
        The location of the beginning of the PAM sequence in the sequence.

    AA_change_site: int
        The location of the beginning of the codon to be nonsynonymously
        changed in the CRISPR/Cas9 sequence.

    Renzyme_list: list
        An enzyme or pair of enzymes that passed all of the filters for
        sequence

    endogenous: int
        The endogenous sequence that the candidate sequence is based on
    
    which_pam: string
        The sequence of the PAM site (for now 'CCN', 'NGG', 'TTTV', and 'BAAA'
        are supported)

    paired: bool
        A boolean value indicating whether or not to run the function on single
        enzymes or enzyme pairs. Default is enzyme pairs (AKA True)
    """
    Rsite_dict = get_Rsite_dict(sequence, Renzyme_list)
    if paired:
        RenzymeA_sites = map_recognition_sites(sequence, Renzyme_list[0], 1)
        RenzymeB_sites = map_recognition_sites(sequence, Renzyme_list[1], 7)
        recognition_sites = [x + y for x, y in zip(RenzymeA_sites, RenzymeB_sites)]
        document.add_paragraph(str(Renzyme_list[0]) + ' (' + Renzyme_list[0].elucidate() + ') and ' + str(Renzyme_list[1]) + ' (' + Renzyme_list[1].elucidate() + ')')
    else:
        recognition_sites = map_recognition_sites(sequence, Renzyme_list[0], 1)
        document.add_paragraph(str(Renzyme_list[0]) + '(' + Renzyme_list[0].elucidate() + ')')
# start a new line for the sequence
    seq = document.add_paragraph()

# add onto the sequence character by character
    for i in range(0, len(sequence)):
        run_i = seq.add_run(sequence[i])
        if i in range(gRNA_start, gRNA_end):
            run_i.font.highlight_color = WD_COLOR_INDEX.GRAY_25
        if i in range(PAM_start, PAM_start + len(which_pam)):
            run_i.font.highlight_color = WD_COLOR_INDEX.YELLOW
        if i in range(AA_change_site, AA_change_site + 3):
            run_i.bold = True
            run_i.italic = True
        if recognition_sites[i]  in [1, 2, 3]:
            run_i.font.color.rgb = RGBColor(255, 0, 0)
        if recognition_sites[i] in [7, 14, 21] :
            run_i.font.color.rgb = RGBColor(0, 0, 255)
        if recognition_sites[i] in [8, 9, 10, 15, 16, 17, 22, 23, 24]:
            run_i.font.color.rgb = RGBColor(127, 0, 255)
        if endogenous is not None:
            if sequence[i] != endogenous[i]:
                run_i.underline = WD_UNDERLINE.THICK

def make_doc(endogenous_sequence, gRNA_start, gRNA_end, PAM_start, AA_change_site, candidate_df, sample, which_pam, paired=True):
    """
    Makes a word document to visualize restriction sites \
    for each pair of enzymes (or single enzyme if paired != True) for each \
    candidate in candidate_dict

    Parameters
    ----------
    endogenous_sequence: string
        A DNA nucleotide sequence. Must be all capital letters.

    gRNA_start: int
        The location of the beginning of the gRNA target region in the sequence.

    gRNA_end: int
        The location of the end of the gRNA target region in the sequence.

    PAM_start:int
        The location of the beginning of the PAM sequence in the sequence.

    AA_change_site: int
        The location of the beginning of the codon to be nonsynonymously
        changed in the CRISPR/Cas9 sequence.

    candidate_df: dataframe
        Pandas dataframe of candidate with either single enzymes or enzyme pairs that passed all filters. If single enzymes, it should have the 
        following columns: {sequence, Renzyme, average distance to cas9 cut, codon frequency change, cutsite distance outside Rsite}. If paired, it
        should have the following columns: {sequence, Renzyme pair, result, RenzymeA average distance to cas9 cut, RenzymeB average distance to cas9
        cut, overall average distance to cas9 cut, codon_frequency_change}

    sample: string
        A unique identifier for this locus

    which_pam: string
        The sequence of the PAM site (for now 'CCN', 'NGG', 'TTTV', and 'BAAA'
        are supported)

    paired: bool
        A boolean value indicating whether or not to run the function on single
        enzymes or enzyme pairs. Default is enzyme pairs (AKA True)
    """
    
    doc = Document()
    key = doc.add_paragraph().add_run('Key:')
    key.bold = True
    gRNA_key = doc.add_paragraph().add_run('gRNA')
    gRNA_key.font.highlight_color = WD_COLOR_INDEX.GRAY_25
    PAM_key = doc.add_paragraph().add_run('PAM site')
    PAM_key.font.highlight_color = WD_COLOR_INDEX.YELLOW
    AA_change_key = doc.add_paragraph().add_run('Amino acid change')
    AA_change_key.bold = True
    AA_change_key.italic = True
    if paired:
        RenzymeA_key = doc.add_paragraph().add_run('Renzyme A recognition site')
        RenzymeA_key.font.color.rgb = RGBColor(255, 0, 0)
        RenzymeB_key = doc.add_paragraph().add_run('Renzyme B recognition site')
        RenzymeB_key.font.color.rgb = RGBColor(0, 0, 255)
        overlapping_sites_key = doc.add_paragraph().add_run('Overlapping Renzyme A and B recognition sites')
        overlapping_sites_key.font.color.rgb = RGBColor(127, 0, 255)
    else:
        Renzyme_key = doc.add_paragraph().add_run('Renzyme recognition site')
        Renzyme_key.font.color.rgb = RGBColor(255,0,0)
    cutsite_key1 = doc.add_paragraph().add_run("^ = cut on (+) strand")
    cutsite_key2 = doc.add_paragraph().add_run("_ = cut on (-) strand")
    substitution_key = doc.add_paragraph().add_run('base pair substitution')
    substitution_key.underline = WD_UNDERLINE.THICK
    header = doc.add_paragraph().add_run('Endogenous sequence:')
    header.bold = True
    Endogenous = doc.add_paragraph()
    for i in range(0, len(endogenous_sequence)):
                run_i = Endogenous.add_run(endogenous_sequence[i])
                if i in range(gRNA_start, gRNA_end):
                    run_i.font.highlight_color = WD_COLOR_INDEX.GRAY_25
                if i in range(PAM_start, PAM_start + len(which_pam)):
                    run_i.font.highlight_color = WD_COLOR_INDEX.YELLOW
    for index, row in candidate_df.iterrows():
        header = doc.add_paragraph().add_run('Option ' + str(index) + ':')
        header.bold = True
        if paired:
            Renzyme_list = list(RestrictionBatch(row['Renzyme pair'].split(' and ')))
        else:
            Renzyme_list = list(RestrictionBatch([row['Renzyme'].strip()]))
        add_sequence_to_doc(doc, row['seq'], gRNA_start, gRNA_end, PAM_start, AA_change_site, Renzyme_list, endogenous_sequence, which_pam,  paired)
        doc.add_paragraph()
    if paired:
        doc.save(sample + '_Renzyme_pairs.docx')
    else:
        doc.save(sample + '_single_Renzymes.docx')


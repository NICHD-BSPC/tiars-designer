__version__='0.1.5'
from Bio.Restriction import *
import Levenshtein
from . import utils
import re

class Designer():
    """The parent class for Endogenous and Candidate, which stores functions and
    variables common to both"""

    def __init__(self, sequence, cutsite, Renzyme_list=None):

        """
        Parameters
        ----------

        sequence : string
            An in frame sequence of DNA nucleotides of max length 273
            nucleotides. It should have a PAM sequence that makes it suitable
            for CRISPR/Cas9 activity.

        cutsite : int
            The index location of the CRISPR/Cas9 cutsite (3 nt after the PAM
            if it's CCN, 3 before if it's NGG).
        Renzyme_list : list
            Optional: provide a list of Renzymes to use. Default is
            Bio.Restriction AllEnzymes.

        """

        self.sequence = sequence
        """string : The input sequence"""
        self.cutsite = cutsite
        """int : The input cutsite"""
        self.cutsite_region = sequence[(cutsite - 15):(cutsite + 15)]
        """string : A subset of sequence that extends 15 bp to either side of cutsite"""
        self.cutsite_frame = self.cutsite % 3
        """int: A number describing the frame of the cutsite region. '0' means it's
        in frame. '1' means the second bp of the cutsite region is the
        first bp of the first codon, and '2' means the third bp of the
        cutsite region is the first bp of the first codon.
        """
        self.Renzyme_list = Renzyme_list
        
        self.Renzymes = utils.get_restriction_sites(sequence, Renzyme_list)
        """dict: A dictionary of restriction enzymes with a recognition sites  in
        self.sequence, with re names as keys and their cutsite as the
        values.
        """

    def get_distance_from(self, pos=0):
        """
        Takes all of the restriction enzymes in self.Renzymes and returns
        a dictionary with REnames as their values and their distances from
        `pos` as the values.

        Parameters
        ----------

        pos : int
            an integer in range(0, len(self.sequence))

        Returns
        -------
        
        dict
            A dictionary of Restriction enzymes that would cut the sequence, with
            the enzymes as keys and the distances between their cutsites and pos as
            values. A positive distance means that the REcutsite is later in
            the string than `pos`, while a negative distance means that the
            REcutsite appears before `pos` in the string.

        """
        REdistances = {}
        for enzyme in self.Renzymes:
            # subtract 1 from the cutsites because they consider the first
            # bp in the sequence to be 1 rather than 0

            REdistances[enzyme] = [
                i - 1 - pos for i in self.Renzymes[enzyme]
            ]

        return REdistances

    def get_splicing_hexamers(self, sequence):
        """
        This function searches the input sequence for splicing hexamers (6 nt
        sequences that enhance premRNA splicing, AKA exonic splicing enhancers)
        from the list provided by RESCUE-ESE
        (http://genes.mit.edu/burgelab/rescue-ese/zebrafish). This list has
        been copied into 'splicing_hexamers.txt' with duplicate hexamers
        removed and returns a dictionary of splicing hexamers found in the
        input sequence and their locations in the sequence.

        Parameters
        ----------
        sequence : string
            A DNA nucleotide sequence. Must be all capital letters.

        Returns
        -------
        dict
            Dictionary of hexamers in the sequence with the hexamer sequences as
            keys and their locations in the input sequence as values
        
        """
        hexamers = [line[:6] for line in open('splicing_hexamers.txt')]
        seq_hexamers = [hexamer for hexamer in hexamers if hexamer in sequence]
        seq_hexamer_sites = {}
        for hexamer in seq_hexamers:
            seq_hexamer_sites[hexamer] = (
                [x.start() for x in re.finditer(hexamer, sequence)]
            )
        return seq_hexamer_sites


class Endogenous(Designer):
    """An Endogenous object is a subclass of the designer class, and an
    Endogenous object is the endogenous sequence for a CRISPR/Cas9 experiment
    before any edits are made to the sequence
    """

    def __init__(
        self,
        sequence,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_codon,
        Renzyme_list=None
    ):
        """
        Parameters
        ----------

        sequence : string
            The endogenous/wt dna sequence without any synonymous changes
            and without the desired CRISPR facilitated amino acid change. Max
            length is 273 bp, must be in frame.

        gRNA_start : int
            the index location in the sequence of the beginning of the
            region the CRISPR/Cas9 guide RNA will target.

        gRNA_end : int
            the index location in the sequence of the end of the region
            the CRISPR/Cas9 guide RNA will target.

        PAM_start : int
            The index location in the sequence of the beginning of the
            PAM site. For now, the PAM region has to be either 'CCN' or 'NGG'.
        
        which_PAM : string
            The sequence of the PAM site (e.g. 'NGG' or 'CCN' for Cas9, or 'TTTN' or 'BAAA' for Cfp1)

        AA_change_site : int
            the index location in the sequence of the beginning of
            the codon that is to be nonsynonymously changed in the CRISPR/Cas9
            process

        new_codon : string
            the codon to replace the existing codon at the
            AA_change_site
        """

        # create dictionary with PAM startsites as keys and the locations of
        # the CRISPR cutsites relative to the PAM startsites as values
        PAM_dict = {
            "CCN": 6,
            "NGG": -3,
            "TTTV": 23,
            "BAAA": -23
        }
        self.gRNA_location = (gRNA_start, gRNA_end)
        """tuple : A tuple of the start and end coordinates of the region in sequence
        to be targeted by the guide RNA."""
        
        self.PAM_start = PAM_start
        """int : The index location of the beginning of the PAM site in sequence.
        For now, the PAM site has to be either 'CCN' or 'NGG'."""
        
        self.which_PAM = which_PAM
        """string : The sequence of the PAM site (e.g. 'NGG' or 'CCN' for Cas9, or 'TTTN' or 'BAAA' for Cfp1)"""
        
        self.cas9_cut_from_PAM = PAM_dict[which_PAM]
        """int : The location of the CRISPR/Cas9 cutsite in the sequence relative to
        the PAM_start site. For 'CCN' it will be 6, and for 'NGG' it will
        be -3. Since the value is different for each PAM, it is also useful
        for indicating whether the PAM sequence is 'CCN' or 'NGG'."""
        super().__init__(sequence, PAM_start + self.cas9_cut_from_PAM, Renzyme_list)
        
        self.AA_change_site = AA_change_site
        """int : The index location in the sequence of the codon that is to be
        nonsynonymously changed in the CRISPR/Cas9 process."""

        self.AA_change_sequence = (
            self.sequence[:self.AA_change_site]
            + str(new_codon)
            + self.sequence[(self.AA_change_site + 3):]
        )
        """string : A sequence identical to the endogenous sequence (self.sequence),
        but it has the new_codon substitution at the AA site.  In other
        words, it's the endogenous sequence with the desired CRISPR/Cas9
        nonsynonymous substitution."""

        self.PAM_can_change = self.PAM_can_change()
        """ bool : A boolean value indicating whether or not the PAM site can be
        interrupted via only synonymous substitutions. This means that one
        or both of the non 'N' basepair's in the PAM site can be changed to
        anything else without changing the resulting amino acid sequence.
        This value will also be True if the desired amino acid change
        interrupts the PAM."""


    def PAM_can_change(self):
        """
        This function will return a boolean value indicating whether or not
        the PAM site can be interrupted via synonymous substitutions or the
        desired nonsynonymous change.  For the PAM to be interrupted, one or
        both of the non 'N' basepair's in the PAM site can be changed to any
        other bp.

        Returns
        -------
        bool
            False if the PAM cannot be interrupted, True if it can
        """
         
        PAM_start = self.PAM_start
        PAM_frame = PAM_start % 3
        AA_change_seq = self.AA_change_sequence
        # if PAM site is 'CCN'
        if self.which_PAM == 'CCN':
            if PAM_frame == 0:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start:PAM_start + 3]
                    )
                )
                for seq in syn_seqs:
                    if seq[0:2] != 'CC':
                        return True
                return False
            if PAM_frame == 1:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start - 1: PAM_start + 2]
                    )
                )
                for seq in syn_seqs:
                    if seq[1:] != 'CC':
                        return True
                return False
            else:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start - 2: PAM_start + 4]
                    )
                )
                for seq in syn_seqs:
                    if seq[2:4] != 'CC':
                        return True
                return False
        # if PAM is 'NGG'
        if self.which_PAM == 'NGG':
            if PAM_frame == 0:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start:PAM_start + 3]
                    )
                )
                for seq in syn_seqs:
                    if seq[1:] != 'GG':
                        return True
                return False
            if PAM_frame == 1:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start - 1: PAM_start + 5]
                    )
                )
                for seq in syn_seqs:
                    if seq[2:4] != 'GG':
                        return True
                return False
            else:
                syn_seqs = list(utils.synonymous_seqs(
                    AA_change_seq[PAM_start + 1: PAM_start + 4])
                )
                for seq in syn_seqs:
                    if seq[0:2] != 'GG':
                        return True
                return False
        if self.which_PAM == 'TTTV' or self.which_PAM == 'BAAA':
            if PAM_frame == 0:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start:PAM_start + 6]
                    )
                )
                for seq in syn_seqs:
                    if self.which_PAM == 'TTTV':
                        if seq[0:3] != 'TTT' or seq[3] == 'T':
                            return True
                    else:
                        if seq[1:4] != 'AAA' or seq[0] == 'A':
                            return True
                return False
            if PAM_frame == 1:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start - 1: PAM_start + 5]
                    )
                )
                for seq in syn_seqs:
                    if self.which_PAM == 'TTTV':
                        if seq[1:4] != 'TTT' or seq[4] == 'T':
                            return True
                    else:
                        if seq[2:5] != 'AAA' or seq[1] == 'A':
                            return True
                return False
            else:
                syn_seqs = list(
                    utils.synonymous_seqs(
                        AA_change_seq[PAM_start -2: PAM_start + 4]
                    )
                )
                for seq in syn_seqs:
                    if self.which_PAM == 'TTTV': 
                        if seq[2:5] != 'TTT' or seq[5] == 'T':
                            return True
                    else:
                        if seq[3:] != 'AAA' or seq[2] == 'A':
                            return True
                return False



    def get_filtered_candidates(self, locusname, single_enzyme_filters, min_Rsite_distance, get_filter_stats, codon_usage_tsv):
        """
        Builds a generator of all potential repair oligo candidates with the
        desired nonsynonymous amino acid change and at least one synonymous
        substitution in self.cutsite_region. Before being run through the
        generator, however, each candidate-Renzyme combination is run through
        all of the single Renzyme filters, and only the combinations that pass
        all of them go through.

        Parameters
        ----------
        locusname : string
            A user created name for the locus for writing the filter stats to file
        single_enzyme_filters: list
            A list of filters to run for this locus (will run all filters if not specified)
        min_Rsite_distance : int
            An integer indicating the minimum allowable distance between any
            novel and endogenous site for a given restriction enzyme
        get_filter_stats : bool
            A boolean value indicating whether or not the function should
            generate filter stats in the form of a percentage of passing
            candidates for each filter.
        codon_usage_tsv : tsv
            A tsv mapping codons to their usage frequency for the organism
            used.
        
        Yields
        ------
        generator
            Generator of all the possible candidate sequence/Renzyme combinations
            that passed all of the single Renzyme filters. A candidate sequence is
            very similar to the endogenous sequence, but it has the amino acid
            change (could be any codon that corresponds to that amino acid) and
            anywhere from one to the maximum possible number of synonymous
            substitutions. If 'get_filter_stats', it will also write out a file
            called '[locusname]_filter_stats.tsv' that includes a percentage of
            candidates that passed for each filter, along with the total number of
            candidates.
        """
        # base the candidates on the endogenous sequence with the desired amino
        # acid change
        codon_freq_dict = utils.get_codon_frequency_dict(codon_usage_tsv)
        AA_change_seq = self.AA_change_sequence
        cutsite = self.cutsite
        # Generate a list of synonymous cutsite regions
        alt_cutsite_regions = utils.generate_exon_seqs(
            cutsite - 15,
            cutsite + 15,
            AA_change_seq
        )


        if get_filter_stats:
            filter_stats = {filter_name: 0 for filter_name in single_enzyme_filters}
        total_cands = 0
        total_passing_cands = 0
        for alt_cutsite_region in alt_cutsite_regions:
            # make the full candidate sequence, which is the endogenous
            # sequence with the cutsite region replaced with
            # alt_cutsite_region. If the cutsite region was extended,
            # restrict it back to 30 bp so no substitutions outsite the 30
            # bp cutsite_region are included.
            candidate_full_seq = (
                AA_change_seq[:cutsite-15]
                + alt_cutsite_region
                + AA_change_seq[cutsite + 15:]
            )
            # make the sequence into a Candidate object and add it to the
            # list!
            candidate = Candidate(self, candidate_full_seq)
            if candidate.PAM_is_disrupted() or candidate.endogenous.which_PAM in ['TTTV', 'BAAA']:

                filter_functions = utils.get_filter_function_dict(single_enzyme_filters, candidate)
                for Renzyme in candidate.sub_dep_Renzymes:
                                       
                    # if getting filter stats, add each filter's result to the tally
                    filter_results = []
                    for Filter in single_enzyme_filters:
                        if Filter == 'no_endogenous_cutsites_within_Nnt_of_endogenous_cutsite':
                            filter_result = filter_functions[Filter](Renzyme, min_Rsite_distance)
                        else:
                            filter_result = filter_functions[Filter](Renzyme)
                        if get_filter_stats:
                            filter_stats[Filter] += filter_result
                        if Filter != 'no_unneeded_substitutions':
                            filter_results.append(filter_result)
                    total_cands += 1
                    # if a candidate passed all the filters, create a candidate
                    # object for it and yield it
                    if all(filter_results):
                        passed_cand_dict = {
                            'seq': candidate.sequence,
                            'Renzyme': str(Renzyme),
                            'codon frequency change': utils.get_codon_frequency_differences(candidate.endogenous.sequence, candidate.sequence, codon_freq_dict),
                            'average distance to Cas cut': candidate.get_average_distance_between_RE_cutsites_and_Cas9_cutsite(Renzyme),
                                                        'cutsite distance outside Rsite': len(Renzyme.elucidate()) - 2 - len(Renzyme.site),
                            'distance to Amino Acid change': candidate.get_distance_from_nonsynonymous_codon_to_Renzyme_cutsite(Renzyme)
                        }
                        if 'no_unneeded_substitutions' in single_enzyme_filters:
                            passed_cand_dict['no unneeded substitutions']= candidate.no_unneeded_substitutions(Renzyme)
                            if filter_functions['no_unneeded_substitutions'](Renzyme): 
                                total_passing_cands += 1
                        else:
                            total_passing_cands +=1
                        yield passed_cand_dict
        if get_filter_stats:
            with open(locusname + '_filter_stats.tsv', 'w') as fout:
                for filtername in filter_stats.keys():
                    fout.write(
                        filtername
                        + '\t'
                        + str(round((filter_stats[filtername]/total_cands) * 100, 2))
                        + '\n'
                    )
                fout.write('total_candidates\t' + str(total_cands) + '\n')
                fout.write('total_passing_candidates\t' + str(total_passing_cands))



    "ENDOGENOUS LEVEL ENZYME FILTERING METHODS"

    def no_endogenous_cutsites_within_Nnt_of_primers(self, Renzyme, N=5):
        """
        Filters out Renzymes for all candidates that have endogenous sites
        within N nucleotides of the primers.

        Parameters
        ----------
        Renzyme : Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction enzyme object

        N: int
            The desired minimum distance in nucleotides an endogenous site must
            be from the primers for the Renzyme to pass this filter. Default
            = 5

        Returns
        -------
        bool
            False if Renzyme has an endogenous site within Nnt of one of the
            primers, True otherwise. This function will be easier to write when
            we've fleshed out exactly what sequences we are giving to Endogenous
            and Candidate

        """
        pass

    def does_not_exceed_max_cutsites(self, Renzyme, N):
        """
        Filters out this Renzyme for all candidates it it has more than
        N cutsites in the endogenous sequence.

        Parameters
        ----------
        N : int
            The maximum number of endogenous cutsites Renzyme is allowed to
            have.

        Returns
        -------
        bool
            False if Renzyme has more than N cutsites in the endogenous sequence,
            True otherwise.
        """
        if len(self.Renzymes[Renzyme]) > N:
            return False
        else:
            return True


class Candidate(Designer):
    """
    The Candidate class is a subclass of the designer class, and
    a candidate is a candidate repair oligo with synonymous substitutions
    to the endogenous sequence within 15 nt of the Cas9 cutsite.
    """

    def __init__(self, endogenous, sequence):

        """
        Parameters
        ----------

        endogenous : Endogenous object
            An Endogenous object, the wildtype sequence that the candidate is
            based on.

        sequence : string
            The candidate sequence, very similar to the endogenous object
            sequence but it has the desired CRISPR/cas9 nonsynonymous
            substitution, as well as a series of synonymous substitutions
            within the cutsite region.

        """

        self.endogenous = endogenous
        """Endogenous object : An Endogenous object, the wildtype sequence that the Candidate is
        based on.
        """
        super().__init__(
            sequence,
            self.endogenous.cutsite,
            self.endogenous.Renzyme_list
        )
        self.sub_dep_Renzymes = self.get_substitution_dependent_Renzymes()
        """dict : A dictionary of Restriction enzymes with cutsites that
        depend on synonymous substitutions in cutsite_region. This means
        that the bp substitutions in the Candidate sequence created recognition
        sites in the candidate sequence for these enzymes that aren't in the
        endogenous sequence. In other words, this is a dictionary of
        enzymes with recognition sites that are in the Candidate sequence
        but not in the endogenous sequence.  If an enzyme has recognition
        sites in the endogenous sequence AND new additional substitution
        dependent sites in the candidate sequence, only the novel
        substitution dependent cutsites will be listed for that enzyme in
        this dictionary.
        """
        self.sub_dep_Renzyme_recognition_sites = self.get_sub_dep_Renzyme_recognition_sites()
        """ dict : very similar to sub_dep_Renzymes, but the dictionary values are the
        coordinates of the enzyme recognition sites rather than the
        cutsites."""

        self.num_subs = self.get_hamming_distance()
        """int : the number of bp substitutions the candidate has compared to the
        wildtypes sequence (endogenous.sequence)"""
        
        self.endogenous_matches = self.get_endogenous_matches()
        """dict : a dictionary of the restriction enzymes with values describing the
        "polarity" of the endogenous site(s) compared to the substitution
        dependent site(s). This dictionary will assign 'None' restriction
        enzymes that have sandwiched endogenous or candidate sites (see
        method 'no_sandwhiched_cutsites' for more details. For polarity,
        '0' means that the restriction enzyme has no endogenous cutsites.
        A '+' value means that the restriction enzyme has one or more
        endogenous sites after the substitution dependent cutsite(s) in the
        sequence string (and none before it). A '-' value means that it has
        one or more endogenous cutsites before the substitution dependent
        cutsite(s) in the string (and none after it)"""

    "METHODS THAT SET ATTRIBUTES"

    def get_hamming_distance(self):
        """
        Returns the hamming distance between the endogenous sequence and the
        candidate sequence, AKA the number of single character substitutions in
        the candidate (self.sequence) compared to the endogenous sequence
        (self.endogenous.sequence)

        Returns
        -------
        int
            Integer indicating the number of substitutions to the endogenous
            sequence are in the candidate sequence.
        """

        return Levenshtein.distance(self.endogenous.sequence, self.sequence)

    def get_substitution_dependent_Renzymes(self):
        """Finds restriction enzyme cutsites that depend on synonymous
        substitutions in the cutsite region (AKA they are in the candidate
        sequence but not the endogenous sequence because one or more
        substitutions created their recognition site) and returns
        a dictionary mapping those cutsites to their enzymes.

        Returns
        -------
        dict
            A dictionary with enzyme names as keys and their substitution dependent
            cutsites as values.  If an enzyme has multiple cutsites in the
            sequence, it will only return the cutsites that are unique to the
            candidate sequence (compared to the endogenous
            self.endogenous.sequence)
        """

        sub_dep_Renzymes = {}
        for Renzyme in self.Renzymes:
            if Renzyme in self.endogenous.Renzymes:
                candidate_unique_sites = list(
                    set(self.Renzymes[Renzyme])
                    .difference(self.endogenous.Renzymes[Renzyme])
                )
                if len(candidate_unique_sites) > 0:
                    sub_dep_Renzymes[Renzyme] = candidate_unique_sites
            else:
                sub_dep_Renzymes[Renzyme] = self.Renzymes[Renzyme]
        return sub_dep_Renzymes

    def get_endogenous_matches(self):

        """This function will look at restriction enzymes with endogenous and
        substitution dependent sites, and it'll return a dictionary of the
        restriction enzymes with values describing the "polarity" of the
        endogenous site(s) compared to the substitution dependent site(s).

        Returns
        -------
        dict
            A dictionary with enzymes as keys and symbols indicating the polarity
            of their cutsites as values. A '0' value means that the restriction
            enzyme has no endogenous cutsites. A '+' value means that the
            restriction enzyme has one or more endogenous sites after the
            substitution dependent cutsite(s) in the sequence string (and none
            before it). A '-' value means that it has one or more endogenous
            cutsites before the substitution dependent cutsite(s) in the string
            (and none after it). A value of None means that that enzyme either has
            novel site(s) sandwiched between endogenous sites, or endogenous
            site(s) sandwiched between novel sites.
        """

        endogenous_matches = {}
        for Renzyme in self.sub_dep_Renzymes:
            if Renzyme not in self.endogenous.Renzymes:
                endogenous_matches[Renzyme] = '0'
            elif (
                max(self.sub_dep_Renzymes[Renzyme])
                < min(self.endogenous.Renzymes[Renzyme])
            ):
                endogenous_matches[Renzyme] = '+'
            elif (
                min(self.sub_dep_Renzymes[Renzyme])
                > max(self.endogenous.Renzymes[Renzyme])
            ):
                endogenous_matches[Renzyme] = '-'
            else:
                endogenous_matches[Renzyme] = None
        return endogenous_matches

    def get_sub_dep_Renzyme_recognition_sites(self):
        """Finds restriction enzyme recognition sites within the candidate seq
        that depend on synonymous substitutions in the cutsite region (AKA they
        are in the candidate sequence but not the endogenous sequence because
        one or more substitutions created them) and returns a dictionary
        mapping those recognition site coordinates to their enzymes.

        Returns
        -------
        dict
            A dictionary with enzyme names as keys and their substitution dependent
            recognition site coordinates as values.  If an enzyme has multiple
            recognition sites in the sequence, it will only return the recognition
            site coordinates that are unique to the candidate sequence (compared to
            the endogenous self.endogenous.sequence)
        """
        unfiltered_candidate_dict = utils.get_Rsite_dict(
            self.sequence,
            self.sub_dep_Renzymes
        )
        endogenous_dict = utils.get_Rsite_dict(
            self.endogenous.sequence,
            self.sub_dep_Renzymes
        )
        filtered_dict = {}
        for Renzyme in unfiltered_candidate_dict:
            filtered_dict[Renzyme] = []
            if len(unfiltered_candidate_dict[Renzyme]) > 0:
                if len(endogenous_dict[Renzyme]) > 0:
                    candidate_unique_sites = list(
                        set(unfiltered_candidate_dict[Renzyme])
                        .difference(endogenous_dict[Renzyme])
                    )
                    if len(candidate_unique_sites) > 0:
                        filtered_dict[Renzyme] = candidate_unique_sites
                else:
                    filtered_dict[Renzyme] = unfiltered_candidate_dict[Renzyme]
        return filtered_dict

    def get_average_distance_between_RE_cutsites_and_Cas9_cutsite(self, Renzyme):
        """
        Returns
        -------
        float
            The average distance in nucleotides between novel cutsite and cas9
            cutsite for each novel cutsite for Renzyme in the candidate sequence.
        """
        REdistances = [
            abs(i - 1 - self.cutsite) for i in self.sub_dep_Renzymes[Renzyme]
        ]
        return sum(REdistances)/len(REdistances)

    def get_distance_from_nonsynonymous_codon_to_Renzyme_cutsite(self, Renzyme):
        AA_change_site = self.endogenous.AA_change_site
        Renzyme_cutsites = [x - 1 for x in self.Renzymes[Renzyme]]
        REdistances = []
        for cutsite in Renzyme_cutsites:
            if cutsite in range(AA_change_site, AA_change_site + 3):
                REdistances.append(0)
            if cutsite < AA_change_site:
                REdistances.append(AA_change_site - cutsite)
            if cutsite >= AA_change_site + 3:
                REdistances.append(cutsite - (AA_change_site + 3))
        return min(REdistances)

    "CANDIDATE FILTERING METHODS"

    def more_than_one_substitution(self):
        """
        Filters out candidates that have less than two substitutions.

        Returns
        -------
        bool
            False if the candidate only has one substitution, True otherwise
        """
        if self.num_subs == 1:
            return False
        else:
            return True

    def no_new_cutsites(self):
        """
        Filters out candidates for which the synonymous changes did not create
        any new enzyme recognition sites.

        Returns
        -------
        bool
            False if the synonymous changes did not create
            any new enzyme recognition sites, True otherwise"""
        if len(self.sub_dep_Renzymes) == 0:
            return False
        else:
            return True

    def no_splicing_interference(self):
        """
        Filters out candidates for which the substitutions interferred with
        splicing hexamers.

        Returns
        -------
        bool
            False if there are any differences between self.splicing_hexamers and
            self.endogenous.splicing_hexamers (which would indicate that the
            candidate's synonymous changes interfered with the splicing hexamers),
            True otherwise.
        """
        if self.splicing_hexamers == self.endogenous.splicing_hexamers:
            return True
        else:
            return False

    def PAM_is_disrupted(self):
        """
        Filters out candidates for which the PAM site could be interrupted by
        a synonymous change or the desired AA change, but isn't.

        Returns
        -------
        bool
            False if the PAM site can be interrupted but isn't, True otherwise.
        """
        PAM_start = self.endogenous.PAM_start
        which_PAM = self.endogenous.which_PAM
        if self.endogenous.PAM_can_change:
            if which_PAM == 'CCN':
                if self.sequence[PAM_start:PAM_start + 2] == 'CC':
                    return False
            
            if which_PAM == 'NGG':
                if self.sequence[PAM_start + 1: PAM_start + 3] == 'GG':
                    return False
            
            if which_PAM == 'TTTV':
                if (self.sequence[PAM_start: PAM_start + 3] == 'TTT'
                and self.sequence[PAM_start + 3] != 'T'):
                    return False
            
            if self.endogenous.which_PAM == 'BAAA':
                if (self.sequence[PAM_start + 1: PAM_start + 4] == 'AAA'
                and self.sequence[PAM_start] != 'A'):
                    return False
        return True

    def miRNA_target_interference(self):
        pass

    "CANDIDATE LEVEL ENZYME FILTERING METHODS"

    # This one is also a helper function for get_endogenous_matches()
    def no_sandwiched_cutsites(self, Renzyme):
        """
        Filters out Renzymes for this candidate that have either one or more
        substitution dependent sites sandwiched between endogenous sites, or
        one or more endogenous sites sandwiched between substitution dependent
        sites.

        Parameters
        ----------
        Renzyme: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        Returns
        -------
        bool
            False if the Renzyme has any sandwiched cutsites, True otherwise
        """

        if self.endogenous_matches[Renzyme] is not None:
            return True
        else:
            return False

    def non_ambiguous_cutsites(self, Renzyme):
        """
        Filters out enzymes for which the exact location of the cutsite is
        not fully known.

        Parameters
        ----------
        Renzyme: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        Returns
        -------
        bool
            False if Renzyme has an unknown/ambiguous cutsite, True otherwise.
        """
        if '?' in Renzyme.elucidate():
            return False
        else:
            return True

    def no_endogenous_cutsites_within_Nnt_of_novel_cutsite(self, Renzyme, N=5):
        """
        Filters out Renzymes for this candidate that have endogenous cutsites
        within N nucleotides of a substitution dependent cutsite.

        Parameters
        ----------
        Renzyme: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        N: int
            The minimum distance in nucleotides any endogenous site can be from
            any substitution dependent site.

        Returns
        -------
        bool
            False if the Renzyme has an endogenous site within N nt of
            a novel site for this candidate, True otherwise.
        """
        if Renzyme not in self.endogenous.Renzymes:
            return True
        for cutsite in self.sub_dep_Renzymes[Renzyme]:
            distances = [abs(x - cutsite) for x in self.endogenous.Renzymes[Renzyme]]
            if min(distances) <= N:
                return False
        return True

    def no_lost_endogenous_sites(self, Renzyme):
        """
        Filters out Renzymes for this candidate that have endogenous sites that
        have been interrupted by substitutions, such that they only appear in
        the endogenous sequence and not in the candidate sequence.

        Parameters
        ----------
        Renzyme: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        Returns
        -------
        bool
            False if any endogenous cutsites for
            Renzyme are lost with the synonymous substitutions, True
            otherwise
        """
        if Renzyme in self.endogenous.Renzymes:
            if Renzyme not in self.Renzymes:
                return False
            for cutsite in self.endogenous.Renzymes[Renzyme]:
                if cutsite not in self.Renzymes[Renzyme]:
                    return False
        return True

    def novel_site_depends_on_at_least_two_substitutions(self, Renzyme):

        """
        Filters out Renzymes for this candidate that have substitution
        dependent recognition sites in the candidate sequence that depend on
        less than two substitutions.

        Parameters
        ----------
        Renzyme: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        Returns
        -------
        bool
            False if any of the substitution dependent recognition sites for this
            Renzyme in this candidate depend on less than two substitutions, True
            otherwise.
        """
        # get the locations of the substitutions that each recognition site for
        # Renzyme depend on
        Renzyme_subs = utils.get_substitutions_Renzyme_depends_on(
            Renzyme,
            self.endogenous.sequence,
            self.sequence,
            self.sub_dep_Renzyme_recognition_sites[Renzyme]
        )

        # If there's a site with only one sub, we don't want that RE!
        if all(
            [
                len(Renzyme_subs[i]) >= 2
                for i in range(0, len(Renzyme_subs))
            ]
        ):
            return True
        return False

    def num_recognition_sites_equals_num_cutsites(self, Renzyme):
        """
        Filters out any Renzymes for this candidate for which the number of
        recognition sites in the whole sequence does not equal the number of
        cutsites for the whole sequence.

        Parameters
        ----------
        Renzyme: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        Returns
        -------
        bool
            False if there are a different number of
            cutsites from recognition sites for Renzyme (which could mean that the
            enzyme cuts outsite the recognitition site and the Renzyme is close
            enough to the edge that the cut would not appear in the sequence, or
            maybe the enzyme cuts twice), True otherwise"""

        full_Rsite_dict = utils.get_Rsite_dict(self.sequence, self.Renzymes)
        if Renzyme in self.Renzymes and Renzyme in full_Rsite_dict:
            if len(self.Renzymes[Renzyme]) == len(full_Rsite_dict[Renzyme]):
                return True
        return False

    "ENZYME PAIR FILTERING METHODS"

    def at_least_one_has_endogenous_site(self, RenzymeA, RenzymeB):
        """
        Filters out Renzyme pairs for this candidate for which neither Renzyme
        in the pair has an endogenous site.

        Parameters
        ----------
        RenzymeA: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        RenzymeB: Bio.Restriction enzyme object
            A restriction enzyme other than RenzymeA, in the form of
            a Bio.Restriction object.

        Returns
        -------
        bool
            False if neither RenzymeA nor RenzymeB have
            any endogenous sites (AKA if both have polarities of '0'), True
            otherwise.
        """
        if (
            self.endogenous_matches[RenzymeA] == '0'
            and self.endogenous_matches[RenzymeB] == '0'
        ):
            return False
        else:
            return True

    def both_have_endogenous_sites(self, RenzymeA, RenzymeB):
        """
        A harsher version of 'at_least_one_has_endogenous_site', this filters
        out Renzyme pairs for this candidate for which at least one of the
        candidates in the pair does not have an endogenous site.

        Parameters
        ----------
        RenzymeA: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        RenzymeB: Bio.Restriction enzyme object
            A restriction enzyme other than RenzymeA, in the form of
            a Bio.Restriction object.

        Returns
        -------
        bool
            False unless both Renzymes have an endogenous site
        """

        if (
            self.endogenous_matches[RenzymeA] == '0'
            or self.endogenous_matches[RenzymeB] == '0'
        ):
            return False
        else:
            return True

    def compatible_polarity(self, RenzymeA, RenzymeB):
        """
        Filters out Renzyme pairs for this candidates if they have endogenous
        sites on opposite sides of their novel sites with respect to the
        beginning of the sequence (AKA one of the primers).

        Parameters
        ----------
        RenzymeA: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        RenzymeB: Bio.Restriction enzyme object
            A restriction enzyme other than RenzymeA, in the form of
            a Bio.Restriction object.

        Returns
        -------
        bool
            This function will return False if Renzyme A and Renzyme B have
            incompatible polarities (e.g. if one has a '-' polarity and the other
            has a '+' polarity), True otherwise.
        """

        if (
            self.endogenous_matches[RenzymeA] == '+'
            and self.endogenous_matches[RenzymeB] == '-'
        ):
            return False

        if (
            self.endogenous_matches[RenzymeA] == '-'
            and self.endogenous_matches[RenzymeB] == '+'
        ):
            return False

        else:
            return True

    def different_recognition_sites(self, RenzymeA, RenzymeB):
        """
        Filters out Renzyme pairs for this candidate that have the same
        recognition site.

        Parameters
        ----------
        RenzymeA: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        RenzymeB: Bio.Restriction enzyme object
            A restriction enzyme other than RenzymeA, in the form of
            a Bio.Restriction object.

        Returns
        -------
        bool
            False if RenzymeA and RenzymeB have the same recognition site"""
        if RenzymeA.is_isoschizomer(RenzymeB):
            return False
        else:
            return True

    def both_Renzymes_depend_on_at_least_one_unique_substitution(
        self,
        RenzymeA,
        RenzymeB
    ):
        """
        Filters out Renzyme pairs for this candidate for which both Renzymes
        depend on the exact same set of substitutions for their substitution
        dependent recognition sites.

        Parameters
        ----------
        RenzymeA: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        RenzymeB: Bio.Restriction enzyme object
            A restriction enzyme other than RenzymeA, in the form of
            a Bio.Restriction object.

        Returns
        -------
        bool
            True only if RenzymeA and RenzymeB each depend on at least one
            substitution that the other does not depend on, False otherwise.
        """
        RenzymeA_subs = utils.get_substitutions_Renzyme_depends_on(
            RenzymeA,
            self.endogenous.sequence,
            self.sequence,
            self.sub_dep_Renzyme_recognition_sites[RenzymeA]
        )
        RenzymeB_subs = utils.get_substitutions_Renzyme_depends_on(
            RenzymeB,
            self.endogenous.sequence,
            self.sequence,
            self.sub_dep_Renzyme_recognition_sites[RenzymeB]
        )

        N = 0
        for Atuple in RenzymeA_subs:
            for Btuple in RenzymeB_subs:
                if len(set(Atuple).difference(Btuple)) > 0 and len(set(Btuple).difference(Atuple)) > 0:
                    N += 1
        if N > 0:
            return True
        else:
            return False

    def no_unneeded_substitutions(self, RenzymeA, RenzymeB=None):
        """
        Filters out Renzyme pairs for this candidate if there are any
        synonymous substitutions that are not needed to create the novel sites
        for Renzyme A and B and the desired AA change. If the PAM site cannot
        be changed, it will not exclude candidates with unneeded substitutions
        that are within the gRNA target region as long as they have no unneeded
        substitutions outside the gRNA target region."

        Parameters
        ----------
        RenzymeA: Bio.Restriction enzyme object
            A restriction enzyme in the form of a Bio.Restriction object.

        RenzymeB: Bio.Restriction enzyme object
            A restriction enzyme other than RenzymeA, in the form of
            a Bio.Restriction object. To run this filter only on RenzymeA, do
            not specify a RenzymeB

        Returns
        -------
        bool
            False if there are any unneeded substitutions (unless the PAM can't be
            changed, in which case there can be unneeded substitutions as long as
            they are within the gRNA target region), True otherwise.
        """

        # get a flattened list of all the substitutions RenzymeA and Renzyme
        # B depend on
        RenzymeA_subs = utils.get_substitutions_Renzyme_depends_on(
            RenzymeA,
            self.endogenous.sequence,
            self.sequence,
            self.sub_dep_Renzyme_recognition_sites[RenzymeA]
        )
        flattened_RenzymeA_subs = list(sum(RenzymeA_subs, ()))
        flattened_RenzymeB_subs = flattened_RenzymeA_subs
        if RenzymeB is not None:
            RenzymeB_subs = utils.get_substitutions_Renzyme_depends_on(
                RenzymeB,
                self.endogenous.sequence,
                self.sequence,
                self.sub_dep_Renzyme_recognition_sites[RenzymeB])
            flattened_RenzymeB_subs = list(sum(RenzymeB_subs, ()))

        # for each character in self.sequence, check if there's a substitution.
        for i in range(0, len(self.sequence)):
            if self.endogenous.sequence[i] != self.sequence[i]:
                # if there is, check if the substitution is outside of the
                # codon yielding the nonsynonymous mutation
                if i not in range(
                        self.endogenous.AA_change_site,
                        self.endogenous.AA_change_site + 3
                ):
                    # then, check to make sure the substitution is not needed
                    # for a Renzyme recognition site
                    if (
                        i not in flattened_RenzymeA_subs
                        and i not in flattened_RenzymeB_subs
                    ):
                        # if the PAM can change, return False if the candidate
                        # has any unneeded substitutions
                        if self.endogenous.PAM_can_change:
                            if self.endogenous.which_PAM == 'CCN':
                                if i not in range(
                                    self.endogenous.PAM_start,
                                    self.endogenous.PAM_start + 2
                                ):
                                    return False
                            if self.endogenous.which_PAM == 'NGG':
                                if i not in range(
                                    self.endogenous.PAM_start + 1,
                                    self.endogenous.PAM_start + 3
                                ):
                                    return False
                            if (
                                    self.endogenous.which_PAM == 'TTTV'
                                    or self.endogenous.which_PAM == 'BAAA'
                            ):
                                if i not in range(
                                    self.endogenous.PAM_start,
                                    self.endogenous.PAM_start + 4
                                ):
                                    return False
                        # if the PAM can't change, only exclude candidates with
                        # unneeded substitutions outside the gRNA
                        else:
                            if i not in range(
                                self.endogenous.gRNA_location[0],
                                self.endogenous.gRNA_location[1]
                            ):
                                return False
        return True

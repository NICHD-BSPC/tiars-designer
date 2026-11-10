import argparse
import pandas as pd
import csv
from Bio import Restriction
from Bio.Restriction import *
import itertools
import gzip
import sys
import xlwt
import openpyxl
import numpy as np
from docx import Document
from docx.enum.text import *
from docx.shared import RGBColor
from xlwt import Workbook
from .__init__ import *



def get_enzyme_pair_stats_filtered(passed_enzymes_df, endogenous, locus, enzyme_pair_filters, codon_frequency_dict, verbose):

    """
    Given the candidates with enzymes that passed all four
    enzyme level filters, runs another set of filters on every possible pair
    of enzymes remaining for each candidate.

    Parameters
    ----------
    passed_enzymes_df: pandas dataframe
        A pandas dataframe that contains the passed candidate-enzyme combos. It
        should have a sequence column, a Renzyme column, a codon frequency
        column and an average distance to cas9 cutsite column.

    endogenous: __init__.Endogenous object
        An endogenous object that all of the remaining candidates are based on.

    locus: string
        The name of the locus; to be used in the csv file name if verbose ==
        True
    enzyme_pair_filters: list
        A list of enzyme pair filters to use
    codon_frequency_dict: dict
        A dictionary with amino acid codons (e.g. 'AGA') as keys and the
        frequency with which they are used.
    verbose: bool
        A boolean value indicating whether or not to create a csv file from the
        output dataframe

    Returns
    -------
    A pandas dataframe with one row for each pair
    of passed enzymes for each candidate in passed_enzymes_df. The columns are
    'seq' (sequence), Renzyme pair (a  pair of Renzymes that each have at least
    one substitution dependent cutsite in its corresponding candidate
    sequence),
    result (boolean value indicating whether or not the Renzyme pair passed the
    filter for that row's candidate), RenzymeA average distance to cas9 cut,
    RenzymeB average distance to cas9 cut, overall average distance to cas9 cut
    (average of previous 2), and codon frequency change. If the verbose option
    is selected, it will also create a csv called
    '[locus]_pair_filter_stats' from the dataframe
    """

    passed_pairs = []
    seq_enzymes_key = pd.Series(
        passed_enzymes_df.Renzyme.values,
        index=passed_enzymes_df.index
    )
    seq_enzymes_key = {
        k: RestrictionBatch(v.strip().split(', '))
        for (k, v) in seq_enzymes_key.items()
    }
    if 'both_have_endogenous_sites' in enzyme_pair_filters and 'at_least_one_has_endogenous_site' not in enzyme_pair_filters:
        enzyme_pair_filters.append('at_least_one_has_endogenous_site')
    if verbose:
        filter_stats = {Filter: 0 for Filter in enzyme_pair_filters}
        total_cands = 0
    for candidate_seq in seq_enzymes_key:
        # make a Candidate object for the sequence and
        # get a list of all possible pairs of passed enzymes,
        # then run the enzyme pair filters on them
        if len(seq_enzymes_key[candidate_seq]) > 0:
            cand = Candidate(endogenous, candidate_seq)
            filter_functions = utils.get_filter_function_dict(enzyme_pair_filters, cand)
            enzyme_pairs = itertools.combinations(
                seq_enzymes_key[candidate_seq],
                2
            )

            for pair in enzyme_pairs:
                enzymeA = pair[0]
                enzymeB = pair[1]
                
                if verbose:
                    for Filter in enzyme_pair_filters:
                        filter_stats[Filter] += filter_functions[Filter](enzymeA, enzymeB)
                    total_cands += 1
                filter_results = {}
                for Filter in enzyme_pair_filters:
                    if Filter != 'both_have_endogenous_sites':
                        filter_results[Filter] = filter_functions[Filter](enzymeA, enzymeB)
                if all(filter_results.values()):
                    RenzymeA_avg_distance_from_cas9_cut = (
                        cand.get_average_distance_between_RE_cutsites_and_Cas9_cutsite(
                            enzymeA
                        )
                    )
                    RenzymeB_avg_distance_from_cas9_cut = (
                        cand.get_average_distance_between_RE_cutsites_and_Cas9_cutsite(
                            enzymeB
                        )
                    )
                    overall_avg_distance_from_cas9_cut = (
                        (
                            RenzymeA_avg_distance_from_cas9_cut
                            + RenzymeB_avg_distance_from_cas9_cut
                        )/2
                    )
                    RenzymeA_cutsite_distance_outside_Rsite = len(enzymeA.elucidate()) - 2 - len(enzymeA.site)
                    RenzymeB_cutsite_distance_outside_Rsite = len(enzymeB.elucidate()) - 2 - len(enzymeB.site)
                    RenzymeA_distance_to_AA_change = cand.get_distance_from_nonsynonymous_codon_to_Renzyme_cutsite(enzymeA)
                    RenzymeB_distance_to_AA_change = cand.get_distance_from_nonsynonymous_codon_to_Renzyme_cutsite(enzymeB)
                    codon_frequency_differences = utils.get_codon_frequency_differences(
                        endogenous.sequence,
                        candidate_seq,
                        codon_frequency_dict
                    )
                    passed_pairs.append(
                        {
                            'seq': candidate_seq,
                            'Renzyme pair': str(enzymeA) + ' and ' + str(enzymeB),
                            'RenzymeA avg distance from cas9 cut': RenzymeA_avg_distance_from_cas9_cut,
                            'RenzymeB avg distance from cas9 cut': RenzymeB_avg_distance_from_cas9_cut,
                            'overall avg distance from cas9 cut': overall_avg_distance_from_cas9_cut,
                            'codon frequency differences': codon_frequency_differences,
                            'RenzymeA cutsite distance outside Rsite': RenzymeA_cutsite_distance_outside_Rsite,
                            'RenzymeB cutsite distance outside Rsite': RenzymeB_cutsite_distance_outside_Rsite,
                            'RenzymeA distance to Amino Acid change': RenzymeA_distance_to_AA_change,
                            'RenzymeB distance to Amino Acid change': RenzymeB_distance_to_AA_change,
                            'both enzymes have endogenous site': cand.both_have_endogenous_sites(enzymeA, enzymeB)
                        }
                    )
    passed_pairs_df = pd.DataFrame(passed_pairs)
    if 'both_have_endogenous_sites' in enzyme_pair_filters:
        try:
            passed_pairs_harsh_df = passed_pairs_df[
                passed_pairs_df['both enzymes have endogenous site'] == True
            ]

        except:
            final_df =passed_pairs_df
        else:
            if len(passed_pairs_harsh_df) == 0:
                final_df= passed_pairs_df
            else:
                final_df = passed_pairs_harsh_df
    else: final_df=passed_pairs_df

    total_passing_cands = len(final_df)
    if verbose:
        with open(locus + '_pair_filter_stats.tsv', 'w') as fout:
            for filtername in filter_stats.keys():
                stat = '0'
                if total_cands != 0:
                    stat = str(round((filter_stats[filtername]/total_cands) * 100, 2))
                fout.write(
                    filtername
                    + '\t'
                    + stat
                    + '\n'
                )
            fout.write('total_candidates\t'
                       + str(total_cands)
                        + '\n')
            fout.write('total_passing_candidates\t'
                       + str(total_passing_cands)
                       )
    return final_df



def get_docx_single(flat_passed_enzymes_df, no_unneeded_substitutions, endogenous, locus):

    """
    Given a dataframe containing all of the candidate-enzyme combos that
    passed all of the filters, it further filters out any candidates with
    unneeded substitutions (if the filter is specified) and creates a docx file called
    '[locus]_single_Renzymes.docx to visualize all of the remaining candidate
    sequences and enzyme pairs. It will also create an xlsx file with more
    in-depth statistics on each candidate-Renzyme combo.

    Parameters
    ----------
    flat_passed_enzymes_df: pandas dataframe
        A pandas dataframe that contains the passed candidate-enzyme pair
        combos, with each combo getting its own row. It should have a sequence
        column, a Renzyme pair column, an average distance to cas9 cut column,
        a codon frequency change column, and a no unnneeded sustitutions
        column.
    no_unneeded_substitutions: bool
        Indicates whether the no_unneeded substitutions filter should be applied.
    endogenous: __init__.Endogenous object
        An endogenous object that all of the remaining candidates are based on.

    locus: string
        The name of the locus; to be used in the tsv file name if verbose ==
        True
    """
    if no_unneeded_substitutions:

        flat_passed_enzymes_df = flat_passed_enzymes_df[
            flat_passed_enzymes_df['no unneeded substitutions']
            == True
        ]
        flat_passed_enzymes_df=(
            flat_passed_enzymes_df
            .drop(columns=['no unneeded substitutions'])
        )
    # Sort the dataframe to put the candidates with the most desirable
    # statistics on top and then only write the first 500000 to file 
    flat_passed_enzymes_df = (
        flat_passed_enzymes_df
        .sort_values(by=['distance to Amino Acid change',
                  'average distance to Cas cut',
                  'codon frequency change',
                  'cutsite distance outside Rsite'])
        .reset_index(drop=True)
    )
    if len(flat_passed_enzymes_df) > 500000:
        flat_passed_enzymes_df = flat_passed_enzymes_df.iloc[0:500000]


    flat_passed_enzymes_df.to_excel(
                locus + '_single_Renzymes.xlsx',
                index=True,
                header=True
            )
    utils.make_doc(
        endogenous.sequence,
        endogenous.gRNA_location[0],
        endogenous.gRNA_location[1],
        endogenous.PAM_start,
        endogenous.AA_change_site,
        flat_passed_enzymes_df,
        locus,
        endogenous.which_PAM,
        False
    )

def full_run(
    endogenous_seq,
    gRNA_start,
    gRNA_end,
    PAM_start,
    which_PAM,
    AA_change_site,
    new_codon,
    locus,
    single_enzyme_filters,
    enzyme_pair_filters,
    codon_usage_tsv,
    min_Rsite_distance=5,
    Renzymes=None,
    verbose=False
):
    """
    Given an endogenous sequence, desired edit + guideRNA, generates all the
    possible synonymous candidates for that sequence and runs all of the
    filters on them to get the final docx and xlsx files (one of each for both
    single enzymes and enzyme pairs) using the above functions. This is the
    only function for which the users will directly input the parameters (via
    commandline argparse).

    Parameters
    ----------

    endogenous_seq: string
        The endogenous/wt dna sequence without any synonymous changes
        and without the desired CRISPR facilitated amino acid change. Max
        length is 273 bp, must be in frame.

    gRNA_start: int
        the index location in the sequence of the beginning of the
        region the CRISPR/Cas9 guide RNA will target.

    gRNA_end: int
        the index location in the sequence of the end of the region
        the CRISPR/Cas9 guide RNA will target.

    PAM_start: int
        The index location in the sequence of the beginning of the
        PAM site. For now, the PAM region has to be either 'CCN' or 'NGG'.

    which_PAM: string
        The sequence of the PAM site (e.g. 'NGG' or 'CCN' for cas9, or 'TTTV' or 'BAAA' for Cfp1)

    AA_change_site: int
        the index location in the sequence of the beginning of
        the codon that is to be nonsynonymously changed in the CRISPR/Cas9
        process

    new_codon: string
        the codon to replace the existing codon at the
        AA_change_site

    locus: string
        The name of the locus; to be used in the csv file name if verbose ==
        True

    single_enzyme_filters: list
        The number of single enzyme filters an enzyme needs to pass in order to
        stay in the running

    enzyme_pair_filters: list
        The number of paired enzyme filters a pair of enzymes needs to pass in
        order to stay in the running.

    codon_usage_tsv: tsv
        A tsv mapping codons to their usage frequency for the organism used
    
    min_Rsite_distance (optional): int
        The minimum allowable distance between the novel restriction site for
        an enzyme and any endogenous restriction sites it may have. This is to
        avoid having the fragment lengths after digestion be too similar to
        detect the difference between a novel site cut and an endogenous site
        cut. Default is 5 bp.

    enzymes: (optional) string
        A string indicating the path to a tsv file containing a list of enzymes
        to consider for the list of candidates. Default is the Bio.Restriction
        'AllEnzymes' restriction batch.

    verbose: (optional) bool
        A boolean value indicating whether or not to create tsv files with filter statistics at each step. Default is False
    """
    if single_enzyme_filters == ['all']:
        single_enzyme_filters = [
            'no_sandwiched_cutsites',
            'no_endogenous_cutsites_within_Nnt_of_novel_cutsite',
            'novel_site_depends_on_at_least_two_substitutions',
            'num_recognition_sites_equals_num_cutsites',
            'no_unneeded_substitutions',
        ]
    if enzyme_pair_filters == ['all']:
        enzyme_pair_filters = [
            'at_least_one_has_endogenous_site',
            'both_have_endogenous_sites',
            'compatible_polarity',
            'both_Renzymes_depend_on_at_least_one_unique_substitution'
        ]
    Renzyme_list = None
    # build a restriction batch out of the enzyme tsv if one was provided
    if Renzymes is not None:
        Renzyme_list = RestrictionBatch([])
        equischizomers = {}
        for line in open(Renzymes):
            Renzyme = getattr(Restriction, line.strip())
            if Renzyme.elucidate() in equischizomers:
                equischizomers[Renzyme.elucidate()].append(str(Renzyme))
            else:
                equischizomers[Renzyme.elucidate()] = [str(Renzyme)]
                Renzyme_list.add(Renzyme)
        with open('equischizomers.tsv', 'w') as fout:
            for equischizomer in equischizomers:
                newline = equischizomer + '\t' + (', ').join(equischizomers[equischizomer]) + '\n'
                fout.write(newline)
    
    # Trim down the endogenous sequence and adjust the coordinates accordingly
    trimmed = utils.get_endogenous_sequence_and_stats(
        endogenous_seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        AA_change_site,
        which_PAM)
    
    endogenous_seq = trimmed['sequence']
    gRNA_start = trimmed['gRNA_start']
    gRNA_end = trimmed['gRNA_end']
    PAM_start = trimmed['PAM_start']
    AA_change_site = trimmed['AA_change_site']

    # Initialize the endogenous object
    endogenous = Endogenous(
        endogenous_seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_codon,
        Renzyme_list
    )
#    # Get a dictionary of the frequency for each codon for output stats
    codon_frequency_dict = utils.get_codon_frequency_dict(codon_usage_tsv)

    # Get a prefiltered list of candidate objects for this endogenous object
    flat_passed_enzymes_df = pd.DataFrame(endogenous.get_filtered_candidates(locus, single_enzyme_filters, min_Rsite_distance, verbose, codon_usage_tsv))
    # If no solutions were found, write out a no solutions error message to the output files
    SE_no_solutions_message = 'Unfortunately, no single enzyme solutions were found for this locus'
    EP_no_solutions_message = 'Unfortunately, no enzyme pair solutions were found for this locus'
    
    if len(flat_passed_enzymes_df) == 0:
        SE_doc = Document()
        SE_no_solutions = SE_doc.add_paragraph().add_run(SE_no_solutions_message)
        SE_doc.save(locus + '_single_Renzymes.docx')


        EP_doc = Document()
        EP_no_solutions = EP_doc.add_paragraph().add_run(EP_no_solutions_message)
        EP_doc.save(locus + '_Renzyme_pairs.docx')
        import os
        print(list(os.walk(locus + '*')))
        
        SE_wb = Workbook()
        SE_sheet=SE_wb.add_sheet('sheet')
        SE_sheet.write(0,0, SE_no_solutions_message)
        SE_wb.save(locus + '_single_Renzymes.xlsx')

        EP_wb = Workbook()
        EP_sheet=EP_wb.add_sheet('sheet')
        EP_sheet.write(0,0, EP_no_solutions_message)
        EP_wb.save(locus + '_Renzyme_pairs.xlsx')
    # If single enzyme solutions were found, join all enzymes that pass all filters to their single candidate
    else:
        passed_enzymes_df = pd.DataFrame(flat_passed_enzymes_df.groupby(['seq'])['Renzyme'].apply(', '.join))

    
        get_docx_single(flat_passed_enzymes_df, 'no_unneeded_substitutions' in single_enzyme_filters, endogenous, locus)

        passed_pairs_df = get_enzyme_pair_stats_filtered(
            passed_enzymes_df,
            endogenous,
            locus,
            enzyme_pair_filters,
            codon_frequency_dict,
            verbose
        )
        # If no enzyme pair candidates pass all filters, write out no solutions messages to the enzyme pair output files        
        if len(passed_pairs_df) == 0:
            EP_doc = Document()
            EP_ = EP_doc.add_paragraph().add_run(EP_no_solutions_message)
            EP_doc.save(locus + '_Renzyme_pairs.docx')

            EP_wb = Workbook()
            EP_sheet=EP_wb.add_sheet('sheet')
            EP_sheet.write(0,0, EP_no_solutions_message)
            EP_wb.save(locus + '_Renzyme_pairs.xlsx')

    # Otherwise, create docx files for passed single Renzymes and Renzyme pairs to map out
    # where the novel restriction sites are, and create excel files to output
    # their statistics 
        else:
            passed_pairs_df = passed_pairs_df.sort_values(
                by = ['RenzymeA distance to Amino Acid change',
                      'RenzymeB distance to Amino Acid change',
                      'overall avg distance from cas9 cut',
                      'codon frequency differences',
                      'RenzymeA cutsite distance outside Rsite',
                      'RenzymeB cutsite distance outside Rsite'],
                      ).reset_index(drop=True)
            if len(passed_pairs_df) > 500000:
                passed_pairs_df = passed_pairs_df.iloc[0:500000]

            passed_pairs_df.to_excel(
                locus + '_Renzyme_pairs.xlsx',
                index=True,
                header=True)
            utils.make_doc(
                endogenous.sequence,
                endogenous.gRNA_location[0],
                endogenous.gRNA_location[1],
                endogenous.PAM_start,
                endogenous.AA_change_site,
                passed_pairs_df,
                locus,
                endogenous.which_PAM
            )

def main():
    """
        A wrapper for arparse that will take in all the command line user inputs and run them through 'full_run()'

    """

    ap = argparse.ArgumentParser()
    ap.add_argument('endogenous_seq')
    ap.add_argument('gRNA_start', type=int)
    ap.add_argument('gRNA_end', type=int)
    ap.add_argument('PAM_start', type=int)
    ap.add_argument('which_PAM', type=str)
    ap.add_argument('AA_change_site', type=int)
    ap.add_argument('new_codon')
    ap.add_argument('locus_name')
    ap.add_argument('--single-enzyme-filters', nargs='+')
    ap.add_argument('--enzyme-pair-filters', nargs='+')
    ap.add_argument('--codon-usage-tsv')
    ap.add_argument('--min-Rsite-distance', type=int)
    ap.add_argument(
        '--enzymes',
        help='name of file containing single column list of enzymes to use'
    )    
    ap.add_argument(
        '--verbose',
        help='if True, this will write csvs/tsvs for each step',
        action="store_true")
    args = ap.parse_args()
    full_run(
        args.endogenous_seq,
        args.gRNA_start,
        args.gRNA_end,
        args.PAM_start,
        args.which_PAM,
        args.AA_change_site,
        args.new_codon,
        args.locus_name,
        args.single_enzyme_filters,
        args.enzyme_pair_filters,
        args.codon_usage_tsv,
        args.min_Rsite_distance,
        args.enzymes,
        args.verbose
    )


if __name__ == "__main__":
    main()

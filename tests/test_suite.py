import Bio
from Bio import Restriction
import tiarsdesigner

single_enzyme_filters = ['no_sandwiched_cutsites',
                         'no_endogenous_cutsites_within_Nnt_of_novel_cutsite',
                         'novel_site_depends_on_at_least_two_substitutions',
                         'num_recognition_sites_equals_num_cutsites',
                         'no_unneeded_substitutions'
                         ]
enzyme_pair_filters = ['at_least_one_has_endogenous_site',
                       'both_have_endogenous_sites',
                       'compatible_polarity',
                       'both_Renzymes_depend_on_at_least_one_unique_substitution',
                       'no_unneeded_substitutions'
                       ]

def get_enzyme(name):
    return getattr(Restriction, name)

def test_designer_class():
    designer = tiarsdesigner.Designer(
        'TCCGCTCACCAACTGTAAGGCCGG',
        16,
        [
            'BshFI',
            'PsuGI',
            'Sth302II',
            'FspEI',
            'HpyCH4III',
            'CviJI',
            'AciI',
            'SsiI',
            'SgrTI',
        ],
    )
    orig = designer.Renzymes


def test_endogenous_class_with_CCN_PAM():
    #in frame PAM and cutsite region

#          0         1         2         3         4         5         6         7         8         9         10
#          012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
#                                                 gRNA_start & PAM      AA_change gRNA_end
#                                                          |                    | |
    seq = 'GCTGGGGGCTCCACGGTGGTCAACATGAGCCATGCGGGCTGGTGGCCCCCAAATTGGTGGTGGTGGTGGTCACTCTGCAACACTCTTTACGTCAACTTCTGC'
    seq_with_codon_change = (
          'GCTGGGGGCTCCACGGTGGTCAACATGAGCCATGCGGGCTGGTGGCCCCCAAATTGGTGGTGGTGGTGGTtACTCTGCAACACTCTTTACGTCAACTTCTGC'
    )
#   NOTE: the AA_change_site coordinate marks the beginning of the codon
#   containing the substitution, not the substitution itself

    gRNA_start = 48
    gRNA_end = 71
    PAM_start = 48
    which_PAM = 'CCN'
    AA_change_site = 69
    new_AA_codon = 'TTA'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'PspOMII',
            'AquII',
            'TaqIII',
            'SetI',
            'Hpy188I'
       ],
    )
    
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    assert not endogenous.PAM_can_change
    assert len(list(endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False,'codon_frequency.tsv'))) == 9

    #out of frame PAM and cutsite region
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start -1,
        gRNA_end - 1,
        PAM_start - 1,
        which_PAM,
        AA_change_site,
        new_AA_codon
    )
    assert endogenous.PAM_can_change

    endogenous.get_filtered_candidates('test2', single_enzyme_filters, 5, True, 'codon_frequency.tsv')
    assert len(list(endogenous.get_filtered_candidates('test', single_enzyme_filters, 5,False, 'codon_frequency.tsv'))) == 104


    #last frame option for PAM and cutsite region

    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start - 2,
        gRNA_end - 2,
        PAM_start - 2,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'PspOMII',
            'AquII',
            'TaqIII',
            'Hpy188I',
            'SetI',
            'Hpy99XIII',
            'BshFI',
            'Bhe175II',
            'PsuGI',
            'Sth302II',
            'FspEI',
            'Cco14983VI',
            'Lpn12272I',
            'Lmo911II',
            'OspHL35III',
            'Bce3081I',
            'Asp103I',
        ],
    )
    assert endogenous.PAM_can_change
    assert len(list(endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv'))) == 96


def test_endogenous_class_with_NGG_PAM():
    #in frame PAM and cutsite region

#          0         1         2         3         4         5         6         7         8         9         10
#          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#                                                    gRNA_start            PAM  gRNA_end
#                                                        |                   |  |
    seq = 'TTGCATTTTGTGTGACAGATGGCTGCGATTCGTAAGAAGCTTGTAATCGTTGGATGGTGGTGGTGGGGGTGGTGGTGGTTACTCATCGTTTTCAGTAAAGACC'
#                                                                            |
#                                                                       AA_change
    seq_with_codon_change =(
          'TTGCATTTTGTGTGACAGATGGCTGCGATTCGTAAGAAGCTTGTAATCGTTGGATGGTGGTGGTGGaGaTGGTGGTGGTTACTCATCGTTTTCAGTAAAGACC'
    )
    gRNA_start = 46
    gRNA_end = 69
    PAM_start = 66
    which_PAM = 'NGG'
    AA_change_site = 66
    new_AA_codon = 'AGA'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'AspAMDIV',
            'EcoE1140I',
            'Sth132I',
            'HinP1I',
            'BshFI',
            'PsuGI',
            'Sth302II',
            'FspEI'
        ],
    )
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv')
    assert len(list(candidates)) == 32


    #different frame
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start - 1,
        gRNA_end -1,
        PAM_start -1,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'AspAMDIV',
            'EcoE1140I',
            'Sth132I',
            'HinP1I',
            'BshFI',
            'PsuGI',
            'Sth302II',
            'FspEI',
            'NmeA6CIII'
        ],
    )
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv')
    assert len(list(candidates)) == 138

    #different frame
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start -2,
        gRNA_end -2,
        PAM_start -2,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'AspAMDIV',
            'EcoE1140I',
            'Sth132I',
            'HinP1I',
            'BshFI',
            'PsuGI',
            'Sth302II',
            'FspEI',
            'CjeNV',
            'Nan12227I'
        ],
    )
    assert len(list(endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv'))) == 186

def test_endogenous_class_with_TTTV_PAM():
    # In frame PAM

#          0         1         2         3         4         5         6         7         8         9         10
#          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#                                                                         AA_change
#                                                                            |  
    seq = 'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGCACTTCCTCCCCATGTTTCTCTGGTCTTGGCCCTGGTGGTGGTGTTGGTGGGTGTGGGTTTAAACAT'
#                                                             |                        | 
#                                                            PAM                    gRNA_end    
#                                                        gRNA_start
    seq_with_codon_change = (
          'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGCACTTCCTCCCCATGTTTCTCTGGTCTTGGaCCTGGTGGTGGTGTTGGTGGGTGTGGGTTTAAACAT'
    )
    gRNA_start = 51
    gRNA_end = 76
    PAM_start = 51
    which_PAM = 'TTTV'
    AA_change_site = 66
    new_AA_codon = 'ACC'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'AmaCSI',
            'AceIII',
            'KpnNH25III',
            'Lcr047I',
            'GsuI',
            'LmnI',
            'BshFI',
            'PsuGI',
            'Sth302II',
            'FspEI'

]
    )
    assert endogenous.cutsite == 74
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv')
    assert len(list(candidates)) == 72

#   New frame

#          0         1         2         3         4         5         6         7         8         9         10
#          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#                                                                         AA_change
#                                                                            |  
    seq = 'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGCACTTCCTCCCCATTTTCCTCTGGTCTTGGCCCTGGTGGTGGTGTTGGTGGGTGTGGGTTTAAACAT'
#                                                            |                        | 
#                                                            PAM                    gRNA_end    
#                                                        gRNA_start
    seq_with_codon_change = (
          'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGCACTTCCTCCCCATTTTCCTCTGGTCTTGGaCCTGGTGGTGGTGTTGGTGGGTGTGGGTTTAAACAT'
    )
    gRNA_start = 50
    gRNA_end = 75
    PAM_start = 50
    which_PAM = 'TTTV'
    AA_change_site = 66
    new_AA_codon = 'ACC'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'AmaCSI',
            'AceIII',
            'KpnNH25III',
            'Lcr047I',
            'GsuI',
            'LmnI',
            'BshFI',
            'PsuGI',
            'Sth302II',
            'FspEI'

        ]
    )
    assert endogenous.cutsite == 73
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv')
    assert len(list(candidates)) == 72

#   New frame

#          0         1         2         3         4         5         6         7         8         9         10
#          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#                                                                         AA_change
#                                                                            |  
    seq = 'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGCACTTCCTCCCCATTTCCCTCTGGTCTTGGCCCTGGTGGTGGTGTTGGTGGGTGTGGGTTTAAACAT'
#                                                           |                        | 
#                                                          PAM                    gRNA_end    
#                                                        gRNA_start
    seq_with_codon_change = (
          'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGCACTTCCTCCCCATTTCCCTCTGGTCTTGGaCCTGGTGGTGGTGTTGGTGGGTGTGGGTTTAAACAT'
    )

    gRNA_start = 49
    gRNA_end = 74
    PAM_start = 49
    which_PAM = 'TTTV'
    AA_change_site = 66
    new_AA_codon = 'ACC'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'AceIII',
            'KpnNH25III',
            'Lcr047I',
            'GsuI',
            'LmnI',
            'BshFI',
            'PsuGI',
            'Sth302II',
            'FspEI'

        ]
    )
    assert endogenous.cutsite == 72
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv')
    assert len(list(candidates)) == 64

def test_endogenous_class_with_BAAA_PAM():

#          0         1         2         3         4         5         6         7         8         9         10
#          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#                                                                         AA_change
#                                                                            |  
    seq = 'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGTGGTGGTGGGATACTTGGTGGTGGTGGTGGCCCTGGTAAAGGTGTTGGTGGGTGTGGGTTTAAACAT'
#                                                             |                    |   | 
#                                                        gRNA_start              PAM  gRNA_end    
#                                                        
    seq_with_codon_change = (
          'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGTGGTGGTGGGATACTTGGTGGTGGTGGTGGaCCTGGTAAAGGTGTTGGTGGGTGTGGGTTTAAACAT'
    )

    gRNA_start = 51
    gRNA_end = 76
    PAM_start = 72
    which_PAM = 'BAAA'
    AA_change_site = 66
    new_AA_codon = 'ACC'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [ 
            'Eco72I',
            'PmaCI',
            'BfaSII',
            'Acc65I',
            'UnbI',
            'BmgT120I',
            'KpnI',
            'HaeI']
    )
    assert endogenous.cutsite == 49
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False,'codon_frequency.tsv')
    assert len(list(candidates)) == 2

#   New Frame
#
#          0         1         2         3         4         5         6         7         8         9         10
#          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#                                                                         AA_change
#                                                                            |  
    seq = 'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGTGGTGGTGGGATACTTGGTGGTGGTGGTGGCCCTGTAAAAGGTGTTGGTGGGTGTGGGTTTAAACAT'
#                                                            |                    |   | 
#                                                        gRNA_start              PAM  gRNA_end    
#                                                        
    seq_with_codon_change = (
          'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGTGGTGGTGGGATACTTGGTGGTGGTGGTGGaCCTGTAAAAGGTGTTGGTGGGTGTGGGTTTAAACAT'
    )
    gRNA_start = 50
    gRNA_end = 75
    PAM_start = 71
    which_PAM = 'BAAA'
    AA_change_site = 66
    new_AA_codon = 'ACC'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'Eco72I',
            'AcvI',
            'PmaCI',
            'BfaSII',
            'Acc65I',
            'UnbI',
            'BmgT120I',
            'KpnI',
            'HaeI'
        ]
    )
    assert endogenous.cutsite == 48
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv')
    assert len(list(candidates)) == 3

#   New Frame
#
#          0         1         2         3         4         5         6         7         8         9         10
#          0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#                                                                         AA_change
#                                                                            |  
    seq = 'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGTGGTGGTGGGATACTTGGTGGTGGTGGTGGCCCTGAAAAAGGTGTTGGTGGGTGTGGGTTTAAACAT'
#                                                           |                    |   | 
#                                                       gRNA_start              PAM  gRNA_end    
#                                                        
    seq_with_codon_change = (
          'AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGTGGTGGTGGGATACTTGGTGGTGGTGGTGGaCCTGAAAAAGGTGTTGGTGGGTGTGGGTTTAAACAT'
    )
    gRNA_start = 49
    gRNA_end = 74
    PAM_start = 70
    which_PAM = 'BAAA'
    AA_change_site = 66
    new_AA_codon = 'ACC'
    endogenous = tiarsdesigner.Endogenous(
        seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'Eco72I',
            'AcvI',
            'BbrPI',
            'PmaCI',
            'BfaSII',
            'Acc65I',
            'UnbI',
            'BmgT120I',
            'KpnI',
            'HaeI'
        ]
    )
    assert endogenous.cutsite == 47
    assert endogenous.AA_change_sequence == seq_with_codon_change.upper()
    candidates = endogenous.get_filtered_candidates('test', single_enzyme_filters, 5, False, 'codon_frequency.tsv')
    assert len(list(candidates)) == 8

def test_candidate_class():
    """I made the following example by putting together a bunch of
    restriction sites in orders that created different cases to test. Here are
    dictionaries of the Renzyme recognition sites for endog_seq and cand_seq
    created below for easier following of the test cases:

    endog_seq:  {
        UnbI: [],
        WviI: [(60, 66)],
        PsuGI: [(3, 8)],
        FmuI: [],
        AluI: [(25, 29)],
        BshFI: [(45, 49), (55, 59)],
        MvnI: [(15, 19)],
        HaeII: [(7, 13), (49, 55)]
    }

    cand_seq: {
        UnbI: [(16, 21)],
        WviI: [(60, 66)],
        PsuGI: [(3, 8), (13, 18)],
        FmuI: [(16, 21)],
        AluI: [(21, 25), (25, 29), (29, 33)],
        BshFI: [(17, 21), (45, 49), (55, 59)],
        MvnI: [],
        HaeII: [(7, 13), (36, 42), (49, 55)]
    }

    """
#                0         1         2         3         4         5         6
#                012345678901234567890123456789012345678901234567890123456789012345
#                          gRNA                 PAM gRNA  AA_change
#                          start                  | end    |
#                             |                   |  |     |
    endog_seq = 'AAATCGGAGCGCCTGCGCGCAAGTTAGCTGGCTAGGAGTGCCAGTGGCCAGCGCTGGCCTCACAAG'
    cand_seq =  'AAATCGGAGCGCCTGCGgGCcAGcTAGCTaGCTAGGAGcGCtAaTGGCCAGCGCTGGCCTCACAAG'
    gRNA_start = 13
    gRNA_end = 36
    PAM_start = 33
    which_PAM = 'NGG'
    AA_change_site = 42
    new_AA_codon = 'AAT'
    endogenous = tiarsdesigner.Endogenous(
        endog_seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'FmuI',
            'BshFI',
            'PsuGI',
            'HaeII',
            'AluI',
            'EaeI',
            'MvnI',
            'WviI',
            'UnbI',
        ],
    )

    candidate = tiarsdesigner.Candidate(
        endogenous,
        cand_seq.upper(),
    )

    FmuI = get_enzyme('FmuI')
    BshFI = get_enzyme('BshFI')
    PsuGI = get_enzyme('PsuGI')
    HaeII = get_enzyme('HaeII')
    AluI = get_enzyme('AluI')
    EaeI = get_enzyme('EaeI')
    MvnI = get_enzyme('MvnI')
    WviI = get_enzyme('WviI')
    UnbI = get_enzyme('UnbI')

    #METHODS THAT SET ATTRIBUTES

    #test get_substitution_dependent_Renzymes

    #no endogenous sites
    assert candidate.sub_dep_Renzymes[FmuI] == [21]
    #has endogenous sites (but they should not be listed here)
    assert candidate.sub_dep_Renzymes[BshFI] == [20]
    #multiple novel sites
    assert candidate.sub_dep_Renzymes[AluI] == [24, 32]
    #no novel sites
    assert EaeI not in candidate.sub_dep_Renzymes

    #Test 'get_endogenous_matches'

    assert candidate.endogenous_matches[BshFI] == '+'
    assert candidate.endogenous_matches[PsuGI] == '-'
    #no endogenous site
    assert candidate.endogenous_matches[FmuI] == '0'
    #sandwiched novel site
    assert candidate.endogenous_matches[HaeII] == None
    #sandwiched endogenous site
    assert candidate.endogenous_matches[AluI] == None

    #Test get_sub_dep_recognition_sites
    assert candidate.sub_dep_Renzyme_recognition_sites[FmuI] ==[(16, 21)]
    assert candidate.sub_dep_Renzyme_recognition_sites[BshFI] == [(17,21)]
    assert candidate.sub_dep_Renzyme_recognition_sites[AluI] == [(21,25),(29,33)]
    assert EaeI not in candidate.sub_dep_Renzyme_recognition_sites

    #CANDIDATE FILTERS
    assert not candidate.PAM_is_disrupted()

    #SINGLE ENZYME FILTERS

    #Test no_sandwiched_cutsites:

    #two endogenous sites on same site of novel site
    assert candidate.no_sandwiched_cutsites(BshFI)
    #sandwiched novel site
    assert not candidate.no_sandwiched_cutsites(HaeII)
    #sandwiched endogenous site
    assert not candidate.no_sandwiched_cutsites(AluI)

    #Test no_endogenous_cutsites_within_Nnt_of_novel_cutsite
    assert candidate.no_endogenous_cutsites_within_Nnt_of_novel_cutsite(HaeII, 5)
    assert not candidate.no_endogenous_cutsites_within_Nnt_of_novel_cutsite(HaeII, 15)
    assert not candidate.no_endogenous_cutsites_within_Nnt_of_novel_cutsite(AluI, 5)

    #Test no_lost_endogenous_sites
    assert not candidate.no_lost_endogenous_sites(MvnI)

    #Test novel_site_depends_on_at_least_two_substitutions

    assert candidate.novel_site_depends_on_at_least_two_substitutions(BshFI)
    assert not candidate.novel_site_depends_on_at_least_two_substitutions(AluI)
    #substitution occurs at IUPAC non ATGC letter in recognition site (e.g. R = 'A' or 'G')
    assert not candidate.novel_site_depends_on_at_least_two_substitutions(HaeII)

    #Test num_recognition_sites_equals_num_cutsites

    #WviI.elucidate() = 'CACRAGNNNNNNNNNNNNNNNNNNN_NN^N'
    assert not candidate.num_recognition_sites_equals_num_cutsites(WviI)

    #PAIR FILTERS

    #Test at_least_one_has_endogenous_site

    #both have endogenous site
    assert candidate.at_least_one_has_endogenous_site(BshFI, PsuGI)
    #only one has endogenous site
    assert candidate.at_least_one_has_endogenous_site(BshFI, FmuI)
    #neither has endogenous site
    assert not candidate.at_least_one_has_endogenous_site(FmuI, UnbI)

    #Test both_have_endogenous_sites

    #both have endogenous site
    assert candidate.both_have_endogenous_sites(BshFI, PsuGI)
    #only one has endogenous site
    assert not candidate.both_have_endogenous_sites(BshFI, FmuI)
    #neither has endogenous site
    assert not candidate.both_have_endogenous_sites(FmuI, UnbI)

    #Test compatible_polarity

    #opposite polarity
    assert not candidate.compatible_polarity(BshFI, PsuGI)
    #only one has endogenous site
    assert candidate.compatible_polarity

    #Test different_recognition_sites

    assert not candidate.different_recognition_sites(FmuI, UnbI)
    assert candidate.different_recognition_sites(FmuI, HaeII)

    #Test both_Renzymes_depend_on_at_least_one_unique_substitution

    assert not candidate.both_Renzymes_depend_on_at_least_one_unique_substitution(FmuI, UnbI)
    assert candidate.both_Renzymes_depend_on_at_least_one_unique_substitution(BshFI, AluI)

    #Test no_unneeded_substitutions
    assert not candidate.no_unneeded_substitutions(PsuGI, BshFI)
    assert not candidate.no_unneeded_substitutions(HaeII, AluI)

def test_candidate_class_on_original_example():
    #The original example they did by hand should pass all of the filters

    endog_seq = 'CAGTCTGTAGCTTTTTCATAAAGGTGTCTGTTATTGCTTGAGTTGCATTATGGGATGCTGTTGTTTGCAGGTGTGTTTATGCCTGTGGGGTTGGTTCTGCAGCCCTGGATGGGCTCTGCTGCTATGGCTCTATCTTCTGTGTCTGTAGTGCTGTCCTCCCTCCTGCTAAAATGGTAAACCTTTCAGATCCTAACATTTAGAGTTGAAATAGATGAAATTAATTAAACATGTTTACATGTGTTAACAGCGTTGATTATACTAATGTCAATTAAT'
    candidate_seq = 'CAGTCTGTAGCTTTTTCATAAAGGTGTCTGTTATTGCTTGAGTTGCATTATGGGATGCTGTTGTTTGCAGGTGTGTTTATGCCTGTGGGGTTGGTcCTtCAGagCTGGATGGGCTCTGCTGCTATGGCTCTATCTTCTGTGTCTGTAGTGCTGTCCTCCCTCCTGCTAAAATGGTAAACCTTTCAGATCCTAACATTTAGAGTTGAAATAGATGAAATTAATTAAACATGTTTACATGTGTTAACAGCGTTGATTATACTAATGTCAATTAAT'
    gRNA_start = 90
    gRNA_end = 113
    PAM_start = 110
    which_PAM = 'NGG'
    AA_change_site = 102
    new_AA_codon = 'AGC',
    endogenous = tiarsdesigner.Endogenous(
        endog_seq,
        gRNA_start,
        gRNA_end,
        PAM_start,
        which_PAM,
        AA_change_site,
        new_AA_codon,
        [
            'HpyAV',
            'Hpy188I',
        ],
    )
    candidate = tiarsdesigner.Candidate(
        endogenous,
        candidate_seq.upper(),
    )

    #candidate filters
    assert candidate.PAM_is_disrupted()

    HpyAV = get_enzyme('HpyAV')
    Hpy188I = get_enzyme('Hpy188I')

    #single enzyme filters
    assert candidate.no_sandwiched_cutsites(HpyAV) and candidate.no_sandwiched_cutsites(Hpy188I)
    assert candidate.non_ambiguous_cutsites(HpyAV) and candidate.non_ambiguous_cutsites(Hpy188I)
    assert candidate.no_endogenous_cutsites_within_Nnt_of_novel_cutsite(HpyAV) and candidate.no_endogenous_cutsites_within_Nnt_of_novel_cutsite(Hpy188I)
    assert candidate.no_lost_endogenous_sites(HpyAV) and candidate.no_lost_endogenous_sites(Hpy188I)
    assert candidate.novel_site_depends_on_at_least_two_substitutions(HpyAV) and candidate.novel_site_depends_on_at_least_two_substitutions(Hpy188I)
    assert candidate.num_recognition_sites_equals_num_cutsites(HpyAV) and candidate.num_recognition_sites_equals_num_cutsites(Hpy188I)

    #enzyme_pair_filters
    assert candidate.at_least_one_has_endogenous_site(HpyAV, Hpy188I)
    assert candidate.compatible_polarity(HpyAV, Hpy188I)
    assert candidate.different_recognition_sites(HpyAV, Hpy188I)
    assert candidate.both_Renzymes_depend_on_at_least_one_unique_substitution(HpyAV, Hpy188I)
    assert candidate.no_unneeded_substitutions(HpyAV, Hpy188I)
 #   assert length(list(candidates

def test_utils():

    #test get_endogenous_sequence_and_stats

    #'CCN' PAM
    full_seq = 'TTTTTTAAAAGCTAAATGATCATTTTTATAAACTTGCAAGTGATTTTAAGTGAGGTATAAGTGAATATATCTGTAATTTTTATATTGATATGCATCTATTTAGATTTATTATTATTTTTTTTTTTTGGTACATATAGTCTGTTGACATTTAAATTGTGCATTTTTTTTGGTCAAACCTAGAGATTTTGTCGAAACGGGTCCATCTTTAGAAGCAGCTCCCTGCCCTGATCTCTGTAGCAACACAACCACATAATGAATGACTTCTCTCGCTCTCTCCATCTCTGCGACTCCACCTATCTGCAGTGTTGTATAAGTAGAGTCTTCATCTTTCTCCTCAGTATTCTTCCAGGAGTTTGGACCGCTGATCTCCAGGATGGATAACGCCACGTACAGCTACGTGAACGACTGCACTCCGCTCACCAACTGTAAGGCCGGCCGTAAGGCTGGGGGCTCCACGGTGGTCAACATGAGCCATGCGGGCAAGAAGCCACCAAATGACTACCTGGTCTGGTCACTCTGCAACACTCTTTACGTCAACTTCTGCTGCCTGGGCTTCATGGCTCTGATCTACTCCATCAAGgtgagcttcagcaccgatatcttcatctgttttgtggctcattcaaaagtttttatttttattttttatgaatactcctgaaggtcacatgatgaaaagctctgattaacgctgaaaatggaccattgccatggaagacatgttttgaacagtgcaggatttttgaagtggtgtcattggagaaccatctttgatccttaaagaatccttaactgaacttaacaaatttaaaagaatccactattaggattttttaaaaataaagggttttatgaatggtaatgtttctgtcaattccgtgggattgttcagaaaagcatctgctattgctttagagaagccagaaacttgttgggaaatttacactcctctcttctgatttgtttgagaaatttacacttctctctcctgttttgtttgagaaatttacacttctcttttctgtt'

    seq_dict = tiarsdesigner.utils.get_endogenous_sequence_and_stats(
        full_seq,
        490,
        513,
        490,
        511,
        'CCN'
    )
    assert seq_dict == {
        'sequence':'TTCCAGGAGTTTGGACCGCTGATCTCCAGGATGGATAACGCCACGTACAGCTACGTGAACGACTGCACTCCGCTCACCAACTGTAAGGCCGGCCGTAAGGCTGGGGGCTCCACGGTGGTCAACATGAGCCATGCGGGCAAGAAGCCACCAAATGACTACCTGGTCTGGTCACTCTGCAACACTCTTTACGTCAACTTCTGCTGCCTGGGCTTCATGGCTCTGATCTACTCCATCAAGGTGAGCTTCAGCACCGATATCTTCATCTGTTTTGTG',
        'gRNA_start': 147,
        'gRNA_end': 170,
        'PAM_start': 147,
        'AA_change_site': 168,
        'which_PAM' : 'CCN'
    }

    #'NGG' PAM
    full_seq = 'atatttagaaattggaagataatacaattaaattcaagcaattcaatattaaacctaaaaaattcaactacatttttccattgttttgcttctcttgatttttttttctattttaaatttgtatttaatatttttctataacataattttgggtgtacttgtttttggactgttattttgttaaataagctccagatttggcttcagtactgactaatctaaagtacaaatataatattgtatgtcttttaaaaatatgaatttaaagatagatttgtgaggggtttacttgcatatactgagcacagtatatacagcttgttattgatgtaaaatgcaggggtatcccactaaataacatcatgagtatttccagaattatttttccatatgaacgtttttatcttggcacacatgcattgacagactgtttggatttatgtaataaaaatggtttccaattctagagtaggaacaatttccacctgcctcctgtaattattggctgatctttgcattttgtgtgacagATGGCTGCGATTCGTAAGAAGCTTGTAATCGTTGGAGATGGAGCCTGTGGAAAGACGTGTTTACTCATCGTTTTCAGTAAAGACCAGTTTCCTGAGGTTTACGTTCCCACAGTGTTTGAGAACTATGTCGCTGATATTGAGGTGGATTCAAAACAGgtgagattatctcttaaacggccatcccacaatcaccatctctctgtaacaatgtgaacttttgcaaaaagtaacttcctcagattgttcacaaacattttaatagtgttcttccttcttttttgtattttttatgcatggtttatcattaatgaatatgttttattgtaattgtgctaacaacattgtgtgtatctgtgaaatatgagtaactaaaaccctcaaaaagccatttgatattttcttctaaaatttgctcctgttggtaggctcaattattttacatttatagtgacggaaaggttttattcattgcctttcaagtgaaataactatacataaccaaagttggatgagaaaaatttgaaaatgtgaaacaggctgtgaaaaatgtgaaaagtttttttttcctccaatttctgtttaaaggagagaactcatttcttaacataatagttttaataattcatctctaataacggatttatttaatctttgccatgatgacagtatataatatttgactagatatttttcaagacacttccatatagcttaaagggcaacttataggcttaactaggttaattaggttaattaggcaggttagggtaactaggcaagttgttgtacaatgatggtttgttg'

    seq_dict = tiarsdesigner.utils.get_endogenous_sequence_and_stats(
        full_seq,
        557,
        580,
        577,
        578,
        'NGG'
    )
    assert seq_dict == {
        'sequence': 'TTCCAATTCTAGAGTAGGAACAATTTCCACCTGCCTCCTGTAATTATTGGCTGATCTTTGCATTTTGTGTGACAGATGGCTGCGATTCGTAAGAAGCTTGTAATCGTTGGAGATGGAGCCTGTGGAAAGACGTGTTTACTCATCGTTTTCAGTAAAGACCAGTTTCCTGAGGTTTACGTTCCCACAGTGTTTGAGAACTATGTCGCTGATATTGAGGTGGATTCAAAACAGGTGAGATTATCTCTTAAACGGCCATCCCACAATCACCATCTC',
        'gRNA_start': 102,
        'gRNA_end': 125,
        'PAM_start': 122,
        'AA_change_site': 123,
        'which_PAM': 'NGG'
    }

    # test get_restriction_sites
    BsrI = get_enzyme('BsrI')
    BseNI = get_enzyme('BseNI')
    Bse1I = get_enzyme('Bse1I')
    BmeDI = get_enzyme('BmeDI')
    RlaI = get_enzyme('RlaI')
    TaqI= get_enzyme('TaqI')
    Hpy188I = get_enzyme('Hpy188I')

    assert tiarsdesigner.utils.get_restriction_sites('ACTGGTT', ['BsrI', 'RlaI']) == {
        BsrI: [7],
        RlaI: [1,3],
    }
    # reverse compliment of recognition site in sequence
    assert tiarsdesigner.utils.get_restriction_sites('AACCAGT', ['BsrI', 'RlaI']) == {
        BsrI: [4],
        RlaI: [3,5]
    }
    # no enzyme list provided
    assert tiarsdesigner.utils.get_restriction_sites('ACTGGTT') == {
        RlaI: [1,3],
        BsrI:[7],
        BmeDI: [4,5,5],
        BseNI: [7],
        Bse1I: [7],
    }

    # test get_Rsite_dict
    assert tiarsdesigner.utils.get_Rsite_dict('ACTGGTT', [BsrI, RlaI]) == {
        BsrI: [(0,5)],
        RlaI: [(0,3),(2,5)],
    }
    # also finds reverse compliments of the sites
    assert tiarsdesigner.utils.get_Rsite_dict('AACCAGT', [BsrI, RlaI]) == {
        BsrI:[(2,7)],
        RlaI: [(2,5),(4,7)],
    }
    # give it a dictionary of Renzymes instead of a list
    assert tiarsdesigner.utils.get_Rsite_dict('ACTGGTT', {BsrI: [7], RlaI: [1,3]}) == {
        BsrI: [(0,5)],
        RlaI: [(0,3),(2,5)],
    }

    # test reformat codon table
    assert tiarsdesigner.utils.reformat_codon_table({'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M'}) =={
        'ATT': ['ATT', 'ATC', 'ATA'],
        'ATC': ['ATT', 'ATC', 'ATA'],
        'ATA': ['ATT', 'ATC', 'ATA'],
        'ATG': ['ATG'],
        'TAA': ['TAA'],
        'TGA': ['TGA'],
        'TAG': ['TAG']}


    # test synonymous_seqs
    assert list(tiarsdesigner.utils.synonymous_seqs('ATCCAATGG')) ==[
        'ATTCAATGG', 'ATTCAGTGG', 'ATCCAATGG', 'ATCCAGTGG', 'ATACAATGG', 'ATACAGTGG'
    ]

    # test map_recognition_sites
    assert tiarsdesigner.utils.map_recognition_sites('AGCTCGATG', TaqI, 1) ==[
            0, 0, 0, 1, 1, 1, 1, 0, 0
    ]

    # test get_substitutions_Renzyme_depends_on
    assert tiarsdesigner.utils.get_substitutions_Renzyme_depends_on(Hpy188I, 'TTCTGCAGCCCTGGATGGGCTCTGCTGCTA', 'TTCTTCAGAGTTGGATGGGTTCTGCTGCTA',[(4,9)])== [(4,8)]

    assert tiarsdesigner.utils.get_codon_frequency_dict('codon_frequency.tsv')=={
            'TTT':18.2,
            'TCT':16.9,
            'TAT':12.6,
            'TGT':11.3,
            'TTC':20.8,
            'TCC':15.2,
            'TAC':17.0,
            'TGC':11.2,
            'TTA':7.0,
            'TCA':13.2,
            'TAA':1.1,
            'TGA':1.4,
            'TTG':12.3,
            'TCG':5.6,
            'TAG':0.6,
            'TGG':11.6,
            'CTT':12.7,
            'CCT':16.6,
            'CAT':10.9,
            'CGT':6.9,
            'CTC':17.0,
            'CCC':12.7,
            'CAC':14.8,
            'CGC':9.6,
            'CTA':6.2,
            'CCA':15.7,
            'CAA':11.8,
            'CGA':6.7,
            'CTG':37.6,
            'CCG':8.2,
            'CAG':33.5,
            'CGG':6.6,
            'ATT':16.5,
            'ACT':14.5,
            'AAT':16.3,
            'AGT':13.2,
            'ATC':23.7,
            'ACC':16.2,
            'AAC':24.1,
            'AGC':18.4,
            'ATA':7.7,
            'ACA':17.0,
            'AAA':29.3,
            'AGA':14.3,
            'ATG':25.5,
            'ACG':7.4,
            'AAG':30.7,
            'AGG':10.2,
            'GTT':14.1,
            'GCT':20.9,
            'GAT':24.8,
            'GGT':13.7,
            'GTC':14.8,
            'GCC':19.5,
            'GAC':27.8,
            'GGC':17.2,
            'GTA':6.7,
            'GCA':16.6,
            'GAA':24.4,
            'GGA':21.5,
            'GTG':28.3,
            'GCG':8.6,
            'GAG':42.8,
            'GGG':10.0,
    }
    
















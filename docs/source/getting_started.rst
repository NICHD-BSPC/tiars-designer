Getting Started
===============

This section describes how to install TIARS Designer and its dependencies.


Installing TIARS Designer
-------------------------

With Python already installed, to use the published version from PyPI:

.. code-block:: bash

    python -m pip install tiarsdesigner

Or to install an editable version from a local copy of the repository, from top level of the repo run:

.. code-block:: bash

   python -m pip install -e .

To test whether installation was successful, run the tool without any inputs. If the installation was successful you should see the following usage message:

.. code-block:: console

   $ tiarsdesigner -h
   usage: tiarsdesigner [-h] [--endogenous_seq ENDOGENOUS_SEQ] [--gRNA_start GRNA_START] [--gRNA_end GRNA_END] [--PAM_start PAM_START] [--which_PAM WHICH_PAM]
                     [--AA_change_site AA_CHANGE_SITE] [--new_codon NEW_CODON] [--locus_name LOCUS_NAME]
                     [--single-enzyme-filters SINGLE_ENZYME_FILTERS [SINGLE_ENZYME_FILTERS ...]] [--enzyme-pair-filters ENZYME_PAIR_FILTERS [ENZYME_PAIR_FILTERS ...]]
                     [--codon-usage-tsv CODON_USAGE_TSV] [--min-Rsite-distance MIN_RSITE_DISTANCE] [--enzymes ENZYMES] [--verbose]

   options:
     -h, --help            show this help message and exit
     --endogenous_seq ENDOGENOUS_SEQ
     --gRNA_start GRNA_START
     --gRNA_end GRNA_END
     --PAM_start PAM_START
     --which_PAM WHICH_PAM
     --AA_change_site AA_CHANGE_SITE
     --new_codon NEW_CODON
     --locus_name LOCUS_NAME
     --single-enzyme-filters SINGLE_ENZYME_FILTERS [SINGLE_ENZYME_FILTERS ...]
     --enzyme-pair-filters ENZYME_PAIR_FILTERS [ENZYME_PAIR_FILTERS ...]
     --codon-usage-tsv CODON_USAGE_TSV
     --min-Rsite-distance MIN_RSITE_DISTANCE
     --enzymes ENZYMES     name of file containing single column list of enzymes to use
     --verbose             if True, this will write csvs/tsvs for each step
    usage: tiarsdesigner [-h] [--endogenous_seq ENDOGENOUS_SEQ]
                     [--gRNA_start GRNA_START] [--gRNA_end GRNA_END]
                     [--PAM_start PAM_START] [--which_PAM WHICH_PAM]
                     [--AA_change_site AA_CHANGE_SITE] [--new_codon NEW_CODON]
                     [--locus_name LOCUS_NAME]
                     [--single-enzyme-filters SINGLE_ENZYME_FILTERS [SINGLE_ENZYME_FILTERS ...]]
                     [--enzyme-pair-filters ENZYME_PAIR_FILTERS [ENZYME_PAIR_FILTERS ...]]
                     [--codon-usage-tsv CODON_USAGE_TSV]
                     [--min-Rsite-distance MIN_RSITE_DISTANCE]
                     [--enzymes ENZYMES] [--verbose]

    options:
        -h, --help      show this help message and exit
    --endogenous_seq ENDOGENOUS_SEQ
                        The full endogenous sequence to be edited. Must extend
                        at least 85 nucleotides upstream (5') of the Cas
                        protein cut and 50 nucleotides downstream (3') of the
                        Cas cut.
    --gRNA_start GRNA_START
                        The zero-indexed location in the sequence of the
                        beginning of the guide RNA (in this case the PAM site
                        is considered to be part of the guide RNA)
    --gRNA_end GRNA_END 
                        The zero-indexed location in the sequence of the end
                        of the guide RNA (again, PAM is considered part of the
                        guideRNA). This should correspond to the location of
                        the first bp after the last bp of the guide RNA, such
                        that slicing the string in python (e.g.
                        endogenous_seq[gRNA_start: gRNA_end]) would yield the
                        full gRNA )
    --PAM_start PAM_START
                        The zero-indexed location in the sequence of the
                        beginning of the PAM site
    --which_PAM WHICH_PAM
                        Specifies which type of PAM your experiment uses.
                        Right now Cas9 (``CCN`` and ``NGG``) and Cpf1
                        (``TVVV`` and ``BAAA``) are supported. NOTE: put the
                        generic sequence of the PAM rather than the exact one
                        the sample uses (e.g. put ``NGG`` rather than ``AGG``)
    --AA_change_site AA_CHANGE_SITE
                        The zero-indexed location in the sequence of the codon
                        that will contain desired CRISPR edit (even if the
                        actual edit changes the second or third basepair in
                        the codon, this value should still correspond to the
                        first basepair in the codon)
    --new_codon NEW_CODON
                        The basepair sequence of the new codon containing the
                        edit (for example, if you wanted to change a GGA
                        glycine to an AGA arginine, you would put AGA)
    --locus_name LOCUS_NAME
                        The name of the locus/sample, which will be used to
                        name the output folder and files
    --single-enzyme-filters SINGLE_ENZYME_FILTERS [SINGLE_ENZYME_FILTERS ...]
                        " A list of single enzyme filters to run. To run all
                        of the single-enzyme filters (recommended), input
                        'all'. The filters names that can be used are: -
                        ``no_sandwiched_cutsites`` -
                        ``no_endogenous_cutsites_within_Nnt_of_novel_cutsite``
                        - ``novel_site_depends_on_at_least_two_substitutions``
                        - ``num_recognition_sites_equals_num_cutsites`` -
                        ``no_unneeded_substitutions``
    --enzyme-pair-filters ENZYME_PAIR_FILTERS [ENZYME_PAIR_FILTERS ...]
                        A list of single enzyme filters to run. To run all of
                        the single-enzyme filters (recommended), input 'all'.
                        The filters names that can be used are: -
                        ``no_sandwiched_cutsites`` -
                        ``no_endogenous_cutsites_within_Nnt_of_novel_cutsite``
                        - ``novel_site_depends_on_at_least_two_substitutions``
                        - ``num_recognition_sites_equals_num_cutsites`` -
                        ``no_unneeded_substitutions``
    --codon-usage-tsv CODON_USAGE_TSV
                        A tsv of codon usage frequency values for the organism
                        used. The tsv should have one codon per line in the
                        format: [codon]>[frequency] (e.g. ``UGA>1.4``) where
                        '>' represents a tab. The tsv in the github repository
                        is for zebrafish. Organism specific codon frequencies
                        can be found at the `Kazusa Codon usage Database
                        <http://kazusa.or.jp/codon>`_
    --min-Rsite-distance MIN_RSITE_DISTANCE
                        The minimum allowable distance (in basepairs) between
                        the novel restriction site and any endogenous
                        restriction sites for a given enzyme. This is to avoid
                        having fragment sizes that are too similar to be able
                        to tell the difference between a novel site cut and
                        and endogenous site cut. This is the value used in the
                        'no Endogenous Cutsites within N nt of Novel Site'
                        filter. (Default: 5)
    --enzymes ENZYMES   
                        A tsv/txt file containing a list of restriction
                        enzymes to use (must be in the Bio.Restriction
                        database, case sensitive) The list should have one
                        enzyme per line. (Default: Bio.Restriction
                        'CommonOnly' list of commonly used enzymes).
    --verbose           If True, this will write out tsv files with passing
                        percentages for each filter as well as raw the total
                        number of candidates before and after filtering.

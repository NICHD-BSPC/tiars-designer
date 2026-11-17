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

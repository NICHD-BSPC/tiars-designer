Usage
=====

This page explains the command line arguments TIARS designer takes to run and a usage example.

Inputs
------

Once TIARS-Designer is installed, the following arguments are required to run the tool.

.. _required-inputs:


Required Inputs
~~~~~~~~~~~~~~~

.. table::
   :widths: 30 70

   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | Argument                | Description                                                                                                        |
   +=========================+====================================================================================================================+
   | endogenous_seq          | The full endogenous sequence to be edited. Must extend at least 85 nucleotides upstream(5') of the Cas protein cut |
   |                         | and 50 nucleotides downstream (3') of the Cas cut.                                                                 |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | gRNA_start              | The location in the sequence* of the beginning of the guide RNA (in this case the PAM site is considered to be part|
   |                         | of the guide RNA)                                                                                                  |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | gRNA_end                | The location in the sequence* of the end of the guide RNA (again, PAM is considered part of the guideRNA). This    |
   |                         | should correspond to the location of the first bp after the last bp of the guide RNA, such that slicing the string |
   |                         | in python (e.g. endogenous_seq[gRNA_start: gRNA_end]) would yield the full gRNA                                    |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | PAM_start               | The location in the sequence* of the beginning of the PAM site                                                     |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | which_PAM               | Specifies which type of PAM your experiment uses. Right now Cas9 (CCN and NGG) and Cpf1 (TVVV and BAAA) are        |
   |                         | supported. NOTE: put the generic sequence of the PAM rather than the exact one the sample uses (e.g. put NGG rather|
   |                         | than AGG)                                                                                                          |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | AA_change_site          | The location in the sequence* of the codon that will contain desired CRISPR edit (even if the actual edit changes  |
   |                         | the second or third basepair in the codon, this value should still correspond to the first basepair in the codon)  |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | new codon               | The basepair sequence of the new codon containing the edit (for example, if you wanted to change a GGA glycine to  |
   |                         | an AGA arginine, you would put AGA)                                                                                |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | sample_name             | The name of the locus/sample, which will be used to name the output folder and files                               |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | \--codon-usage-tsv      | A tsv of codon usage frequency values for the organism used. The tsv should have one codon per line in the format: |
   |                         | [codon]>[frequency] (e.g. UGA>1.4) where '>' represents a tab. The tsv in the github repository is for zebrafish.  | 
   |                         | Organism specific codon frequencies can be found at the `Kazusa Codon usage Database <http://kazusa.or.jp/codon>`_ |           
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | \--single-enzyme-filters| A list of single enzyme filters to run. To run all of the single-enzyme filters (recommended), input 'all'.        |
   |                         | The filters names that can be used are:                                                                            |
   |                         |                                                                                                                    |
   |                         | -  no_sandwiched_cutsites                                                                                          |
   |                         | -  no_endogenous_cutsites_within_Nnt_of_novel_cutsite                                                              |
   |                         | -  novel_site_depends_on_at_least_two_substitutions                                                                |
   |                         | -  num_recognition_sites_equals_num_cutsites                                                                       |
   |                         | -  no_unneeded_substitutions                                                                                       |
   |                         |                                                                                                                    |
   |                         | See :ref:`Filters` for a complete description of the filters.                                                      |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+
   | \--enzyme-pair-filters  | A list of enzyme pair filters to run. To run all of the enzyme-pair filters(recommended), input 'all'.             |
   |                         | The filter names that can be used are:                                                                             |
   |                         |                                                                                                                    | 
   |                         | * at_least_one_has_endogenous_site                                                                                 |
   |                         | * both_have_endogenous_sites                                                                                       |
   |                         | * compatible_polarity                                                                                              |
   |                         | * no_unneeded_substitutions                                                                                        | 
   |                         | * both_Renzymes_depend_on_at_least_one_unique_substitution                                                         |
   |                         |                                                                                                                    |
   |                         | See :ref:`Filters` for a complete description of the filters.                                                      |
   +-------------------------+--------------------------------------------------------------------------------------------------------------------+


* NOTE: any time the 'location in the sequence' of an element is mentioned, this is referring to its coordinate in the sequence, where the first basepair in the sequence is considered 0, the second basepair is 1, the third is 2, etc. 


Optional Inputs
~~~~~~~~~~~~~~~

The tool has the following optional arguments to customize use: 

.. table::
   :widths: 30 70

   +----------------------+---------------------------------------------------------------------------------------------------------------------------+
   | Argument             | Description                                                                                                               |
   +======================+===========================================================================================================================+
   | \--min_Rsite_distance| The minimum allowable distance (in basepairs) between the novel restriction site and any endogenous restriction sites     |
   |                      | for a given enzyme. This is to avoid having fragment sizes that are too similar to be able to tell the difference between |
   |                      | a novel site cut and and endogenous site cut. This is the value used in the 'no Endogenous Cutsites within N nt of Novel  |
   |                      | Site' filter (see :ref: `Filters`). Default=5                                                                             |
   +----------------------+---------------------------------------------------------------------------------------------------------------------------+
   | \--enzymes           | A tsv/txt file containing a list of restriction enzymes to use (must be in the Bio.Restriction database, case sensitive). |
   |                      | The list should have one enzyme per line                                                                                  |
   |                      | Default= Bio.Restriction 'CommonOnly' list of commonly used enzymes.                                                      |
   +----------------------+---------------------------------------------------------------------------------------------------------------------------+
   | \--verbose           | If this flag is included, output tsv files with passing percentages for each filter as well as raw the total number of    |
   |                      | candidates before and after filtering.                                                                                    |
   +----------------------+---------------------------------------------------------------------------------------------------------------------------+

Example Usage
-------------

The following example shows a sequence with its coordinates mapped out for easy viewing of the PAM_start, and gRNA start and end sites (in reality this sequence would have to be longer to go in the tool, it's been shortened for easier display). 

.. code-block:: rest

   # coordinate 0         1         2         3         4         5
   #            012345678901234567890123456789012345678901234567890

   # sequence   GCTGGTGGCCCCCAAATATTCAGGTTTAATGTTCACTCTGCAACACTCTTT
   #                       |                    | | 
   #               gRNA_start and       AA_change gRNA_end
   #                  PAM_start

For this example, the codon to be changed is the TCA starting at coordinate 32, and we want to change it to TTA. The gRNA sequence is CCAAATATTCAGGTTTAATGTTC. This example would have the following values:

| sample: example_sequence
| AA_change_site: 32
| new_codon: TTA
| PAM_start: 11
| gRNA_start: 11
| gRNA_end: 34
| which_PAM: CCN
| sequence: GCTGGTGGCCCCCAAATATTCAGGTTTAATGTTCACTCTGCAACACTCTTT


An example command to run this sample would be:

.. code-block:: console

   tiarsdesigner \
      --endogenous_seq=GCTGGTGGCCCCCAAATATTCAGGTTTAATGTTCACTCTGCAACACTCTTT   \
      --gRNA_start=11 \
      --gRNA_end=34  \
      --PAM_start=11 \
      --which_PAM=CCN \
      --AA_change_site=32 \
      --new_codon=TTA \
      --locus_name=example_sequence \
      --single-enzyme-filters=all \
      --enzyme-pair-filters both_have_endogenous_sites compatible_polarity \
      --codon-usage-tsv /path/to/codon_usage.tsv \
      --enzymes /path/to/enzymes.tsv \
      --verbose


Running Tiars Designer without Endogenous Sites
-----------------------------------------------

For uses cases where pseudo-endogenous sites can be incorporated into the primer and endogenous recognition sites in the sequence are not needed, TIARS Designer can be run without including them. This can be done by removing the 'at least one has endogenous site', 'both have endogenous sites',  and 'compatible polarity' filters (AKA run it with all enzyme pair filters EXCEPT for those three). This will stop the tool from requiring the restriction enzymes to have an endogenous recognition site in order to make it to the output.

In order to only output options that don't have endogenous sites, the optional '--min-Rsite-Distance' parameter can be set to a number longer than the length of the sequence. This effectively filters out all options with endogenous recognition sites. 

The usage example below demonstrates how to run the tool in this way. 


.. code-block:: console

   tiarsdesigner \
     --endogenous_seq=GCTGGTGGCCCCCAAATATTCAGGTTTAATGTTCACTCTGCAACACTCTTT \
     --gRNA_start=11 \
     --gRNA_end=34 \
     --PAM_start=11 \
     --whic_PAM=CCN \
     --AA_change_site=32 \
     --new_codon=TTA \
     --locus_name=example_sequence \
     --single-enzyme-filters=all \
     --enzyme-pair-filters no_unneeded_substitutions both_Renzymes_depend_on_at_least_one_unique_substitution \
     --codon-usage-tsv /path/to/codon_usage.tsv  \
     --enzymes /path/to/enzymes.tsv \
     --verbose \
     --min-Rsite-distance=500


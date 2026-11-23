Usage
=====

This page contains usage examples with the sequence and guide RNA out for visualization

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


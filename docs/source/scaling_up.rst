Scaling up with Snakemake
=========================

If you have multiple samples to run, you can run them all at once using the provided Snakefile. This will require consolidating the sample information into a 'sampletable.tsv' file with a 'config.yaml' file for the less sample specific arguments. Template versions of these files can be found in the 'scaling up' folder of this repository. Once the sampletable and config file are set up, you can run the Snakefile with the following command:

.. code-block:: 

    snakemake --cores N


Sampletable.tsv
---------------

For the sampletable.tsv file, each row will represent a sample and will have the following columns. Descriptions for each of these columns can be found on the :ref:`Required Inputs <required-inputs>` page. 

- sample
- AA_change_site
- new_codon
- PAM_start
- gRNA_start
- gRNA_end
- which_PAM
- sequence

Example sampletable.tsv
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    sample	AA_change_site	new_codon	PAM_start	gRNA_start	gRNA_end	which_PAM	sequence
    sample_A	69	TTA	48	48	71	CCN	GCTGGGGGCTCCACGGTGGTCAACATGAGCCATGCGGGCTGGTGGCCCCCAAATTGGTGGTGGTGGTGGTCACTCTGCAACACTCTTTACGTCAACTTCTGC
    sample_B	66	AGA	66	46	69	NGG	TTGCATTTTGTGTGACAGATGGCTGCGATTCGTAAGAAGCTTGTAATCGTTGGATGGTGGTGGTGGGGGTGGTGGTGGTTACTCATCGTTTTCAGTAAAGACC
    sample_C	66	ACC	51	51	76	TTTV	AGTGGGAAACAGAGACCGAACGGCACTGAGTTCATGCACTTCCTCCCCATGTTTCTCTGGTCTTGGCCCTGGTGGTGGTGTTGGTGGGTGTGGGTTTAAACAT

Config.yaml
-----------

The config.yaml file will include the following arguments. Descriptions can be found on the :ref:`Required Inputs <required-inputs>` page.

- sampletable: path to 'sampletable.tsv file'
- enzyme list: path to tsv file of desired enzymes. If default full enzyme list is desired write 'default'
- codon_usage_tsv: path to tsv file of codon usage frequencies for organism used
- min_Rsite_distance: put '5' if you want to use the default
- single_enzyme_filters: list of desired single enzyme filters. To use all filters put 'all'
- enzyme_pair_filters: list of desired enzyme pair filters. To use all filters put 'all'
- verbose: boolean value

Example config.yaml
~~~~~~~~~~~~~~~~~~~

.. code-block:: rest

    sampletable: 'sampletable.tsv'

    enzyme_list: 'default'

    codon_usage_tsv: 'codon_frequency.tsv' 

    min_Rsite_distance: 5

    single_enzyme_filters: 'all'

    enzyme_pair_filters: 
     - 'compatible_polarity'
     - 'both_Renzymes_depend_on_at_least_one_unique_substitution'
     - 'no_unneeded_substitutions'

    verbose: TRUE



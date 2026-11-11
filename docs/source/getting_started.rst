Getting Started
===============

This section describes how to install TIARS Designer and its dependencies.

Dependencies
------------

Use `conda <https://docs.conda.io/en/latest/>`_ to create an environment and install the dependencies listed in the 'requirements.txt' file in the top level of the github repository.

.. code-block:: bash

    conda create -p ./env --file requirements.txt

Installing TIARS Designer
-------------------------

With the conda environmnent activated, install TIARS designer using pip

.. code-block:: bash

    conda activate ./env
    pip install -i https://test.pypi.org/simple/ tiarsdesigner

To test whether installation was successful, run the tool without any inputs. If the installlation was successful you should see the following usage message:

.. code-block:: console

   $tiarsdesigner
   usage: tiarsdesigner [-h] [--single-enzyme-filters SINGLE_ENZYME_FILTERS [SINGLE_ENZYME_FILTERS ...]] [--enzyme-pair-filters ENZYME_PAIR_FILTERS [ENZYME_PAIR_FILTERS ...]]
                     [--codon-usage-tsv CODON_USAGE_TSV] [--min-Rsite-distance MIN_RSITE_DISTANCE] [--enzymes ENZYMES] [--verbose]
                     endogenous_seq gRNA_start gRNA_end PAM_start which_PAM AA_change_site new_codon locus_name




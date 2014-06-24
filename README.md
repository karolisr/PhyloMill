# PhyloMill

<!--Available on [GitHub](https://github.com/karolisr/krpy)-->

Using sequence alignments of multiple loci to build phylogenies is a common practice. Aligning sequences and reconstructing phylogeny are mostly automated steps, however the acquisition and filtering of raw sequence data is tedious. Often sequences at public databases are mislabeled in some way. Organism names may be inaccurate, sequence direction may be misannotated, or worse, sequences may be labeled entirely incorrect. Even when none of these problems occur, the increasing number of sequences available at public databases makes it very difficult to manually download and inspect all of them.

PhyloMill automates the process of multi-locus dataset construction so the researcher can spend more time focusing on the biological questions and not the frustrating technical details.

## How it works

Basic use of the pipeline requires little configuration. Only thing you will need is a set of taxon and locus names that you want to use.

It is not necessary to list all species names you want to use. If you  wanted to cover all species in the genus *Solanum*, you could simply put *Solanum* in the configuration file under [Main Taxa] heading. If you want to exclude some species, simply put them under [Excluded Taxa] heading. You can explicitly specify outgroups under [Outgroup Taxa] heading.

To tell the pipeline which loci you want to use in the study simply put any of the predefined locus names under the [Loci] heading in the configuration file. Currently predefined loci are:

ATP6, ATP8, atpB, CCR2, COI, cytb, EGR1, FGB, GBSSI, GGR, ITS1-5.8S-ITS2, matK, MYC, NADH2, ndhF, RAG1, RAG2, rbcL, TGFB2, trnS-trnG

The locus definitions are simple text files and can be found in the "search_strategies" directory in your project. In case you want to use some other locus, you can simply edit and save one of these text files with a new name. I will write a short description of how to do this soon.

To run the pipeline first you would run the command listed below, which will create a new project directory:

    phylomill -p PATH_TO_PROJECT_DIR

This clean project directory contains a sample configuration file that will allow to construct a phylogeny of a select taxa from *Cetacea* using three loci: COI, RAG1, cytb.

Now you can simply run:

    phylomill -p PATH_TO_PROJECT_DIR -c autopilot

This process takes about 7 minutes on a fairly modern MacbookPro.

## What it does

1. PhyloMill will download records for the taxa and loci of interest from Genbank.

2. Optionally, organism names may be checked:

    * It is possible to drop hybrid species and unknown species, e.g. "*Solanum sp.*".
    * It is possible to predefine simple species-to-species mappings or GI-to-species mappings if you know an accession is mislabeled.
    * It also possible to define a more complex synonymy table, that I will document later.

3. Loci of interest will be extracted from the downloaded records, as these may contain several loci. At this point PhyloMill will try to identify and blacklist misannotated records under the assumption that majority of the downloaded records are annotated correctly.

4. There may be several sequences per locus, per organism. At this step, all the sequences for the locus at hand from the particular species will be aligned together using [MAFFT](http://mafft.cbrc.jp/alignment/software) by default. [Muscle](http://www.drive5.com/muscle) and [Clustal Omega](http://www.clustal.org/omega) are also available. If non-overlapping set of sequences for the same locus is detected (as it may happen when several different exons are sequenced, etc.), PhyloMill will use longest sequences in the dataset to find their position.

5. The alignments from the previous step will be saved in files under "output/01_locus_files" directory. These can be reviewed for accuracy and edited manually. During the next run PhyloMill will incorporate these changes into its database.

6. The alignments from the previous step will be collapsed into consensus sequences and will be used for between-species locus alignments.

7. Alignments for each locus will be concatenated. All alignments produced can be edited manually, PhyloMill will incorporate these changes.

8. Species tree (using partitioned concatenated alignment) and, optionally, gene trees will be constructed using [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html). You can use your own favorite phylogeny program at this point, simply point it to the alignments in the "output/02_alignments" directory. Trees are saved in the "output/03_trees" directory.

## Installation

Currently (as of June 24th, 2014) I am working on the detailed installation instructions. This should be a relatively simple process on MacOSX and Linux, however, PhyloMill depends on a few other Python packages and a couple of external alignment/phylogeny programs which may confuse the matter. I will make installation instructions available as soon as possible.

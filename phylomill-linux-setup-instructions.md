# PhyloMill Setup on a GNU Linux System

PhyloMill is a Python 2.7 script that pools several other tools together to
automatically produce multi-locus alignments for a set of taxa. Because of this,
all the pieces of software that PhyloMill depends on need to be present for the
pipeline to work.

This document will explain how to install PhyloMill pipeline and all the
dependencies needed for PhyloMill to function properly on a relatively modern
GNU Linux distribution. Instructions below specifically address installation on
Ubuntu, but it should be trivial to modify the commands for other distributions.
You can skip over things that you know are already installed on your system.

1. Install dependencies

    First, let's update Ubuntu repositories:

        sudo apt-get -y update

    Git

        sudo apt-get -y install git

    This will install [Biopython](http://biopython.org). First, let's install
    the dependencies:

        sudo apt-get -y build-dep python-biopython

    It is possible to install Biopython straight from Ubuntu repositories:

        sudo apt-get -y install python-biopython

    However, depending on your Ubuntu version, this may be hoplessly outdated.
    Let's install straight from their website:

        wget http://biopython.org/DIST/biopython-1.64.tar.gz
        tar xzf biopython-1.64.tar.gz
        cd biopython-1.64
        python setup.py build
        python setup.py test

    The previous step will take a bit...

        sudo python setup.py install
        cd ..
        rm -rf biopython-1.64*

    This will install [SciPy](http://www.scipy.org).

        sudo apt-get -y install python-scipy

    These commands will install a few Python libraries required by PhyloMill

        sudo apt-get -y install python-Levenshtein
        sudo apt-get -y install python-unidecode
        sudo apt-get -y install python-bs4
        sudo apt-get -y install python-bs4-doc

2. Install [MAFFT](http://mafft.cbrc.jp/alignment/software)

    Follow the commands below to install the latest version of MAFFT:

        wget http://mafft.cbrc.jp/alignment/software/mafft-7.158-with-extensions-src.tgz
        tar xzf mafft-7.158-with-extensions-src.tgz
        cd mafft-7.158-with-extensions/core
        make clean
        make
        sudo make install
        cd ..
        cd extensions
        make clean
        make
        sudo make install
        cd ..
        cd ..
        rm -rf mafft-7.158-with-extensions*

3. Install [Muscle](http://www.drive5.com/muscle)

    Follow the commands below to install the latest version of Muscle:

        wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
        tar xzf muscle3.8.31_i86linux64.tar.gz
        sudo mv muscle3.8.31_i86linux64 /usr/local/bin/muscle
        rm -rf muscle3.8.31_i86linux64.tar.gz

4. Install [RAxML](https://github.com/stamatak/standard-RAxML)

    Follow the commands below to install the latest version of RAxML:

        git clone https://github.com/stamatak/standard-RAxML.git
        cd standard-RAxML
        make -f Makefile.AVX.PTHREADS.gcc

    Depending on your system the previous step may fail, then try this:

        rm *.o
        make -f Makefile.PTHREADS.gcc

    Then, if the first version worked:

        sudo mv raxmlHPC-PTHREADS-AVX /usr/local/bin/raxml

    Else:

        sudo mv raxmlHPC-PTHREADS /usr/local/bin/raxml

    Clean up:

        cd ..
        rm -rf standard-RAxML

5. Install BLAST+ command line tools

        sudo apt-get -y install ncbi-blast+

6. At this point, we can install PhyloMill

    Type these commands to download and setup PhyloMill:

        cd ~
        git clone https://github.com/karolisr/krpy.git

    To update PhyloMill at any point you can do this:

        cd ~/krpy
        git pull

7. Run PhyloMill by typing:

        cd ~/krpy/krpy/workflows
        ./phylomill -p ~/pm_test_project -c autopilot

    This will create a sample project with name "pm_test_project" in your home
    directory and produce a sample phylogeny.

<!--
        cd ~
        git clone https://github.com/karolisr/krpy.git

        echo -e "\n# PhyloMill" | tee -a ~/.bashrc
        phylomill_bin_path=PATH=\"${HOME}/krpy/krpy/workflows:'${PATH}'\"
        echo -e $phylomill_bin_path | tee -a ~/.bashrc
        echo -e "export PATH" | tee -a ~/.bashrc

        phylomill_python_path=PYTHONPATH=\"${HOME}/krpy:'${PYTHONPATH}'\"
        echo -e $phylomill_python_path | tee -a ~/.bashrc
        echo -e "export PYTHONPATH" | tee -a ~/.bashrc

    Restart the terminal and you are done.

    To update PhyloMill at any point you can do this:

        cd ~/krpy
        git pull
-->

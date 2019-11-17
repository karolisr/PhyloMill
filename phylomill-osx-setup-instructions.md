# PhyloMill Setup on a completely clean OS X 10.9 (Mavericks) system

PhyloMill is a Python 2.7 script that pools several other tools together to
automatically produce multi-locus alignments for a set of taxa. Because of this,
all the pieces of software that PhyloMill depends on need to be present for the
pipeline to work.

This document will explain how to install PhyloMill pipeline and all the
dependencies needed for PhyloMill to function properly on a clean OS X 10.9
install. You can skip over things that you know are already installed on
your system.

1. Latest Apple updates

    Click on the Apple icon in the top-left corner and select
    "Software Update...". On the screen that comes up click "Update All" if any
    updates are available.

2. Xcode and command line developer tools

    Go to the "App Store", find and install Xcode. This will take a bit. Xcode
    provides basic tools required to develop for the Mac and is required to get
    many PhyloMill dependencies to work. If you had done any bioinformatics work
    on your system, it is very likely you already have Xcode, if not grab a cup
    of coffee this will take a while.

    Once Xcode is installed, you will be able to install command line developer
    tools. Open the Terminal app and type in this:

        xcode-select --install

    A window will pop up, click "Install", then click "Agree". Once that is
    done, type in this:

        sudo xcodebuild -license

    You will be asked for your password, type it in and press enter. A license
    text will pop up, press 'q' to quit. Type the word "agree" and press enter.
    That was easy.

3. MacPorts

    MacPorts provides a centralized way of managing a lot of popular unix/linux
    software that is not available at the App Store. Go to this
    [website](https://www.macports.org/install.php) and download MacPorts
    installation file. Double click the file and install it. Should take about a
    minute.

    Open the Terminal app again and type this in:

        export PATH="/opt/local/bin:/opt/local/sbin:$PATH"

    Then update the list of available MacPorts software, type in this (you will
    be asked for your password):

        sudo port -d selfupdate

4. Install dependencies using MacPorts

    This will install a newer version of [git](http://git-scm.com) than is
    available on OS X. In the Terminal app run these lines:

        sudo port install git

    This will install version of Python required by PhyloMill.

        sudo port install python27
        sudo port select --set python python27

    This will install [Biopython](http://biopython.org).

        sudo port install py27-biopython

    This will install [SciPy](http://www.scipy.org).

        sudo port install py27-scipy

    These commands will install a few Python libraries required by PhyloMill

        sudo port install py27-datrie
        sudo port install py27-levenshtein
        sudo port install py27-unidecode
        sudo port install py27-beautifulsoup4

5. Install [MAFFT](http://mafft.cbrc.jp/alignment/software)

    Follow the commands below to install the latest version of MAFFT:

        curl -o mafft-7.158-with-extensions-src.tgz http://mafft.cbrc.jp/alignment/software/mafft-7.158-with-extensions-src.tgz
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

6. Install [Muscle](http://www.drive5.com/muscle)

    Follow the commands below to install the latest version of Muscle:

        curl -o muscle3.8.31_i86darwin64.tar.gz http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86darwin64.tar.gz
        tar xzf muscle3.8.31_i86darwin64.tar.gz
        sudo mv muscle3.8.31_i86darwin64 /usr/local/bin/muscle
        rm -rf muscle3.8.31_i86darwin64.tar.gz

7. Install [RAxML](https://github.com/stamatak/standard-RAxML)

    Follow the commands below to install the latest version of RAxML:

        git clone https://github.com/stamatak/standard-RAxML.git
        cd standard-RAxML
        make -f Makefile.AVX.PTHREADS.mac
        rm *.o
        sudo mv raxmlHPC-PTHREADS-AVX /usr/local/bin/raxml
        cd ..
        rm -rf standard-RAxML

8. Install latest BLAST+ command line tools

    Download the OS X installation package here:
    [ncbi-blast-2.2.29+.dmg](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+.dmg).
    Double click on the downloaded file and follow instructions to install. OS X
    may not allow to install this package because of the very strict security
    settings. If that happens, go to "System Preferences -> Security & Privacy".
    Unlock the little lock icon on the bottom left of the window. Select "Allow
    apps downloaded from: Anywhere". Then try installing again.

9. At this point, we can install PhyloMill

    Type these commands to download and setup PhyloMill:

        cd ~
        git clone https://github.com/karolisr/krpy.git

    To update PhyloMill at any point you can do this:

        cd ~/krpy
        git pull

10. Run PhyloMill by typing:

        cd ~/krpy/krpy/workflows
        ./phylomill -p ~/pm_test_project -c autopilot

    This will create a sample project with name "pm_test_project" in your home
    directory and produce a sample phylogeny.

<!--
        cd ~
        git clone https://github.com/karolisr/krpy.git

        echo -e "\n# PhyloMill" | tee -a ~/.profile
        phylomill_bin_path=PATH=\"${HOME}/krpy/krpy/workflows:'${PATH}'\"
        echo -e $phylomill_bin_path | tee -a ~/.profile
        echo -e "export PATH" | tee -a ~/.profile

        phylomill_python_path=PYTHONPATH=\"${HOME}/krpy:'${PYTHONPATH}'\"
        echo -e $phylomill_python_path | tee -a ~/.profile
        echo -e "export PYTHONPATH" | tee -a ~/.profile

    Restart the Terminal app and you are done.

    To update PhyloMill at any point you can type this in the Terminal app:

        cd ~/krpy
        git pull
-->

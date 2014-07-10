#!/usr/bin/env bash
echo
################################################################################
echo "--- Will install PhyloMill"
echo
################################################################################
echo "--- First, let's update Ubuntu repositories."
sudo apt-get -y update
echo
################################################################################
echo "--- Now, let's install PhyloMill dependencies:"
echo
echo "  1. Git"
echo
sudo apt-get -y install git
echo
echo "  2. Biopython"
echo
sudo apt-get -y build-dep python-biopython
sudo apt-get -y install python-biopython
echo
echo "  3. SciPy"
echo
sudo apt-get -y install python-scipy
echo
echo "  4. Levenshtein"
echo
sudo apt-get -y install python-Levenshtein
echo
echo "  5. MAFFT"
echo
sudo apt-get -y install mafft
echo
echo "  6. Muscle"
echo
sudo apt-get -y install muscle
echo
echo "  7. RAxML"
echo
sudo apt-get -y install raxml
sudo ln -s /usr/bin/raxmlHPC /usr/bin/raxml
echo
echo "  8. NCBI-BLAST+"
echo
sudo apt-get -y build-dep ncbi-blast+
sudo apt-get -y install ncbi-blast+
echo
################################################################################
echo "--- Now, let's download PhyloMill:"
echo

install_dir=$1

if [ ! -d "$install_dir" ] ; then
   mkdir "$install_dir"
fi

cd $install_dir
rm -rf krpy

git clone https://github.com/karolisr/krpy.git
################################################################################

echo -e "\n# PhyloMill" | tee -a ~/.bashrc
phylomill_bin_path=PATH=\"$install_dir/krpy/krpy/workflows:'${PATH}'\"
echo -e $phylomill_bin_path | tee -a ~/.bashrc
echo -e "export PATH" | tee -a ~/.bashrc

phylomill_python_path=PYTHONPATH=\"$install_dir/krpy:'${PYTHONPATH}'\"
echo -e $phylomill_python_path | tee -a ~/.bashrc
echo -e "export PYTHONPATH" | tee -a ~/.bashrc

echo "Done."
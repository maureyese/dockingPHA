# Download libraries for Molecular Docking and Molecular Dynamics

1. Open PowerShell application.
2. Write ```wsl``` and press Enter.

On WSL Terminal:

**Step 1.** Update package list.

```shell
sudo apt update
```

```shell
sudo apt upgrade
```

**Step 2.** Install Python.

```shell
sudo apt-get install python3-pip
```

Verify installation

```shell
pip --version
```

**Step 3.** Install Anaconda.

```shell
wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
```

```shell
bash Miniconda3-latest-Linux-x86_64.sh
```

```shell
nano ~/.bashrc
```

It will open a text editor.

Copy and paste the following line at the end of the file.

```
export PATH="$HOME/miniconda3/bin:$PATH" 
```

To exit the text editor:

1. Ctrl + x
2. Press Y
3. Press Enter

```shell
source ~/.bashrc
```

Verify installation:

```shell
conda --version
```

(OPTIONAL) Only if not installed:

```shell
~/miniconda3/bin/conda init
```

**Step 4.** Create Conda-based Python environment

A Python environment is a place where you install Python libraries.

```shell
conda create -n vina python=3
```

```shell
conda activate vina
```

**Step 5.** Install Autodock Vina

```shell
conda install -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme
```

```shell
sudo apt install autodock-vina
```

```shell
git clone https://github.com/ccsb-scripps/AutoDock-Vina
```

```shell
cd AutoDock-Vina
```

```shell
git checkout meeko_wrap
```

```shell
cd build/python
```

```shell
python setup.py clean
```

```shell
python setup.py build
```

```shell
sudo pip python3 setup.py install
```

```shell
pip install .
```

**Step 6.** Verify Autodock Vina installation

```shell
ls ~/miniconda3/envs/vina/lib/python3.12/site-packages/vina
```

Go back to home directory:

```shell
cd
```

Ensure Python recognizes packages

```shell
python -m vina
```

**Step 7.** Install other packages

```shell
conda install numpy pandas biopython MDAnalysis
```

```shell
conda install conda-forge::gromacs
```

**Step 8.** Install AutoDockTools

```shell
git clone https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3.git
```

```shell
cd AutoDockTools_py3/
```

```shell
conda install conda-build
```

```shell
conda develop .
```

## Set Conda-based Python Environment on VS Code

Go back to home directory:

```shell
cd
```

Access to folder of GitHub Repository

```shell
cd dockingPHA/
```

Open Visual Studio Code:

```shell
code .
```

On the bottom-left corner you will see the version of Python. Click on it and select the conda environment ```vina```.
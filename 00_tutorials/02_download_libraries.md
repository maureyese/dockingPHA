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

## Conda-based Python environment

**Step 4.** Create Conda-based Python environment

A Python environment is a place where you install Python libraries.

```shell
conda create -n vina python=3.12
```

```shell
conda activate vina
```

## Installing Autodock Vina (Two methods)

### Method 1. Installing through source code

**Step 5.** Install Autodock Vina

Install these libraries Autodock Vina needs

```shell
conda install -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme
```


Install Autodock Vina Program on WSL.

```shell
sudo apt install autodock-vina
```

Retrieve the Autodock Vina Python Library Source code from GitHub.

```shell
git clone https://github.com/ccsb-scripps/AutoDock-Vina
```

Access to Autodock Vina Python library folder.

```shell
cd AutoDock-Vina
```

Start installation.

```shell
git checkout meeko_wrap
```

A message should appear:

Branch 'meeko_wrap' set up to track remote branch 'meeko_wrap' from 'origin'.
Switched to a new branch 'meeko_wrap'

```shell
cd build/python
```

```shell
python setup.py clean
```

A message should appear. At the end should say:

running clean

```shell
python setup.py build
```

It should print a long trash message with some colors.

```shell
python setup.py install
```

```shell
pip install .
```

A message should appear:

Successfully built vina
Installing collected packages: vina
  Attempting uninstall: vina
    Found existing installation: vina 1.2.7.dev38+g7b00bc5
    Uninstalling vina-1.2.7.dev38+g7b00bc5:
      Successfully uninstalled vina-1.2.7.dev38+g7b00bc5
Successfully installed vina-1.2.7.dev38+g7b00bc5

**Step 6.** Verify Autodock Vina installation

```shell
ls ~/miniconda3/envs/vina/lib/python3.12/site-packages/vina
```

You should see the Python Scripts from the Autodock Vina Python library:

__init__.py  __pycache__  _vina_wrapper.cpython-312-x86_64-linux-gnu.so  utils.py  vina.py  vina_wrapper.py

Remove the Autodock Vina folder:

```shell
cd ../../../
```

You should now be on DockingPHA folder. Now remove the "Autodock-Vina" folder.

THIS STEP IS VERY IMPORTANT BECAUSE WE DON'T WANT TO UPLOAD THE FOLDER TO GITHUB.

```shell
rm -rf AutoDock-Vina/
```

Ensure Python recognizes packages

```shell
python -m vina
```

There two possible messages:

_Message 1 (SUCESSFUL):_

```
/miniconda3/envs/vina/bin/python: No module named vina.__main__; 'vina' is a package and cannot be directly executed
```

If you have a other message which is not that one, try the second method.

### Method 2. Installing through pip

Although this method is simpler, I have found some errors with Autodock Vina through this one, so it could potentially be unstable than the first one.

Remove the current vina library from your Conda Environment.

```shell
pip uninstall vina
```

Install directly from pip

```shell
pip install vina
```

## Installing other packages

### Molecular docking

**Step 7.** Install other packages

```shell
conda install numpy pandas biopython MDAnalysis
```

### (Optional) Flexible docking

In flexible docking you select specific amino acids, instead of simulating the whole structure.

Although this project does not mandate flexible docking analysis, I still added the steps to install AutodockTools: A library used as part of flexible docking workflow.

THIS IS OPTIONAL AND YOU MAY FACE MULTIPLE ERRORS WHICH MAY BE DIFFICULT FOR YOU TO SOLVE.

I encourage you to only install it if I'm available to help you. Otherwise it will be difficult and ChatGPT (or any other AI) may not be as useful as expected.

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
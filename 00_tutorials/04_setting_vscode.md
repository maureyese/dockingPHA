# Set VSCode for Analysis

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

**Step 2**. Access to folder of GitHub Repository

```shell
cd dockingPHA/
```

Open Visual Studio Code:

```shell
code .
```

## Install Python extension on VSCode 

On VS-Code left side you will find multiple icons:

1. "Explorer" (Icon: Two stacked files)
2. "Search" (Icon: Magnifying glass)
3. "Source control" (Icon: Dots connected through lines)
4. "Run and Debug" (Icon: Triangle with a bug on the bottom corner).
5. "Extensions" (Icon: Three squares with a rhombus)

### Explorer Section

The most important is "Explorer" where you can see the files and folders from our GitHub Repository.
(Note: Please ensure you do not have the "Autodock-Vina" folder. Otherwise, REMOVE IT).

The files/folder names can have three different colors:

1. "White": The file/folder has been untouched.
2. "Green": A new file/folder has been created for the first time.
3. "Orange": A file/folder has been modified.

These colors correspond depending on the current status of the GitHub Repository.

In other words, your folder will have differences with the GitHub Repository on the website.

If you want your files to be up to date between your folder and GitHub Repository. Go to the "Source Control" Section.

### Source Control Section

In this section, you can manage git (The program which helps you upload new or modified files to your GitHub repository).

You will see three main parts:

1. A text box with title: "Message (Ctrl + Enter to commit on "branch-name")

In branch-name should say your name.

2. A blue botton saying "Commit".

3. A list with title "Changes": Files modified (with letter "M") or created (with letter "U") are listed here, which are the ones that are not up to date with GitHub Repository.

To get GitHub repository up to date:

1. Write a short message on the text box (For example: "docking script added")
2. Click on "Commit" button.
3. A new button will appear: "Sync Changes". Click on it.
4. Now go to GitHub on your browser: https://github.com/maureyese/dockingPHA.git
5. Select your branch.
6. Now you will see your files up to date.

### Extensions Control Section

We will need to install some extensions to work properly on Visual Code.

Click on the Extensions Control Section. Search:

1. "Python" and install: "Python Debugger", "Python", "Pylance" and "Python Environment". All of them are from Microsoft with a verified blue icon at the beginning of "Microsoft" name.
2. "Jupyter notebook and install: "Jupyter", "Jupyter Notebook Renderers", "Juypter Keymap". All of them are from Microsoft with a verified blue icon at the beginning of "Microsoft" name.

### Connect Conda Environment with VS Code

This is the last part of this journey.

Create a Python file (Do not save it):

1. Click on "File".
2. Click on "New file..."
3. Select "Python File"

On the bottom-left corner you will see the version of Python. Click on it and select the conda environment ```vina```.

Close the file without saving.

Now create a new Jupyter Notebook file:

1. Click on "File".
2. Click on "New file..."
3. Click on "Jupyter Notebook".

On the upper-right corner you will see a button "Select Kernel".

1. Click on "Select Kernel".
2. Two options will appear: "Python environments..." or "Existing Jupyter server..."
3. Click on "Python environments...".
4. Select the conda environment ```vina```.

You're good to go.

## End
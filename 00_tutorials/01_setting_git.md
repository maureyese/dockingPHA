# Setting up git on WSL

**Step 1.** Open WSL terminal.

Search "PowerShell" Application, write ```wsl``` and press Enter.

**Step 2.** Update package list.

Copy/paste the following commands on WSL terminal.

```shell
sudo apt update
```

```shell
sudo apt upgrade
```

**Step 3.** Install git on WSL terminal.

```shell
sudo apt install git
```

**Step 5.** Verify git is installed.

```shell
git --version
```

**Step 6.** Set your GitHub username and email on git.

Remove ```username``` from the following line and write your GitHub username instead.

```shell
git config --global user.name "username"
```

Remove ```your_email@example.com``` from the following line and write your email (registered on GitHub) instead.

```shell
git config --global user.email "your_email@example.com
```

**Step 7.** Set default branch as ```main```.

```shell
git config --global init.defaultBranch main
```

**Step 8.** Review your settings.

```shell
git config --list
```

## Connect git with Github

**Step 9.** Create a ssh key on WSL terminal.

Ensure you're on home directory. Write:

```shell
cd
```

Remove ```your_example@example.com``` from the following line and write your email (registered on GitHub) instead.

```shell
ssh-keygen -t ed25519 -C "your_email@example.com" -f .ssh/gogec
```
It will ask you: ```Enter passphrase (empty for no passphrase):```. Just press Enter twice (Do not write anything).

**Step 10.** Ensure SSH Key was created.

```shell
ls -la ~/.ssh/
```

Search for files ```gogec``` and ```gogec.pub```.

**Step 11.** Start SSH agent.

```shell
eval "$(ssh-agent -s)"
```

**Step 12.** Add your SHH PRIVATE key to the agent (PRIVATE!!!!!)

```shell
ssh-add ~/.ssh/gogec
```

**Step 13.** Open GitHub.

1. Go to your browser (Google Chrome, Microsoft edge) and access to GitHub.
2. Login.
3. Click on your image at the right-upper corner.
4. Click on "Settings".
5. Click on "SSH and GPG keys".
5. Click on "New SSH key".
6. Write on "Title": GOGEC

**Step 13.** Copy the content of SSH key

Go back to the WSL terminal.

THIS IS THE MOST IMPORTANT PART OF THE PROCESS.

```shell
cat ~/.ssh/gogec.pub
```

It will print a code. For exmaple: ```ssh-ed25519 FO2R093N...EINVIV your_mail@example.com```.

Copy it.

**Step 14.** Add your SSH PRIVATE key on GitHub.

1. Go back to GitHub on your browser.
2. Paste the SSH key on "Key".
3. Click on "Add SSH key".

DO NOT SHARE YOUR SSH PRIVATE KEY TO ANYONE.

**Step 15.** Test your conection.

Go back to WSL Terminal.

```shell
ssh -T git@github.com
```

It should print a message like:

```Hi username! You've successfully authenticated, but GitHub does not provide shell access.```

## Clone this repository

**Step 17.** Pass me your username.

Pass your username to me (maureyese), so I can make you a collaborator of this repository.

A GitHub mail will be sent to your email to confirm be part of the repository.

**Step 16.** Git clone.

The following line will copy the content of the GitHub repository to your WSL.

```shell
git clone git@github.com:maureyese/dockingPHA.git
```

**Step 17.** Git status.

The following line will check if you have any untracked file.

```shell
git status
```

**Step 18.** Create a personal branch.

This repository stores files on a branch named ```main```. It is the branch displayed on the GitHub repository URL.

You can create other secondary branches (so you do not modify the ```main``` one).

For GOGEC, each one will create its own branch.

Remove ```branch-name``` and write your name instead (example: ```mau```).

```shell
git checkout -b branch-name
```

You will now work on your own branch.

Verify your working on your own branch.

```shell
git branch
```

**Step 19.** Create a file.

To ensure you're working on your own branch, we will create a file.

```shell
nano test.txt
```

It will open a text editor. Write on it:

```
This is a message to check I'm working on my own branch.
```

To exit text editor:

1. Press Ctrl + x. 
2. Press "Y". 
3. Press Enter.

**Step 20.** Git add.

Give git the files you want to upload to the GitHub repository.
There are two options:

- _Option 1_. Add all your files.

```shell
git add .
```

- _Option 2_. Otherwise you can add files manually.

```shell
git add test.txt
```

**Step 20.** Git commit.

Each time git modifies your GitHub repository it will need a message.
This will help you track your changes on the repository.

```shell
git commit -m "test.txt created"
```

In this case, we will use the message ```test.txt created```, but don't forget to change it each time you use the command.

The message should be as short as possible.

**Step 23.** Check what you're about to commit.

```shell
git status
```

It should appear your ```test.txt``` file.

**Step 22.** Git push.

This final step will finally update your GitHub repository (on your own branch).

Remove ```branch-name``` and write your name (branch name).

```shell
git push -u origin branch-name
```
**Step 23.** Check your push was successful.

```shell
git status
``` 

```shell
git log --oneline -3
```

## If you want to practice git...

**DO NOT USE THIS REPOSITORY**.

**Step 24.** Go back to your home directory

```shell
cd
```

**Step 25**. Create a new directory (a new folder).

Remove ```folder_name``` and write the name of the folder. DO NOT REMOVE THE ```/``` AT THE END!

```shell
mkdir folder_name/
```

**Step 26**. Access to the directory.

```shell
cd folder_name/
```

**Step 27**. Create a README file.

All repositories need a README file.

```shell
echo "# Project Title\n\n*" >> README.md
```

**Step 28**. Activate git.

```shell
git init
```
**Step 29** git add, git commit and git branch.

```shell
git add README.md
```

```shell
git commit -m "Initial commit"
```

Set default branch (main).

```shell
git branch -M main
```
**Step 30.** Create GitHub repository.

1. Go back to your browser and enter to GitHub.
2. Click on the "+" button on the upper-right corner.
3. Click on "New repository"
4. On "Repository name" write the same name as your folder.
5. On "Description" write anything you want.
6. Click on "Create repository".

It will open your GitHub repository page.

Copy the last two commands from the section "…or create a new repository on the command line". Go back to WSL terminal.

The first one should look like:

```git remote add origin git@github.com:maureyese/tesst.git``` (This is an example, do not copy/paste it). Copy the one GitHub gives you and paste it on WSL terminal.

The second one is:

```shell
git push -u origin main
```
You can copy this one and paste it on WSL terminal.

## Overall workflow

Each time you will work with this repository (or any other), you should always follow these steps.

1. Open PowerShell application.
2. Write ```wsl``` and press Enter.
3. Access to the repository using ```cd folder_name/``` (Change ```folder_name``` to the actual folder name).
4. Check if your folder and GitHub repository are up to date using ```git status```
5. Write ```code .``` to access to Visual Code Studio.
6. Do some fun bioinformatics stuff, create or modify files, etc...
7. Once you're done, update your modifications using ```git add .``` or write manually the files ```git add file.txt script.py```.
8. Commit the changes with a message using ```git commit -m "message"``` (Change ```message```to the actual message).
9. Verify what you will commit using ```git status```
10. Push changes ```git push origin branch-name``` (Modify ```branch-name``` to the actual branch-name. For this repository is your name. If your working with other repository, it would be ```main```).
11. Verify everything's up to date using ```git status``` and ```git log --oneline -3```.
12. You're done.

## Finish of the document
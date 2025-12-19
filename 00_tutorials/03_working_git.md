# Collaborative Workflow with Git: Pulling Updates from Main

This guide will show you how to work collaboratively on a shared GitHub repository, keeping your personal branch up-to-date with changes from the main branch.

## Prerequisites

- You have completed the Git setup tutorial (`01_setting_git.md`)
- You have cloned the repository and created your personal branch
- You are comfortable with basic Git commands (add, commit, push)

## Overview

When working collaboratively:
- The `main` branch contains the official, stable version
- Each collaborator works on their own branch
- Periodically, you need to **pull** changes from `main` to keep your branch updated
- This prevents merge conflicts and ensures everyone is working with the latest code

## Step-by-Step Workflow

### 1.Start Your Work Session

Open Powershell and write:
```shell
wsl
```

Navigate to your repository

```
cd dockingDPA/
```

### 2. Check Your Current Branch

```shell
git branch
```

You should see:
- Your branch (with a star `*` next to it)
- The `main` branch
- Possibly other collaborators' branches

If you're not on your branch please run this:

Change branch-name to your name.

```shell
git checkout branch-name
```

Check if you're up to date

```shell
git status
```

A message should appear:

On branch main
Your branch is up to date with 'origin/main'.

(OPTIONAL PART =======================) 

If the message continues with:

Untracked files:
  (use "git add <file>..." to include in what will be committed)
        00_tutorials/03_working_git.md

nothing added to commit but untracked files present (use "git add" to track)

Please run:

1. ```git add .```
2. ```git commit -m "Add a message"``` (Don't forget to change the message name).
3. ```git push origin branch-name``` (Change to your branch).

(OPTIONAL PART =======================) 

Pull the latest status of GitHub repository

```shell
git fetch origin
```

Check if you view all branches

```shell
git branch -a
```

Make sure you're on your branch

### 3. SAFEST WAY: Retrieve files from main to your branch

Change ```{file/folder}``` to the actual filename or foldername. This command will help you retrieve a file or folder from main. 

For safety reasons, please retrieve files and folder one-by-one.

```shell
git checkout main -- {file/folder}
```

### 4. Upload newest changes to your branch

Now you can work with files and folders in your own branch. 

Once you finish remember to run the following commands:

```shell
git add .
```
```shell
git commit -m "Brief description of changes"
```
```shell
git push origin your-branch-name
```
Verify
```shell
git status
```
```shell
git log --oneline -5
```

## Daily commands workflow

1. Access to the folder:

```shell
cd dockingPHA/
```

2. Check if your working on your branch:

```shell
git branch
```

3. Check if you're up to date:

```shell
git status
```

(git add. ; git commit -m "message", git push origin branch-name; if needed)

4. Get files from main branch if needed

```shell
git checkout main -- {file/folder}
```

5. Work in VS Code:

```shell
code .
```

6. Stage your changes

```shell
git add .
```

7. Commit with a message

```shell
git commit -m "message"
```

8. Push to GitHub Repository

```shell
git push origin branch-name
```

## Other ways to retrieve files from main branch

### 3. OTHER WAY WITH POTENTIAL ISSUES: Switch to Main Branch

```shell
git checkout main
```

**Important:** Always make sure you have committed or stashed your changes before switching branches.

### 4. Pull Latest Changes from GitHub

```shell
git pull origin main
```

This command:
- Fetches the latest changes from the remote GitHub repository
- Merges them into your local `main` branch

### 5. Switch Back to Your Branch

```shell
git checkout your-branch-name
```

Replace `your-branch-name` with your actual branch name.

### 6. Merge Main into Your Branch

```shell
git merge main
```

This brings the latest changes from `main` into your personal branch.

### 7. Handle Merge Conflicts (If Any)

If you see a message like:
```
Auto-merging file.py
CONFLICT (content): Merge conflict in file.py
Automatic merge failed; fix conflicts and then commit the result.
```

Don't panic! This happens when both you and someone else modified the same part of a file.

**To resolve conflicts:**

1. Open the conflicting file in VS Code:
   ```shell
   code file.py
   ```

2. Look for conflict markers:
   ```
   <<<<<<< HEAD
   Your changes
   =======
   Changes from main
   >>>>>>> main
   ```

3. Decide which changes to keep (or combine them)
4. Remove the conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`)
5. Save the file

6. Mark the conflict as resolved:
   ```shell
   git add file.py
   ```

7. Complete the merge:
   ```shell
   git commit -m "Resolved merge conflict"
   ```

### 8. Push Your Updated Branch

```shell
git push origin branch-name
```

### 9. Alternative: Using Git Rebase

Some teams prefer `rebase` instead of `merge`:

```shell
# While on your branch
git rebase main
```

**Rebase vs Merge:**
- **Merge**: Creates a merge commit, preserves branch history
- **Rebase**: Rewrites history, creates linear timeline
- For beginners, `merge` is recommended as it's simpler to understand

## Best Practices for Collaboration

### 1. **Commit Frequently**
- Small, focused commits are better than one huge commit
- Each commit should represent one logical change

### 2. **Write Good Commit Messages**
```
Add: New script for data preprocessing
Fix: Corrected error in energy calculation
Update: Improved documentation for docking workflow
```

### 3. **Pull Before You Start**
Always pull from `main` at the beginning of each work session.

### 4. **Communicate with Your Team**
- Let others know what files you're working on
- Discuss major changes before implementing them
- Use GitHub Issues or pull requests for discussion

### 5. **Keep Your Branch Focused**
- One branch = one feature/task
- Don't mix unrelated changes in the same branch

## Using VS Code for Collaboration

### 1. **Source Control Panel**
- Shows changes between your branch and `main`
- Visual conflict resolution tools
- Easy staging and committing

### 2. **Branch Management**
- Click on branch name in bottom-left corner
- Switch branches with one click
- Create new branches easily

### 3. **Pull Changes Directly**
In VS Code:
1. Click Source Control panel ( )
2. Click "..." (More Actions)
3. Select "Pull" or "Pull from..."

## Troubleshooting

### "Your local changes would be overwritten"
```shell
# Stash your changes
git stash

# Pull from main
git pull origin main

# Apply your changes back
git stash pop
```

### "Already up to date"
Your branch is already synchronized with `main`. No action needed.

### "Permission denied"
Make sure:
- You're using SSH (not HTTPS)
- Your SSH key is properly added to GitHub
- You're a collaborator on the repository

## Pull Request Workflow (Advanced)

When your work is complete:
1. Push your final changes
2. Go to GitHub.com
3. Create a "Pull Request" from your branch to `main`
4. Request review from teammates
5. Address feedback if needed
6. Merge into `main` after approval

## Summary

Working collaboratively requires:
1. **Regular synchronization** with `main`
2. **Clear communication** with your team
3. **Organized workflow** (pull → work → commit → push)
4. **Conflict resolution skills**

By following this guide, you'll minimize conflicts and ensure smooth collaboration on the shared repository.

---
**Remember:** Always pull before you start working, and push when you're done!
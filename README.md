# TAESOPT.jl

Transport Aircraft and the Environment System OPTimization - 

Julia based TASOPT++ includes PCEC, hydrogen fuel options etc. Perhaps need some cooler name but can always change later - how about ZITA: Zero Impact Transport Aircraft

## NPSS integration

Currently NPSS is being used to model the turboshaft engine by simply "file wrapping" - Julia writes to a file, NPSS reads that, computes and outputs to a file that julia can then read.

> Update (Feb 12th 2021):

> Using NPSS to do the propulsion system calculations becomes slow if NPSS has to be started/ compiled and then the engine model solved for each call to the propulsion system. 
> To circumvent this issue Julia now starts up an asynchronous process on your system that runs NPSS and then writes engine inputs to the `STDIN` of the NPSS process and waits till NPSS returns a success/ failure code on its `STDOUT` before continuing onto the next calcualtions. 
> This ideally means that we can get rid of file-based transfer but, that remains to be implemented


## Collaboration guide

Some guidelines to add to this repo

### Work in branches

Usually you should be working on your own fork, **exceptions** are if you have been added as a collaborator to this repo. If you think you are going to be writing large chunks and would like to be added to this repo as a collaborater contact Prashanth, so you can directly modify this repo, without having to create pull requests (you should still be working on separate branches and then you can just do pull requests).
Don't commit anything to the master branch. Here's how you create your own fork and branch.  
First create a fork of this repo by clicking on "fork" on the top right hand of the github page. This creates a copy of this repo that is separate form this one. This ensures that any chagnes made to your fork will not affect other's forks.

You will now need to clone your forked version of the repo to your machine where you will be writing code. This [page](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) has a good overview of how to do this whole process.
After this I highly recommend using branches within your own fork and I prefer the git command line interface over the web interface.
Use this command to create a new branch for example `dev_prash`

    git checkout -b <name_of_branch>

as you add files, use the following command to tell git to track your file

    git add <path to file, or use . to track all files>

once you are happy with the edits you have made, you have to *commit* the edits by doing

    git commit -a -m "COMMIT MESSAGE"

make your commit messages useful. These messages along with frequent commits will create a rich history to the repo so as this grows you will be able to track changes and debug more efficiently

After committing to your local machine, you need to *push* these changes to the remote github server. The _first_ time you do this for each branch you need to tell the remote repo what to track
```bash
    git push -u origin <name_of_branch>
```
for subsequent pushes, you only need to do `git push`.

### Sync your repo with this main repo periodically

Keep syncing your fork to this repo regularly so that you have all the latest bug fixes and changes that others have done. 

To sync your fork, `cd` into the right folder on your terminal and do the following
```bash
    git fetch upstream
```
This fetches all the changes made to the master repo. Then switch to your master branch (you should be working in branches even in your own repo)
```bash
   git checkout master
```
Once you are in your master branch do:
```bash
   git merge upstream/master
```
This brings your fork's `master` branch into sync with the upstream repository, without losing your local changes.

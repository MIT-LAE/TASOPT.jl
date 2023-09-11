# TASOPT.jl

Transport Aircraft and the Environment System OPTimization - 

Julia based TASOPT++ includes PCEC, hydrogen fuel options etc. Perhaps need some cooler name but can always change later - how about ZITA: Zero Impact Transport Aircraft

## Getting started

### Simple install

The easiest way to run `TASOPT.jl` would be to add the package using the julia package manager using the github repository.

You can do this by starting a Julia session and then activating the package manager by typing `]` and then entering:
```julia-repl
pkg> add "git@github.mit.edu:LAE/TAESOPT.jl.git"
```

You can then import `TASOPT` as you would with any Julia package:
```julia-repl
julia> using TASOPT
```
### Local development

If you are going to develop the source code of `TASOPT.jl` you might benefit from a local clone of the git repository which
can then fit into a workflow using [`Revise.jl`](https://timholy.github.io/Revise.jl/stable/) for example.

Step 1: Clone the git repo locally
```bash
git clone git@github.mit.edu:LAE/TAESOPT.jl.git
```

Step 2: `cd` to the folder where TASOPT is cloned

Step 3: Use `Pkg` to install/ develop the package

```julia
pkg> dev .
```

You should now be able to import TASOPT from within any Julia script in your base environment.

If you are using `Revise.jl` be sure to first import `Revise` before importing `TASOPT`

```julia
using Revise
using TASOPT
```


Within the `run` folder is an example set of files to get started. To run the code on hex first `cd` into the run directory. Then you have 2 choices:

1. `julia example_run.jl` : runs the example file in the current session
2. `sbatch tas.sh example_run.jl` : submits the job to slurm

The slurm script `tas.sh` has further details on options (including emailing you when your job is done etc.) that you can set while submitting to slurm.

## NPSS integration

Currently NPSS is being used to model the turboshaft engine by simply "file wrapping" - Julia writes to a file, NPSS reads that, computes and outputs to a file that julia can then read.

> Update (Feb 12th 2021):

> Using NPSS to do the propulsion system calculations becomes slow if NPSS has to be started/ compiled and then the engine model solved for each call to the propulsion system. 
> To circumvent this issue Julia now starts up an asynchronous process on your system that runs NPSS and then writes engine inputs to the `STDIN` of the NPSS process and waits till NPSS returns a success/ failure code on its `STDOUT` before continuing onto the next calcualtions. 


## Collaboration guide

**Important**
Before submitting a new pull request (PR), go to the `test` folder and run the regression test `run_regression_test.jl` and the unit test `run_unit_test.jl`
If there is an error, it is your responsibility to edit your code and make it work.
The PR will not be reviewed if the regression or the unit test fails.

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

## Convert Fortran to Julia

The file `fort2julia.vim` was created to help this conversion. 
At first, make sure the target file (e.g., `tar_file.f`) and the vim file (`fort2vim.vim`) are in the same directory.
Then, open the target file in vim.
```bash
    vim tar_file.f
```
Next, type
```
    :source fort2vim.vim
```
Now, the file has been edited by vim. 
But we do not want it to overwrite the original file.
Thus, we save it to a different file (e.g., `tar_file.jl`), and we also need to ensure that the original file is not changed.
```
    :w tar_file.jl
```
The updated file is saved in `tar_file.jl`.
Now, we exit vim without changing the original f77 file by
```
    :q!
```

To check it, you should have two files now: `tar_file.jl` and `tar_file.f` where the former is the edited Julia file and the latter is the original f77 file.

Running `fort2julia.vim` **cannot solve all the problems**, though.
You need to debug the code to make sure: 
* The syntax is correct. Several items to check:
  * The array indexing is still using `()` inherting from f77. You need to change it to Julia's `[]`.
  * Some data type such as `parameter`, `data` are from f77 but not used in Julia. They need to be changed.
  * Some function is with a type, e.g., `real function dot(x)`, this needs to be modified.
  * Some function will modify input parameters, (say it is named as f(x)). 
    A regular Julia function cannot do it. 
    You need to let Julia know it by adding a `!` mark behind the function. 
    So the new function will be: `f!(x)` in Julia.
* It gives the same output with the original f77 code. 

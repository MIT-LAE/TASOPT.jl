# TASOPT.jl
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://mit-lae.github.io/TASOPT.jl/dev/) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mit-lae.github.io/TASOPT.jl/stable/) [![CI](https://github.com/MIT-LAE/TASOPT.jl/actions/workflows/CI.yml/badge.svg?branch=main&job=test&version=lts)](https://github.com/MIT-LAE/TASOPT.jl/actions/workflows/CI.yml) [![codecov](https://codecov.io/github/MIT-LAE/TASOPT.jl/graph/badge.svg?token=J1FNXGO3SD)](https://codecov.io/github/MIT-LAE/TASOPT.jl) [![DOI](https://joss.theoj.org/papers/10.21105/joss.08521/status.svg)](https://doi.org/10.21105/joss.08521)


Transport Aircraft and the Environment System OPTimization (TASOPT) implemented in Julia. Originally based on Mark Drela's [FORTRAN code](https://web.mit.edu/drela/Public/web/tasopt/) of the same name.

## Getting started

### Simple install

The easiest way to run `TASOPT.jl` would be to add the package using the julia package manager using the github repository.

You can do this by starting a Julia session and then activating the package manager by typing `]` and then entering:
```julia-repl
pkg> add TASOPT
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
git clone git@github.com:MIT-LAE/TASOPT.jl.git
```

Step 2: `cd` to the folder where TASOPT is cloned

Step 3: Use `Pkg` to install/ develop the package

```julia
pkg> dev .
```

You should now be able to import TASOPT from within any Julia script in your base environment.

Note: If you clone another version of TASOPT, `using TASOPT` will always use the directory where `dev .` was used.

If you are using `Revise.jl` be sure to first import `Revise` before importing `TASOPT`

```julia
using Revise
using TASOPT
```

## Collaboration guide

**Important**
Before submitting a new pull request (PR), go to the `test` folder and run the tests by doing 
```julia
using TASOPT, Pkg
Pkg.test("TASOPT")
```
If there is an error, it is your responsibility to edit your code and make it work.
The PR will not be reviewed if the regression or unit test fails. If you find that the tests do not capture the right behavior or are flawed, please raise an issue.

### Work in branches

Don't commit anything to the main branch. Here's how you create your own fork and branch.  
First create a fork of this repo by clicking on "fork" on the top right hand of the github page. This creates a copy of this repo that is separate form this one. This ensures that any changes made to your fork will not affect other's forks.

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
This fetches all the changes made to the main repo. Then switch to your main branch (you should be working in branches even in your own repo)
```bash
   git checkout main
```
Once you are in your main branch do:
```bash
   git merge upstream/main
```
This brings your fork's `main` branch into sync with the upstream repository, without losing your local changes.

# Observatory
[![Build status](https://travis-ci.org/Freshwater-Initiative/Observatory.svg?branch=master)](https://travis-ci.org/Freshwater-Initiative/Observatory)

Tools for observing the terrestrial and aquatic surfaces of Earth 

"You belong among the wildflowers, you belong somewhere you feel free." Tom Petty

An observatory is a location used for observing terrestrial or celestial events (thank you Wikipedia). Observatories have been as simple as containing an astronomical sextant, and as complicated as modern academic supported observatories containing multi-million dollar instruments, tools, with institutions supporting long term research and education programs.  While observatories are usually thought of as star-gazing investments in the field of astronomy, observatories have also been constructed in climatology/meteorology, geophysical, oceanography and volcanology communities, in order to investigate and coordinate their research efforts.  

This repository is intended for the sharing and distribution of open-source Python based code useful for model and data integration that improves access to large datasets, reduces computational burden, reinvent the wheel less often, and share and communicate more about how to synthesize earth surface observations in useful ways.

## Installing Version on conda-forge
Install package with conda:

```bash
conda install -c conda-forge ogh
```

Execute from Jupyter Notebook

```
!conda install -c conda-forge --yes ogh
import ogh
```

## Installing Latest Master Version

Linux/OSX:

```bash
wget https://raw.githubusercontent.com/Freshwater-Initiative/Observatory/master/requirements.txt
wget https://raw.githubusercontent.com/Freshwater-Initiative/Observatory/master/requirements-dev.txt
conda create -n oghenv -c conda-forge python=2.7 --file requirements.txt --file requirements-dev.txt
source activate oghenv
pip install git+https://github.com/Freshwater-Initiative/Observatory.git
```

# Work with a git-versioned-folder in hydroshare to develop your own Utilities
0) Make a fork of Freshwater-Initiative/Observatory
1) In HydroShare, get to JupyterHub and open up a terminal instance.
2) Change the working directory to notebooks/utilities
## If you haven't cloned this repository into Hydroshare yet:
3) type/copy in "git clone <github link>" that is available from your fork (eg., https://github.com/username/Observatory.git)
4) then type in your github username and password to then download the git clone.
  $ git config --global user.name "your git username"
  $ git config --global user.email "your email that you used to setup the git account"
5) you should now have notebooks/utilities/Observatory subdirectory with this README.md and the observatory_gridded_hydrometeorology.py (OGH) within.
## If you have previously cloned the git folder to notebooks/utilities, update this to the latest file
3) Change the working directory to notebooks/utilities/Observatory
4) Pull the latest (your updated fork - pull from Freshwater-Initiative/Observatory master, before you do this) files from the repository
  $ git pull
5) you should now have notebooks/utilities/Observatory subdirectory that matches your fork on github.com (which should match Freshwater-Initiative/Observatory - if you pulled from the master)
## Either way:
6) Change the working directory to notebooks/utilities/Observatory
7) Copy observatory_gridded_hydromet.py to notebooks/utilities/
  $ cp observatory_gridded_hydromet.py ../
  Now, the file '/notebooks/utilities/observatory_gridded_hydromet.py' is updated to the latest state (from git) 

## Saving changes back to git repository

8) Work on the file and save the changes.
9) Change the working directory to '/notebooks/utilities'
10) Copy the modified file back to git versioned folder
  $ cp 'observatory_gridded_hydromet.py' Observatory/
11) Change the working directory to '/notebooks/utilities/Observatory'
12) Check file changes 
  $ git status
13) Commit the changes
  $ git add observatory_gridded_hydromet.py
  $ git commit -m 'Add message to describe these changes, if any'
  $ git push
14) Check that the files changed
  $ git status

 “Use only that which works, and take it from any place you can find it.” Bruce Lee

# notes

just trying to scrach down some notes before i forget everything about python again.

## virtualenv

I named it `scip`, could be anything on your end. BUT you need to specify python 3 at createon eg.: `virtualenv scip -p python3`

Usage:

* start with `source ./scip/bin/activate`
* end with `deactive` when done
* `pip` should point to the virtualenv binary (`./scip/bin/pip`) and thus install everything locally.
* after installing a package pls. run `pip freeze -l > requirements.txt`. Probably better to add the actual package manually since pip and virtualenv is full of old unfixed bugs... ( https://bugs.launchpad.net/ubuntu/+source/python-pip/+bug/1635463 )

## file structure and git

### git

pls ignore every large files and user specific configs (eg. virtual env and .vscode/ folders) also ignore large input data files (i dont know how to setup git LFS yet, so we share those with some other method probably gdrive or 1drive) .probably everthing will be good if we stick to the same folder structure described below.

### folder structure

defalt folders and what files belong to them (:

* root: only global configs and notes
* src: all source code should be under this folder, sub folder structre will be determined on the fly.
* data: large data files to be ignored.
* sample: small sample data files which can be shared to test correctness.

## List of stuff to ASK from Zoli

* numerical precision ? currently we have a 5.748208131706178e-09 difference when summing up cummulative probabilities. which is 1 order smaller than the average probability. Tried to manually sum up probabilities and we still have about the same error. (the difference beetwen numpy.sum and my safe sum is in the serveral order of magnitudes smaller ~1e-13)

* how do we attach amino acids together (?). I tried to reexamine the task and concluded that we need the psi angle from the 'left' amino acid and the phi angle of the 'right' amino acid.

TODO:

* database building
* database query
* next aa structure (gets an aa type and previous/next aa type, and returns a random structure with the help of the distribution draw and database query)
* structure builder (adds a new aa, with given 3d structure into an existing structure)
* sequence iterator (iterates over the sequence call the random draw and database query and ther structure builder)
* collision detector (to be called at the end)

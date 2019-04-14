# ensemble-superset-generator

Project for the scientific python course
[![Build Status](https://github.drone.morasz.hu/api/badges/morbalint/ensemble-superset-generator/status.svg)](https://github.drone.morasz.hu/morbalint/ensemble-superset-generator)

## Setup steps

1. Clone the repo and open it. ( `git clone https://github.com/morbalint/ensemble-superset-generator`)
1. Install python 3, pip 3 and virtualenv (`sudo apt install python3 virtualenv`)
1. create a virtualenv with the python 3 interpreter, give it whatever name you want, we will use `scip` in this document (`virtualenv -p python3 scip`)
1. add the folder created by virtualenv to git ignore, if its not already ignored. (virtualenv creates a single folder named exactly as the environment name given in the above command, if you used scip it is already ignored)
1. Activate the virtual environment (`source ./scip/bin/activate`)
1. install packages with pip from requirements.txt (`pip install -r requirements.txt`)
1. run tests, write code, interpret it, (this is still a work in progress)
1. git commit & push

(TODO: move to docker)
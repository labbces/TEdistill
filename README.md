# panTE

This is our script to generate a panTE library based on genome-specific EarlGrey TE libraries. It is heavily inspired in the panEDTA.sh from EDTA

# Installing

Create a python environment in which to install some dependencies, as biopython

```
python -m venv .venv
. ./.venv/bin/activate 
pip install biopython
```

Then, clone out repo:

```
git clone git@github.com:labbces/panTE.git
```

Then cd into the repo and add submodules

```
git submodule init
git submodule update
```

This will add out own copy of the EDTA repository, as we are using some of their scripts to generate the panTE

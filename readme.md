# Genetic diagnosis of Mendelian disorders

Template for the Genetic diagnosis of Mendelian disorders project.

## Load python environment

* make sure our `wBuild` repo is cloned at the same level as this project 
* load latest anaconda `module load i12g/anaconda`
* check if snakemake command is available `which snakemake`

if not

* install it via `pip install snakemake`


## Run Snakemake
Available parameters:

- all (default): `snakemake all`
- show (opens output in chrome)
- publish (rsyncs output to webdir)
- graph (makes nice dependency graph; requires graphviz)

## FAQ

* What will be rendered?

A: All R scripts with a `wbuild` header in `./Scripts`. 
The folder `./src/` will not be searched for report files.

* What is a valid `wbuild` header?

```
#' wb:
#'   input: "Data/table0.tsv"
#'   output: [
#'     "Data/table1.tsv", 
#'     "Data/table2.tsv"
#'   ]
```

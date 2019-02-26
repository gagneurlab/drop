# Genetic diagnosis of Mendelian disorders

Goals:

* Setup pipeline from Kremer, Bader paper.
* Plug and Play for new RNA samples
    * using aberrant expression pkg
    * using splicing pkg
    * develop MAE functions/pkg
* extend pipeline by proteomics

## Diagnosis tools:
- [Webserver](https://i12g-gagneurweb.in.tum.de/project/genetic_diagnosis/)
- [TABLE of disease associated genes](https://i12g-gagneurweb.in.tum.de/project/genetic_diagnosis/#Scripts_diagnosis_tools_disease_associated_genes.html)
- [Patient reports](https://i12g-gagneurweb.in.tum.de/shinyserver/)

# Getting started with wBuild

If you have already `wBuild` installed please make sure you have the newest version

```
wbuild upgrade
```

If you do not have `wBuild` installed make sure everything is commited and then start the following commands

```
wbuild init
git checkout readme.md wbuild.yaml
```

# Literature
* Kremer, Bader genetic diagnosis [paper](https://www.nature.com/articles/ncomms15824)
* [OUTRIDER](https://www.cell.com/ajhg/pdf/S0002-9297(18)30401-4.pdf)
* Excellent review on [variant prioritization](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5935497/)

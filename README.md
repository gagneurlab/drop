# Detection of RNA Outlier Pipeline

installation

```
# clone repos
git clone git@i12g-gagneurweb.informatik.tu-muenchen.de:mumichae/drop_analysis.git
git clone git@i12g-gagneurweb.informatik.tu-muenchen.de:salazar/drop.git

# init submodules
cd drop
git submodule update --init --recursive

# install drop package
pip install -e .

# init wbuild (still uses repo folder not user folder!)
(cd drop/modules/aberrant-splicing-pipeline;
    mkdir .wBuild;
    wbuild update;)

# how to run one analysis? (this should be through drop_analysis itself!)
cd drop/modules/aberrant-splicing-pipeline
snakemake -n --configfile ../../../../drop_analysis/config.yaml 

# something like
# cd ../drop_analysis
# snakemake -n splicing
```


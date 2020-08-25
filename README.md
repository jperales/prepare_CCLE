# Prepare CCLE
> Author: Javier Perales-Patón - javier.perales@bioquant.uni-heidelberg.de - ORCID: 0000-0003-0780-6683

> NOTE: This repository is supposed to be used as `git submodule` in another repository using branches as tagged versions.

CCLE multi-omics data is ensembled into a [`MultiAssayExperiment`](http://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) (R BioC).

Current omics added are: see [raw/README.md](raw/README.md) for details.

Visit [INSTALL.md] to setup the environment. Then run
```
conda activate envs/mae
make
# 1. Download multi-omics data files
# 2. Build the MultiAssayExperiment (see data.md)
```

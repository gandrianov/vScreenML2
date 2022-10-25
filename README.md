VScreenML_v2
==============================

Modification of original [VScreenML](https://www.pnas.org/content/117/31/18477) module


**Required packages**

- [PyRosetta](https://www.pyrosetta.org/downloads)
- [RDKit](https://anaconda.org/conda-forge/rdkit)
- [XGBoost](https://anaconda.org/conda-forge/xgboost)
- [ODDT](https://anaconda.org/oddt/oddt)

The most packages should be installed automatically with VScreenML2. But PyRosetta should be installed manually because of licensing. To obtain credentials for PyRosetta, you need to apply for licence on [UW website](https://els2.comotion.uw.edu/product/pyrosetta). In the end you will get login and password that will be used for the installation:

```
conda install -c https://<LOGIN>:<PASSWORD>conda.graylab.jhu.edu pyrosetta
```

pip install https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python310.mac/PyRosetta4.Release.python310.mac.release-333.tar.bz2


**Installation**

```
conda install -c gandrianov vscreenml_v2 
```



**Usage**

The package contains several built-in scripts for different operations

**Ligand parameterization**
```
mol2genparams.py -s <LIGAND>.mol2 --prefix <LIGAND> --comment_bonds=true
```

**Protein-ligand complex minimization**
```
pyrosetta_minimizer.py -ligand <LIGAND>.pdb -params <LIGAND>.params -protein <PROTEIN>.pdb -output <COMPLEX>.pdb
```

**Features calculation**
```
calculate_features.py -complex <COMPLEX>.pdb -params <LIGAND>.params -output <FEATURES>.csv
```

**VScreenML score prediction**
```
predict_class.py -features <FEATURES>.csv -output <PREDICTIONS>.csv
```

### Copyright

Copyright (c) 2022, Grigorii Andrianov

## Dalton Lab - UNICAMP - *trRosetta_distance-constraints*
Script to filter distance constraints from the *NpzFile* generated by trRosetta and create a Rosetta constraint file using the Spline function. This constraint file can be used to model protein structure from sequence data with [Rosetta](http://new.rosettacommons.org/docs/latest/Home).

This script uses the ***trRosetta*** inter-residue contacts/distances/orientations predictions from the protocol developed in: [Improved protein structure prediction using predicted inter-residue orientations](https://www.pnas.org/content/117/3/1496) and available in the [trRosetta](https://yanglab.nankai.edu.cn/trRosetta/) server.

### Usage
```
python3 ./npz_converter.py --npz file.npz --fasta file.fasta --threshold 0.9 --target 'target_name' --path 'path'
```
By default, this script will save a Rosetta constraint file (Spline function) and spline files for all residue pairs inside the threshold (0.9 by default, but customizable) in the specified path. The spline files must remain in the path passed as argument because the path is specified in the constraint file for each spline file. Optional arguments to save the spline function plots and not saving the constraint/spline files are available:
```
--constraints False --graph True
```

### Example Files
The fasta file and *NpzFile* (distribution matrix) of the CASP13 T0957s1 target is being provided as an example.

#### Note
- Rosetta is freely available to academic laboratories. 
- Obtaining a license is free for academic users. 
- Information on how to download and install Rosetta can be found [here](https://www.rosettacommons.org/demos/latest/tutorials/install_build/install_build). 
- Unfortunately, there is no support for the whole Rosetta on Windows. Dual booting or virtual machines running Linux/MacOS are options.

###### For bugs and suggestions, you can reach me by <a href="mailto:amiquilini@id.uff.br?subject=Github%20Contact">E-mail</a>

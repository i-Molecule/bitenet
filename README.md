# BiteNet

The default copyright laws apply. To request commercial or academical license please contact : p.popov@skoltech.ru.

If you use or refer to BiteNet please cite: 
Igor Kozlovskii and Petr Popov, "Spatiotemporal identification of druggable binding sites using deep learning", Communications Biology, accepted.

If you use training set or our source code, please also cite: 10.5281/zenodo.4043664

# Abstract

Identification of novel protein binding sites expands «druggable genome» and opens new opportunities for drug discovery.

Generally, presence or absence of a binding site depends on the three-dimensional conformation of a protein, making binding site identification resemble to object detection problem in computer vision. 

Here we introduce a computational approach for the large-scale detection of protein binding sites, that considers protein conformations as 3D-images, binding sites as objects on these images to detect, and conformational ensembles of proteins as 3D-videos to analyze.

BiteNet is suitable for spatiotemporal detection of hard-to-spot allosteric binding sites, as we showed for conformation-specific binding site of the epidermal growth factor receptor, oligomer-specific binding site of the ion channel, and binding site in G protein-coupled receptor.

BiteNet outperforms state-of-the-art methods both in terms of accuracy and speed, taking about 1.5 minutes to analyze 1000 conformations of a protein with ~2000 atoms. 

## Installing

BiteNet v1.0 works with tensorflow==1.14, so we advise to make virtual environment first:
```bash
python3 -m virtualenv bitenet-env
```
and then install tensorflow (it can be non single-line for gpu)
```bash
python3 -m pip install tensorflow-gpu==1.14
```
or
```bash
python3 -m pip install tensorflow==1.14
```

You can install BiteNet directly from git
```
python3 -m pip install https://github.com/i-Molecule/bitenet.git
```

## Usage 
It is possible to run BiteNet inside python or using scripts, visualization using Pymol is also provided.

### From python

For running binding site predictions:
```python
from bitenet import BiteNet
bitenet = BiteNet()

# single .pdb file
predictions, residues = bitenet("./examples/temp.pdb")

# dataset .ds file
names, predictions, residues = bitenet("./examples/temp.ds")

# folder 
names, predictions, residues = bitenet("./examples/temp/")

# alter minibatch size if memory allocation error occurs
bitenet.minibatch_size = 4

# or other parameters
bitenet.dataloader.stride = 32  # smaller stride for grid splitting
bitenet.dataloader.rotation_eval = True # to run rotations during predict

bitenet.prediction_processer.score_threshold = 0.01 # smaller score threshold for more predictions
bitenet.prediction_processer.distance_threshold = 4 # more predictions as less predictions are filtered in non max suppression
bitenet.prediction_processer.distance_residues = 4 # distance threshold for protein residues to be considered to be on predictions interface
```

For trajectory predictions clustering:
```python
from bitenet import read_predictions
from bitenet.clustering import Clustering_MeanShift, Clustering_DBSCAN, \
    Clustering_Agglomerative, Clustering_Agglomerative_Residues

_, predictions, residues = read_predictions("predictions.log", 
    get_residues=True)

# mean shift
clustering = Clustering_MeanShift(distance_merge=5)
clustering.cluster(predictions)

# DBSCAN
clustering = Clustering_DBSCAN(eps=0.5, min_samples=5)
clustering.cluster(predictions)

# Agglomerative clustering on predictions coordinates
clustering = Clustering_Agglomerative()
clustering.cluster(predictions)
clustering.refit(n_clusters=10) # to get different number of clusters

# Agglomerative clustering on predictions residues
clustering = Clustering_Agglomerative_Residues()
clustering.cluster(predictions, residues)
clustering.refit(n_clusters=10)


print(clustering.get_summary_str(all=False))    # print clusters info: scores, coordinates, top score frames
with open("clusters.txt", "w") as file:
    file.write(clustering.get_summary_str())    # write to file
clustering.export_summary("clusters.csv", all=False)    # csv with clusters
clustering.export_summary("clusters_all.csv", all=True) # csv with not filtered clusters
clustering.plot("clusters.png") # plot cluster scores across trajectory; however better use your custom plotting for more accurate images
```

### From pymol
```python
from bitenet import BiteNet_Draw
cmd.load("temp.pdb")
model = BiteNet_Draw()
model("temp.pdb") # predict and draw predictions for file
model("temp")     # or for pymol protein object (it will just write the same pdb file)

from bitenet import read_predictions
from bitenet.clustering import Clustering_MeanShift
from bitenet.pymol_draw import draw_clusters_predictions, draw_clusters_density
_, predictions, residues = read_predictions("predictions.log", 
    get_residues=True)  # read predictions
clustering = Clustering_MeanShift()
clustering.cluster(predictions) # cluster
draw_clusters_predictions(clustering)   # draw colored predictions
draw_clusters_density(clustering)       # draw colored densities for clusters
```

### Scripts

*Predict* script runs prediction for single or multiplt pdbs:
```bash
bitenet [path] [out]
```
where *path* can be single *.pdb* file, folder with *.pdb* files or *.ds* file with list of pdbs on separate lines; *out* is path to output file with predictions or folder with separate predictions if *--sep* is provided. More options are available with:
```bash
bitenet -h
```

*Test* script runs prediction and evalutate model accuracy at the same time:
```bash
bitenet-test [path]
```

*Train* script runs training for a new model. Better look at the script itself.

*Clustering* script for getting prediction clusters for predictions run on a trajectory.
```bash
bitenet-cluster [path] [out]
```
where *path* is path to prediction log file output from bitenet predict; *out* is path to new folder to output to. Output files contain clusters scores and coordinates csv files and plot of per frame cluster scores. More options are available:
```bash
bitenet-cluster -h
```

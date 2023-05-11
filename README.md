# Pangenome-construction-using-VG
Steps to the creation of a pangenome in VG using Singularity -- CSIRO INTERNSHIP 2023

#Installing and Using Singularity


This step requires conda.
installing singularity :
```sh
conda install -c conda-forge singularity 
```
creating an environment in conda :
```sh
conda create pang
conda activate pang
```
building an image in singularity :
```sh
singularity help build
// this is to see the different options to create the image
singularity build image.sif docker://quay.io/vgteam/vg:v1.48.0   
//this is the image used to launch vg
singularity shell image.sif
```
now vg can be lanuched :
```sh
vg help 
```

# Pangenome-construction-using-VG
Steps to the creation of a pangenome in VG using Singularity -- CSIRO INTERNSHIP 2023
The goal is to understand how singularity and vg work, in order to build a Pangenome.


# Installing and Using Singularity


This step requires conda.
To install singularity, use the following :
```sh
conda install -c conda-forge singularity 
```
Install SAMTools via conda :

```sh
conda install -c bioconda samtools
```
SAMTools will be used later to edit the VCF files to make subgroups with 0.1% minor allele frequency or above.

Creating an environment, here ```pang```, in conda :
```sh
conda create pang
conda activate pang
```
Building an image in singularity from an image in the docker library :
```sh
singularity help build
// this is to see the different options to create the image
singularity build image.sif docker://quay.io/vgteam/vg:v1.48.0   
singularity shell image.sif
```
Here, the image ```quay.io/vgteam/vg:v1.48.0``` is the one used to be able to launch vg
Now vg can be launched :
```sh
vg help 
```

However, it is necessary to bind the singularity container with the VCF files. 
By default Singularity bind mounts ```/home/$USER```
Assume the files are in the directory ```data``` :
```sh
singularity shell --bind data:/mnt image.sif
```

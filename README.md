# Pangenome-construction-using-VG
Steps to the creation of a pangenome in VG using Singularity -- CSIRO INTERNSHIP 2023

requires conda
installing singularity :
conda install -c conda-forge singularity  
creating an environment in conda :
conda create pang
conda activate pang
building an image in singularity :
singularity build image.sif docker://quay.io/vgteam/vg:v1.48.0   
//this is the image used to launch vg
singularity shell image.sif
now vg can be lanuched :
vg help 


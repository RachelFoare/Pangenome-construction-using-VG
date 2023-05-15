# Pangenome-construction-using-VG
Steps to the creation of a pangenome in VG using Singularity -- CSIRO INTERNSHIP 2023

The goal is to understand how singularity and vg work, in order to build a Pangenome.


# Installing and Using Singularity


This step requires conda.
To install singularity, use the following :
```sh
conda install -c conda-forge singularity 
```
Install SAMTools (Sequence Alignment/Map) via conda :

```sh
conda install -c bioconda samtools
```
SAMTools will be used later to edit the VCF files to make subgroups with 0.1% minor allele frequency or above, and to merge all VCF files into one.

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
Here, the image ```quay.io/vgteam/vg:v1.48.0``` is the one used to be able to launch vg.

Let's check if it worked properly :
```sh
singularity inspect image.sif
```
It should return something like the following :
```sh
org.label-schema.build-arch: amd64
org.label-schema.build-date: Thursday_11_May_2023_10:16:5_ACST
org.label-schema.schema-version: 1.0
org.label-schema.usage.singularity.deffile.bootstrap: docker
org.label-schema.usage.singularity.deffile.from: quay.io/vgteam/vg:v1.48.0
org.label-schema.usage.singularity.version: 3.8.6
```

However, it is necessary to bind the singularity container with the VCF files. 
By default Singularity bind mounts ```/home/$USER```
Assume the files are in the directory ```data``` :
```sh
singularity shell --bind data:/mnt image.sif
```
This will open a temporary sandbox, that you can exit at anytime.
If needed, it is possible to run several instances of the image at the same time, see https://docs.sylabs.io/guides/3.0/user-guide/running_services.html


Now vg can be launched :
```sh
vg help 
```
Note that before starting to build a graph, an index of the VCF and FASTA files are required. 

The VCF index file can be generated using the tabix command provided by SAMtools : ``` tabix -p vcf your-vcf-file.vcf.gz``` 

The FASTA .gz file needs to be converted to a .bgz before indexing, like so : ```zcat reference.fna.gz | bgzip -c > reference.fna.bgz```
And then indexed :  ```samtools faidx  reference.fna.bgz```

To start building a graph, still in the ```data``` directory, with your ```reference.fna.bgz``` file in FASTA and ```your-vcf-file.vcf.gz```, the following will construct a graph in ```x.vg``` :

```sh
vg construct -r data/reference.fna.bgz -v data/your-vcf-file.vcf.gz >x.vg
```
With several VCF files, either merge them into one file using SAMtools or run your command with the following syntax :

```sh
vg construct -r data/reference.fa -v data/your-vcf-file1.vcf.gz -v data/your-vcf-file2.vcf.gz >x.vg
```

Then, you can get the reference, the VCF and the index. 
```sh
# get the reference
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# get the HGDP vcfs by changing the command with the correct files found here depending on the requested chromosome : https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/
 wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz
```


# Creating Subsets for the target population
The goal now is to filter the VCF files by keeping the individuals from the target population only. This can be done with bcftools, a SAMtools project.

There are 2 criterion to make these subsets : the population and the Minor Allele Frequency (MAF). 
The targets are Mozambite individuals and MAF > 0.01%.

First, let's try it out on a single vcf file, here chromosome 21.
Here is the command line to exclusively keep the Mozambite population and put it in a new and compressed vcf file :
```sh
cd data
bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  hgdp_21.vcf.gz > gzip > moza21.vcf.gz
```
When creating this subset with this population, one individual is apparently not in the VCF file but was listed in HGDP documentation 
```sh
Warn: subset called for sample that does not exist in header: "HGDP01273"... skipping
```

The reference for the population IDs can be found here : https://www.internationalgenome.org/data-portal/population/MozabiteHGDP 

# Pangenome trial run
Let's try and build a pangenome graph using chromosome 21 for the Mozambite population. 
Let's start by indexing the reference and the vcf files :
```sh
# index the ref file
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

# zip and index the vcf 
bgzip -c data/moza21.vcf.gz > data/moza21.vcf.gz                                                                                                                                                                   
tabix -p vcf data/moza21.vcf.gz 

# run singularity
singularity shell --bind data:/mnt image.sif
vg construct -r data/GRCh38_full_analysis_set_plus_decoy_hla.fa -v data/moza21.vcf.gz
# error: no VCF header
```
There is an error, apparently with the VCF file, probably bc of the zipping process. Once the file is loaded, it has a header amd all the information we need. However, it cannot be indexed (?) so the solution is to zip it in .gz, but this seems to delete the content of the file.








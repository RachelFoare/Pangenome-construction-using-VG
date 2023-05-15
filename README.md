# Pangenome-construction-using-VG
Steps to the creation of a pangenome in VG using Singularity -- CSIRO INTERNSHIP 2023

The goal is to understand how singularity and vg work, in order to build a Pangenome.
In Somalis, this is the 3rd highest cause of death per year (4.57% as per 2020 WHO statistics). What's more, when it comes to genetic studies, Somalis (and East Africa in general) is severely understudied. In fact there is not one single paper with East African whole genomes. Because the Somalis have unique indigenous ancestral components that are distinct from the rest of Africa, let alone the rest of the world (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004393#pgen.1004393-Pagani1), most whole genomes for which data is available (including Southern and Western Africa) will not be reliable for constructing a pangenome. In fact the Fst values between many of these Northern, Southern and Western Africans is as high (in some cases higher) than what you see between Europeans and Asians.
The only exception is of course, the North African Mozabites, who have the smallest Fst and highest ancestry contributions to Somali populations amongst all source populations outside of East Africa (22% on average, the best fitting model according to qpAdm), likely due to shared drift and past admixture events. The Bedouin act as a good control (as they have the lowest Fst, ~0.01 and the most ancestral affinity to the Mozabite), to evaluate the performance of the Pangenome that will be constructed.


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
This will open a temporary sandbox, that you can ```exit``` at anytime.
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
Here is the command line to exclusively keep the target population and MAF and put it in a new and compressed vcf file :
```sh
cd data
bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_21.vcf.gz > gzip > moza21.vcf.gz
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
cd
# index the ref file
samtools faidx data/GRCh38_full_analysis_set_plus_decoy_hla.fa

# index the vcf using tabix                                                                                                                                                                  
tabix -p vcf data/moza21.vcf.gz 
# or using bcftools
bcftools index moza21.vcf.gz 

# run singularity
singularity shell --bind data:/mnt image.sif
vg construct -r data/GRCh38_full_analysis_set_plus_decoy_hla.fa -v data/moza21.vcf.gz
```

![Alt text](http://full/path/to/img.jpg "Population-Specific Pangenome Graph of Chromosome 21 for Mozabite People with MAF>0.01
")


# Using the HPC
The HPC uses Slurm.
Here are some basic but helpful commands :

Running a job : ```sbatch jobscript-name.sh```
Seeing the queue : ```squeue```
Seeing details of the job using its ID : ```scontrol show job jobID```
Deleting the job using its ID : ```scancel jobID```

Changing the HPC resource configs :
```sh
#SBATCH --job-name=HelloWorld
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=500m
#SBATCH --account=OD-012345
```

# Creating a Pangenome on the HPC
Let's start with writting a complete .sh file that the HPC will be able to run.

It can be useful to use 2 jobscripts, since all the vcfs need to run in the same singularity container. Let's write a script that will write another script containing the commands needed to run vg in singularity.
Open a note pad, e.g. ```micro```.
In the notepad, specify the resources needed for the job (i.e. job name, wall time, nodes, number of tasks, cpu, RAM, etc).

Let's do a dummy script called ```helloworld_script.sh``` to test that, without specified resources :
```sh
#!/bin/bash
echo "Hello, world!"
```

And run it : 
```sh
sbatch helloworld_script.sh
```

Printing HelloWorld is quite easy, now let's write a script called ```helloworld-advanced_script.sh```that will create and execute another script :
```sh
#!/bin/bash
cat > my-2nd-script.sh <<EOF
#!/bin/bash
echo "Hello, world!"
EOF

chmod +x my-2nd-script.sh
./my-2nd-script.sh
```

And run it :
```sh
sbatch helloworld-advanced_script.sh
``` 
Let's download all the VCF files, the reference file and the container image in a new directory called ```data``` :

```sh
cd
mkdir data
cd data
singularity build image.sif docker://quay.io/vgteam/vg:v1.48.0 
singularity shell --bind data:/mnt image.sif

# get the reference
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# get the HGDP vcfs by changing the command : https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr1.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr2.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr3.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr4.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr5.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr6.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr7.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr8.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr9.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr10.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr11.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr12.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr13.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr14.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr15.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr16.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr17.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr18.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr19.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr20.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr22.vcf.gz
```

Let's do the subsets for each VCF file :
```sh
bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data1.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr2.vcf.gz > gzip > data2.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr3.vcf.gz > gzip > data3.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz

bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_wgs.20190516.full.chr1.vcf.gz > gzip > data21.vcf.gz
```


And now with our actual commands to create the graph, considering all the VCFs and the reference files are in ```data``` and have not been indexed yet :
```sh
#!/bin/bash
#SBATCH --job-name=myjob
#SBATCH --output=myjob_output.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --partition=your_partition_name

# go to the correct directory
cd
# form the subsets of Mozabite individuals with MAF>0.01

# index all the files

# load singularity
module load singularity
# create a script file with the commands to run in the same singularity container
cat > sing-commands.sh <<EOF
#!/bin/bash

EOF
chmod +x sing-commands.sh

# Run the Singularity container with the script file as an argument
singularity exec data/image.sif bash sing-commands.sh
```
Now, let's run the job.

```sh
sbatch pangenome_script.sh
```


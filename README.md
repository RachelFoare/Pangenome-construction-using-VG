# Pangenome-construction-using-VG  
Steps to the creation of a pangenome in VG using Singularity -- CSIRO INTERNSHIP 2023

Note : The pangenome will be constructed with autosomal chromosomes only. The data and values are cited from an Ongoing Research by Johar A. et al from Mayo Clinic Laboratory, 2023

The goal is to understand how singularity and vg work, in order to build a Pangenome.
In Somalis, CHD is the 3rd highest cause of death per year. What's more, when it comes to genetic studies, Somalis (and East Africa in general) is severely understudied. [Because the Somalis have unique indigenous ancestral components that are distinct from the rest of Africa, let alone the rest of the world](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004393#pgen.1004393-Pagani1), most whole genomes for which data is available (including Southern and Western Africa) will not be reliable for constructing a pangenome. The Fst values between many of these Northern, Southern and Western Africans is at least as high as what you see between Europeans and Asians.
The only exception is the North African Mozabites, who have the smallest Fst and highest ancestry contributions to Somali populations amongst all source populations outside of East Africa (22% on average, the best fitting model according to qpAdm), likely due to shared drift and past admixture events. The Bedouin act as a good control (as they have the lowest Fst, ~0.01 and the most ancestral affinity to the Mozabite), to evaluate the performance of the Pangenome that will be constructed. 




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
This will open a temporary sandbox, that you can ```exit``` at anytime. If bound correctly, it is possible to naviguate through the existing directories.
If needed, it is possible to run several instances of the image at the same time, see [here](https://docs.sylabs.io/guides/3.0/user-guide/running_services.html).



Then, you can get the reference, the VCF and the index for each file. 
```sh
# get the reference
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# get the HGDP vcfs by changing the command with the correct files found here depending on the requested chromosome : https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/
 wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz
```

Now vg can be launched :
```sh
vg help 
```
Note that before starting to build a graph, an index of the VCF and FASTA files are required, referenced in the following pipeline as ```.fai``` and ```.tbi``` :

![Alt text](https://github.com/vgteam/vg/raw/master/doc/figures/pipeline_flowchart.svg "VG Pipeline")



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


# Creating Subsets for the target population

The goal now is to filter the VCF files by keeping the individuals from the target population only. This can be done with bcftools, a SAMtools project.

There are 2 criterion to make these subsets : the population and the Minor Allele Frequency (MAF). 
The targets are Mozabite individuals and MAF > 0.01%.

First, let's try it out on a single vcf file, here chromosome 21.
Here is the command line to exclusively keep the target population, called by their IDs (e.g. HGDP01275) and MAF and put it in a new and compressed vcf file :
```sh
cd data
bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' hgdp_21.vcf.gz > bgzip > sub21.vcf.gz
```
When creating this subset with this population, one individual is apparently not in the VCF file but was listed in HGDP documentation 
```sh
Warn: subset called for sample that does not exist in header: "HGDP01273"... skipping
```

The reference for the population IDs can be found here : https://www.internationalgenome.org/data-portal/population/MozabiteHGDP 




# Using the HPC
The HPC uses Slurm.
Here are some basic but helpful commands :

Running a job : ```sbatch jobscript-name.sh```

Seeing the queue : ```squeue```

Seeing details of the job using its ID : ```scontrol show job jobID```

Seeing how resource efficient the job is : ```seffx jobID.batch```

Deleting the job using its ID : ```scancel jobID```

Seeing the differennt partitions : ```sinfo```

Changing the HPC resource configs :
```sh
#SBATCH --job-name=HelloWorld
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=500m
#SBATCH --account=OD-012345
```

Access the HPC using PuTTY to have an SSH connection (note : this part was done on Windows)

Let's start with writting a complete .sh file that the HPC will be able to run.

Open a note pad, e.g. ```nano```.
In the notepad, specify the resources needed for the job (i.e. job name, wall time, nodes, number of tasks, cpu, RAM, etc). Note : optional.

Let's do a dummy script called ```helloworld_script.sh``` to test that, without specified resources :
```sh
#!/bin/bash
echo "Hello, world!"
```

And run it : 
```sh
 sbatch -A OD-221017 helloworld_script.sh
 # the HPC required to specify a project ID
 # to check the status of the jobs 
 squeue -u username
 # to see the output
 less slurm-jobID.out
```

Before going to the next section, check if the singularity module is installed on the HPC :
```sh
 module avail
```

# Accessing CSIRO's HPC -- data manager 
The documentation about accessing the data manager to download files with a large amount of data, for CSIRO's HPC can be found [here](https://confluence.csiro.au/display/SC/CSIRO+SC+Shared+Cluster+-+Petrichor).
It is important to make sure the data is stored in the right place so i dooesn't get 'flushed' (i.e. deleted).
Preferably, access the petrichor h4 petrichor-login-h4.hpc.csiro.au  (giving access to 4TB nodes) with PuTTY or an ssh session on a linux terminal.
Otherwise, access the petrichor-dm interactive session using this command : ```sinteractive -p io -A JOBID  -t 2:00:00```
Access the user's datastore directory using : ```cd /datastore/username```


If possible, transfer your files on the HPC. Otherwise, download the image and the files again (cf. next step).
Note : When building the image in Singularity with ```singularity inspect image.sif```, checking the image returns the following :
```sh
org.label-schema.build-arch: amd64
org.label-schema.build-date: Wednesday_17_May_2023_14:4:21_AEST
org.label-schema.schema-version: 1.0
org.label-schema.usage.singularity.deffile.bootstrap: docker
org.label-schema.usage.singularity.deffile.from: quay.io/vgteam/vg:v1.48.0
org.label-schema.usage.singularity.version: 3.7.3

```
# Creating a Pangenome on the HPC

File system convention :
- the files containing the data are stored in ```/datastore/username/data```
- the job scripts are created and ran in ```/scratch3/username```

It can be useful to use 2 jobscripts, since all the vcfs need to run in the same singularity container. It is possible to write a script that will write another script containing the commands needed to run vg in singularity. Otherwise, a simple loop can do the trick. 
Let's download all the VCF files, the reference file and the container image in a new directory called ```data``` :

Note : All the data has actually been retrieved using a batch job, using the ```wget``` commands in a shell script, so that the downloads could run effectively and in parallel.  

```sh
cd /datastore/username
mkdir data
cd data
# get the image
module load singularity
singularity build image.sif docker://quay.io/vgteam/vg:v1.48.0 
# get the reference
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa && mv GRCh38_full_analysis_set_plus_decoy_hla.fa ref.fa

# get the HGDP vcfs by changing the command : https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/
# and rename the files as chr#.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr1.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr2.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr3.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr4.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr5.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr6.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr7.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr8.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr9.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr10.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr11.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr12.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr13.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr14.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr15.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr16.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr17.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr18.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr19.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr20.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr22.vcf
```

Create the subsets.
It is possible to do the subsets with each VCF file in a loop :
```sh
subset() {
 chr=$1
 vcfLIST=$(while read iid; do echo "-v ./"chr${iid}".vcf.gz"; done )
 cd data
 module load bcftools
 bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' chr${chr}.vcf.gz > bgzip > sub-${chr}.vcf.gz
}
```
It can also be done by simply changing the name of the file as follows :
```sh
module load bcftools
bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' chr1.vcf.gz > bgzip > sub-chr1.vcf.gz
bcftools view --force-samples -s HGDP01275,HGDP01282,HGDP01256,HGDP01263,HGDP01268,HGDP01270,HGDP01276,HGDP01257,HGDP01264,HGDP01272,HGDP01277,HGDP01258,HGDP01260,HGDP01265,HGDP01254,HGDP01259,HGDP01261,HGDP01266,HGDP01273,HGDP01280,HGDP01255,HGDP01262,HGDP01267,HGDP01279,HGDP01274  -i 'MAF > 0.01' chr2.vcf.gz > bgzip > sub-chr2.vcf.gz
```

It is possible that the HPC does not have the required tools such as bgzip or tabix to index the files in their module (look for  ```htslib```), so make sure it is, otherwise install conda, then pbgzip and tabix from bioconda. If ```htslib``` is available, simply load the module and use ```bgzip```.
```sh
 wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.3.1-0-Linux-x86_64.sh
 conda create -n environment
 conda activate environment
 conda install -c bioconda pbgzip
 conda install -c bioconda tabix
```
Then use the ```pbgzip``` command and the ```tabix``` command to index all the vcfs, and use ```samtools faidx``` to index the reference.

```sh
 module load samtools
 samtools faidx ref.fa
 
 # and run a job containing the following commands for each vcf file :
 pbgzip sub-chr21.vcf
 tabix -p vcf sub-chr21.vcf.gz
 
 # if you have trouble using samtools faidx, try zipping the ref.fa file with pbgzip and then unzip it and try again
 pbgzip ref.fa
 gunzip ref.fa.gz
 samtools faidx ref.fa
```

And now with our actual commands to create the graph :
Note : these commands are in a script called ```pangenome_script.sh``` and for chromosome 21 specifically to test it.
```sh
#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G


cd /datastore/username
module load singularity

singularity exec --bind /datastore/user/data:/datastore/user/data image.sif vg construct -r /datastore/user/data/ref.fa -v /datastore/user/data/sub-chr21.vcf.gz > /datastore/user/data/p21.vg
```
Create all the graphs using this command, either in a loop or by repeating the singularity command in the script like so :
```sh
singularity exec --bind /datastore/username/data:/datastore/username/data image.sif vg construct -r /datastore/username/data/ref.fa -v /datastore/username/data/sub-chr9.vcf.gz >/datastore/username/data/p9.vg
singularity exec --bind /datastore/username/data:/datastore/username/data image.sif vg construct -r /datastore/username/data/ref.fa -v /datastore/username/data/sub-chr8.vcf.gz >/datastore/username/data/p8.vg
singularity exec --bind /datastore/username/data:/datastore/username/data image.sif vg construct -r /datastore/username/data/ref.fa -v /datastore/username/data/sub-chr7.vcf.gz >/datastore/username/data/p7.vg
```
Then the nodes need to be coordinated :

```sh
module load singularity
singularity exec --bind /datastore/username/data:/datastore/username/data image.sif vg ids -j $(for i in $(seq 1 22); do echo p$i.vg; done)
```

And then indexed as ```wgx.xg``` **and** ```wg.gcsa``` files. 

The following lines index the pangenomes as ```wgx.xg```. Note that this batch was run with 1TB of RAM, but 1.5TB is ideal. It can take several hours.
```sh
#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=1000G

cd /datastore/user
module load singularity


singularity exec --bind /datastore/user/data:/datastore/foa003/data image.sif vg index -j /datastore/user/data/snarls.txt -b /datastore/user/data -p /datastore/user/data/p1.vg /datastore/user/data/p2.vg /datastore/user/data/p3.vg /datastore/user/data/p4.vg /datastore/user/data/p5.vg /datastore/user/data/p6.vg /datastore/user/data/p7.vg /datastore/user/data/p8.vg /datastore/user/data/p9.vg /datastore/user/data/p10.vg /datastore/user/data/p11.vg /datastore/user/data/p12.vg /datastore/user/data/p13.vg /datastore/user/data/p14.vg /datastore/user/data/p15.vg /datastore/user/data/p16.vg /datastore/user/data/p17.vg /datastore/user/data/p18.vg /datastore/user/data/p19.vg /datastore/user/data/p20.vg /datastore/user/data/p21.vg /datastore/user/data/p22.vg

```

Before indexing as ```wg.gcsa``` it is necessary to prune the graphs, i.e. to mask out highly complex regions.

```sh
#SBATCH --array=1-22

export APPTAINER_CACHEDIR=/scratch3/user/singularity_cache
export S_IMG=/scratch3/user/image.sif
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
singularity exec ${S_IMG} vg prune -t ${SLURM_CPUS_PER_TASK} -k 45 -r /scratch3/user/data/p${SLURM_ARRAY_TASK_ID}.vg > /scratch3/user/data/pr${SLURM_ARRAY_TASK_ID}.vg
```

To run the following job properly, make sure to add ```-z 20480``` to the sbatch command, and make sure the temp storage is higher than 20TB. 

Note that constructing this index takes 7 days to construct.
```sh
singularity exec ${S_IMG} vg index -t ${SLURM_CPUS_PER_TASK} --temp-dir /scratch3/foa003/temp -g /scratch3/foa003/wg.gcsa /scratch3/foa003/data/pr{1..22}.vg -p -Z 20480
```

# Assessing the Pangenome - Statistics for the internship report

In order to properly assess the performance of the graph, it is essential to examine criteria reflecting the quality of the alignment of the graph with the reads of an external group remaining genetically close, and compare it with the same alignment with the original linear reference used to build the pangenome. 
These criteria can be found and measured using vg, the ```Picard``` module and SAMtools.
The stats will be calculated by the results of the alignment of Bedouin reads to the graph, and an alignment of the same individuals on the linear reference.
The linear alignments can be found [here](https://www.internationalgenome.org/data-portal/sample), with the specified population (Bedouin for HGDP or Bedouin B for SGDP).
However in this case it was decided to redo the alignments in order to keep track of the methodology using ```bowtie 2.5.1``` AND ```bwa 0.7.17-r1188``` using the FASTA files of the reads.
Due to lack of time and storage, 5 sequences of Bedouin people from HGDP and 2 from SGDP are used for the alignments. The 5 samples from HGDP were picked randomly. 
Let's start by downloading the FASTA file for this individual (HGDP00616), and indexing the reference file :

```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR141/002/ERR1419132/ERR1419132_2.fastq.gz
module load bowtie
bowtie2-build ref.fa ref-index
module load bwa
bwa index ref.fa 
```
The samples used for the following statistics are the following :
from SGDP
```sh
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR141/001/ERR1419141/ERR1419141_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR141/002/ERR1419132/ERR1419132_2.fastq.gz
```
from HGDP
```sh
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/005/ERR1451355/ERR1451355_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/000/ERR1451680/ERR1451680_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/001/ERR1451031/ERR1451031_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR145/009/ERR1451329/ERR1451329_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR144/000/ERR1443460/ERR1443460_2.fastq.gz
```

 The sequences will be aligned using both bwa and bowtie, depending on their database of origin. The resulting alignments sould be :
 - bowtie alignment with the 5 HGDP samples, built with ```bowtie2 -q -x ref-index -U B1.fastq.gz,B2.fastq.gz,B3.fastq.gz,B4.fastq.gz,B5.fastq.gz -S output-bow-5.sam```
 - bwa alignment with the 5 HGDP samples ```bwa mem ref.fa <(cat B1.fastq.gz B2.fastq.gz B3.fastq.gz B4.fastq.gz B5.fastq.gz) > output-bwa-5.sam```
 - bowtie alignment with the 2 SGDP samples, built with ```bowtie2 -q -x ref-index -U B6.fastq.gz,B7.fastq.gz -S output-bow-7.sam```
 - bwa alignment with the 2 SGDP samples ```bwa mem ref.fa <(cat B6.fastq.gz B7.fastq.gz) > output-bwa-7.sam```

These alignments are assessed using ```vg stats``` and ```picard CollectAlignmentSummaryMetrics``` in order to have information about indels, MQ0, mismatches and soft clipping. 
```sh
samtools stat out7.sam > statbwa7.txt
picard CollectAlignmentSummaryMetrics -I out7.sam -O picbwa7.txt
```

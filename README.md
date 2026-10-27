<h1 align="center">

PLATE: Plasmid Atlas of Food Ecosystems

</h1>

If you use this repository, please cite:

Author list, PLATE: Plasmid Atlas of Food Ecosystems, *bioRxiv*, 2026
[Link].

------------------------------------------------------------------------

# Table of contents

-   [Overview](#overview)
-   [Step-by-Step Workflow](#step-by-step-workflow)
-   [References](#references)

------------------------------------------------------------------------

# Overview

**PLATE** is a comprehensive resource of **food-related plasmidomes**,
built with an end-to-end pipeline that retrieves metagenomic assemblies,
identifies plasmid sequences and annotates them.

The PLATE datasets are publicly available through the links below:

| Catalog | Catalog link |
|------------------------------------|------------------------------------|
| PLATE - Plasmids | [Link](https://plate-db.com/downloads) |
| PLATE - Metadata | [Link](https://plate-db.com/samples) |

PLATE is also accessible as a searchable catalog at
<https://plate-db.com/>

------------------------------------------------------------------------

# Step-by-Step Workflow

PLATE was constructed by:

1)  Retrieving **metagenomic assemblies** (or assembling them),
2)  Predicting **plasmid contigs** with Plasmer,
3)  Removing **identical contigs** with seqkit and grouping the unique
    contigs into **plasmid-like clusters** with MMseqs2,
4)  Predicting **genes and smORFs** on plasmid contigs,
5)  Assigning **hosts** to plasmid contigs,
6)  Running **functional annotation** tools (protein families, resistance
    genes, CAZymes, CRISPR-Cas, AMPs and ACPs),
7)  Constructing a **refined plasmid set** by passing the Plasmer-predicted
    plasmid contigs through five successive filters (PlasClass,
    viralVerify, PlasFlow, Infernal and minimap2).

------------------------------------------------------------------------

## Part 1 — Data acquisition

### LOGAN

Metagenomic contigs are retrieved from the LoganSRA repository. The
assemblies are compressed and must be downloaded and extracted.

#### Usage

``` bash
wget https://s3.amazonaws.com/logan-pub/c/$accession/$accession.contigs.fa.zst
zstd -d $accession.contigs.fa.zst
```

#### Notes

-   accession: The SRA run accession ID (e.g., ERR123456).
-   zstd -d: Decompresses the `.zst` file into a standard FASTA file
    (`.fa`).

### Assembly (optional)

Required only if your data is not available via LOGAN. Unitigs are
first built with Cuttlefish2, then assembled into contigs with Minia3.

#### Cuttlefish2

##### Setup

``` bash
conda create -n cuttlefish -c bioconda -c conda-forge -c defaults python=3.10 sra-tools seqkit "cmake>=3.16" make gxx_linux-64 zlib git -y
conda activate cuttlefish

mkdir -p vendor
cd vendor
git clone https://github.com/rchikhi/cuttlefish.git
cd cuttlefish
git checkout 9401ef5

mkdir -p build
cd build
export CMAKE_POLICY_VERSION_MINIMUM=3.5
cmake -DCMAKE_INSTALL_PREFIX=../ \
      -DPYTHON_EXECUTABLE="$(which python3)" \
      ..
make -j "$(nproc)" install

mkdir -p bin
ln -sf vendor/cuttlefish/bin/cuttlefish bin/cuttlefish
```

##### Usage

``` bash
cuttlefish build --ref -s input.fa -k 31 -c 2 -t 8 -o output.unitig.fa
```

##### Notes

-   --ref : Builds a compacted reference de Bruijn graph (for FASTA
    input)
-   -s : Input FASTA file
-   -k : k-mer length
-   -c : Minimum k-mer abundance cutoff
-   -t : Number of threads
-   -o : Output file (unitigs)

#### Minia3

##### Setup

``` bash
conda create -n minia3 python=3.10 -y
conda activate minia3

git clone https://github.com/GATB/minia.git
cd minia
git checkout 71484e8
git submodule update --init --recursive

mkdir build && cd build
CMAKE_POLICY_VERSION_MINIMUM=3.5 cmake ..
make -j$(nproc)
```

##### Usage

``` bash
minia \
  -in input.unitig.fa \
  -kmer-size 31 \
  -abundance-min 1 \
  -traversal contig \
  -out output
```

##### Notes

-   -in : Input file (unitigs from Cuttlefish2)
-   -kmer-size : k-mer length used for assembly
-   -abundance-min : Minimum k-mer abundance threshold
-   -traversal : Graph traversal mode (builds contigs from the unitig
    graph)
-   -out : Output prefix

------------------------------------------------------------------------

## Part 2 — Plasmid prediction

Plasmer is run on the metagenomic contigs from Part 1. Its predicted
plasmid contigs (`predPlasmid.fa`) are the input of Parts 3–7.

### Plasmer

#### Setup

``` bash
conda create -n plasmer -c bioconda -c conda-forge
conda activate plasmer

conda install -y seqkit=2.2.0 python=3.10.4 biopython prodigal=2.6.3 hmmer=3.3.2 "blast>=2.10.1" infernal=1.1.4 diamond=2.0.8 parallel=20220722 kraken2=2.1.2 r-base=4.2.0 r-randomforest
conda install -y -c bioconda kmer-db=1.9.4

git clone https://github.com/nekokoe/Plasmer.git
cd Plasmer
export PATH=$PATH:$(pwd)
```

#### Install the databases

The pre-built database are available at [Zenodo](https://zenodo.org/records/7030675) and [Google Drive](https://drive.google.com/drive/folders/1C1XLqTuPVfwyB7N3b5QAU09Rn-Qq7RdD?usp=sharing).

Verify the downloaded databases:

``` bash
$ sha1sum plasmerMainDB.tar.xz 
0b08f5c30d60b137f54de6024ab7557031850db6  plasmerMainDB.tar.xz

$ sha1sum customizedKraken2DB.tar.xz 
b14efdd9232fd5f6d066716bd8e3e6ca80c9c0de  customizedKraken2DB.tar.xz
```

Extract the databases:

``` bash
mkdir -p db_Plasmer
tar -xJvf plasmerMainDB.tar.xz -C db_Plasmer
tar -xJvf customizedKraken2DB.tar.xz -C db_Plasmer
```

#### Usage

``` bash
./Plasmer/scripts/Plasmer -g input.fa -p $accession -d db_Plasmer -t 3 -m 500 -l 500000 -o outpath
```

#### Notes

-   -g : Input FASTA file (contigs)
-   -p : Output file prefix
-   -d : Path to the Plasmer database
-   -t : Number of CPU threads
-   -m : Minimum contig length (bp)
-   -l : Maximum contig length (bp)
-   -o : Output directory. Predictions are saved as `.predClass.tsv` and separate FASTA files.

------------------------------------------------------------------------

## Part 3 — Dereplication and clustering

Identical contigs are first collapsed with seqkit, then the unique
contigs are grouped into plasmid-like clusters with MMseqs2.

### seqkit

Removes identical contigs from the Plasmer-predicted plasmid contigs
and from the refined plasmid set (Part 7). **Contigs with the same
sequence on either strand (reverse complements count as identical) are
collapsed into one; the first record is kept.**

#### Setup

``` bash
conda create -n seqkit -c conda-forge -c bioconda seqkit
conda activate seqkit
```

#### Usage

``` bash
seqkit rmdup -s -j 4 -D duplicates.txt input.fa > output.unique.fa
```

#### Notes

-   rmdup : Removes duplicated sequences
-   -s : Compares sequences instead of names; both strands are compared
-   -j : Number of threads
-   -D : Output file of the duplicated contigs (number of copies and
    their names)
-   output.unique.fa : Output FASTA file with one copy of every sequence

### MMseqs2

Clusters the unique contigs of the refined plasmid set (seqkit output)
into plasmid-like clusters. **A contig joins a cluster when it is ≥ 95%
identical to the cluster representative and the alignment covers ≥ 85%
of the contig; the longest contig of each cluster is its
representative.**

#### Setup

``` bash
conda create -n mmseqs2 -c conda-forge -c bioconda mmseqs2
conda activate mmseqs2
```

#### Usage

``` bash
mmseqs easy-cluster input.fa output_prefix tmp_dir \
    --min-seq-id 0.95 \
    -c 0.85 \
    --cov-mode 1 \
    --cluster-mode 2 \
    --alignment-mode 3 \
    --kmer-per-seq-scale 0.2 \
    --max-seq-len 500000 \
    --threads 4 \
    --split-memory-limit 9G
```

#### Notes

-   output_prefix : Output prefix (`_cluster.tsv`, `_rep_seq.fasta`,
    `_all_seqs.fasta`)
-   tmp_dir : Temporary directory
-   --min-seq-id : Minimum sequence identity to the cluster
    representative
-   -c : Minimum alignment coverage
-   --cov-mode : Coverage is measured on the cluster member
-   --cluster-mode : Greedy incremental clustering by length (the
    longest contig represents its cluster)
-   --alignment-mode : Computes the alignment score, coverage and
    sequence identity
-   --kmer-per-seq-scale : Number of k-mers per sequence, scaled with
    sequence length
-   --max-seq-len : Maximum sequence length (the default, 65,535 bp, is
    shorter than the longest contigs)
-   --threads : Number of threads
-   --split-memory-limit : Maximum memory used per database split

------------------------------------------------------------------------

## Part 4 — Gene prediction

Prodigal and smORFinder are run in meta mode on the predicted plasmid
contigs (`predPlasmid.fa`, the output of Plasmer).

### Prodigal

#### Setup

``` bash
git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal/
make install
```

#### Usage

``` bash
prodigal -p meta \
         -i input.fa \
         -o output.gff \
         -a output.faa \
         -d output.fna \
         -f gff
```

#### Notes

-   -p : Mode (single/meta)
-   -i : Input FASTA file
-   -o : Gene coordinates (.gff)
-   -a : Predicted protein sequences (.faa)
-   -d : Predicted gene nucleotide sequences (.fna)
-   -f : Output format

### smORFinder

#### Setup

``` bash
conda create -n smorfinder python=3.8.20
conda activate smorfinder

pip install smorfinder
```

#### Usage

``` bash
smorf meta input.fa -o outpath -t $threads
```

#### Notes

-   meta : Runs smORFinder in metagenomic mode
-   input.fa : Input FASTA file (plasmid contigs)
-   -o : Output directory
-   -t : Number of threads

------------------------------------------------------------------------

## Part 5 — Host assignment

HOTSPOT is run on the predicted plasmid contigs (`predPlasmid.fa`).

### HOTSPOT

#### Setup

``` bash
git clone https://github.com/Orin-beep/HOTSPOT
cd HOTSPOT/
conda env create -f environment.yaml -n hotspot
conda activate hotspot
```

Install the databases:

``` bash
sh prepare_db.sh        # download and unzip the database folder
sh prepare_mdl.sh       # download and unzip the model folder
```

#### Usage

``` bash
python preprocessing.py \
        --fasta input.fa \
        --database folder/of/database \
        --model_path folder/of/model \
        --threads $threads
```

``` bash
python hotspot.py --threads $threads --out outpath --device cpu 
```

#### Notes

-   preprocessing.py : Prepares the input sequences for prediction
-   --fasta : Input FASTA file (plasmid contigs)
-   --database : Path to the downloaded database folder (PC proteins,
    MOB/MPF proteins and replicons)
-   --model_path : Path to the downloaded model folder
-   --threads : Number of threads
-   hotspot.py : Predicts the hosts of the preprocessed sequences
-   --out : Output directory for the prediction results
-   --device : Device used for prediction ('gpu' or 'cpu')

------------------------------------------------------------------------

## Part 6 — Functional annotation

Protein-level tools use the Prodigal proteins (`.faa`), Macrel and
MetaPepticon use the smORFinder peptides, and dbCAN (CGC mode) and
CRISPR-Cas Typer use the plasmid contigs directly.

### eggNOG-mapper

Identifies protein family profiles of plasmid contigs.

#### Setup

``` bash
conda create -n emapper -c bioconda -c conda-forge python=3.10.19
conda activate emapper
conda install eggnog-mapper
```

Install the databases:

``` bash
aria2c -c -x 8 -s 8 -k 1M http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
aria2c -c -x 8 -s 8 -k 1M http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
aria2c -c -x 8 -s 8 -k 1M http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
```

#### Usage

``` bash
emapper.py \
    -i input.faa \
    -o output_prefix \
    --output_dir output_dir \
    --data_dir eggnog_db_dir \
    --scratch_dir tmp_dir \
    --cpu 8 \
    -m diamond \
    --override \
    --dmnd_ignore_warnings
```

#### Notes

-   -i : Input protein FASTA file
-   -o : Output file prefix
-   --output_dir : Output directory
-   --data_dir : Path to EggNOG database files (eggnog.db,
    eggnog_proteins.dmnd, eggnog.taxa.db)
-   --scratch_dir : Temporary directory
-   --cpu : Number of CPU threads
-   -m diamond : Uses Diamond for homology search
-   --override : Overwrites existing output files
-   --dmnd_ignore_warnings : Tolerates ambiguous amino acids (e.g., X)

### AMRFinderPlus

Detects resistance genes (antimicrobial, stress, heat, heavy-metal,
etc.).

#### Setup

``` bash
conda create -y -c conda-forge -c bioconda -n amrfinder --strict-channel-priority ncbi-amrfinderplus
conda activate amrfinder
```

#### Usage

``` bash
amrfinder -p input.faa --plus --threads $threads -o output.tsv
```

#### Notes

-   -p : Protein FASTA file
-   --plus : Adds "Plus" genes (virulence factors, stress-response
    genes, etc.)
-   -o : Output file path
-   --threads : Number of threads

### Resistance Gene Identifier

Detects antimicrobial resistance (AMR) genes using the CARD database.

#### Setup

``` bash
conda create --name rgi -c conda-forge -c bioconda -c defaults 
conda activate rgi
conda install -c conda-forge -c bioconda -c defaults rgi 
```

CARD database setup (system-wide):

``` bash
rgi load --card_json /path/to/card.json 
```

CARD database setup (local/working directory):

``` bash
rgi load --card_json /path/to/card.json --local
```

#### Usage

``` bash
rgi main -a DIAMOND \
         -i input.faa \
         -o outpath \
         --local \
         --clean \
         -t protein \
         -n $threads
```

#### Notes

-   load : Loads the CARD database
-   --card_json : Path to the CARD database JSON file
-   main : Runs AMR gene detection
-   -a : Alignment tool
-   -i : Input file path
-   -o : Output file path
-   --local : Uses the local database
-   --clean : Removes temporary files
-   -t : Input type (contig or protein)
-   -n : Number of threads

### dbCAN

Identifies CAZymes and CAZyme gene clusters.

#### Setup

``` bash
git clone https://github.com/bcb-unl/run_dbcan.git
cd run_dbcan
conda create -f envs/environment.yml
conda activate run_dbcan_env
pip install dbcan
```

Install the databases:

``` bash
mkdir -p db
run_dbcan database --db_dir db
```

#### Usage

``` bash
# for CAZyme Annotations
run_dbcan CAZyme_annotation \       
    --mode protein \
    --input_raw_data input.faa \
    --db_dir db \
    --output_dir outpath
```

``` bash
# for CAZyme Gene Clusters
run_dbcan easy_CGC \                
   --mode meta \
   --input_raw_data input.fasta \
   --db_dir db \
   --output_dir outpath
```

#### Notes

-   --mode : Input type. `protein` for CAZyme_annotation, `meta` (plasmid
    contigs) for easy_CGC
-   --input_raw_data : Input file (.faa or contig FASTA)
-   --db_dir : Database directory
-   --output_dir : Output directory

### CRISPR-Cas Typer

Detects CRISPR-Cas genes and arrays.

#### Setup

``` bash
conda create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper
conda activate cctyper
pip install "setuptools<70"
```

#### Usage

``` bash
cctyper input.fasta output_dir --prodigal meta
```

#### Notes

-   --prodigal : Prodigal mode used by CCTyper (single/meta)

### Macrel

Identifies antimicrobial peptides (AMPs).

#### Setup

``` bash
conda create -n macrel python=3.10 -y
conda activate macrel
conda install -y -c bioconda -c conda-forge megahit paladin ngless
pip install macrel
```

#### Usage

``` bash
macrel peptides \
  --fasta input.faa \
  --output output_dir \
  --threads 8 
```

#### Notes

-   peptides : Predicts AMPs directly from peptide sequences
-   --fasta : Amino acid FASTA file (output of smORFinder)
-   --output : Output directory
-   --threads : Number of threads

### MetaPepticon

Identifies anticancer peptides (ACPs).

#### Setup

``` bash
conda create -c conda-forge -c bioconda -n snakemake snakemake=8.25.5 python=3.11
git clone --recursive https://github.com/arikanlab/MetaPepticon.git
```

#### Usage

``` bash
bash workflow/scripts/prepare_config.sh
conda activate snakemake
snakemake -s workflow/Snakefile --resources slot=1 --cores 16 --use-conda
```

#### Notes

-   --resources slot=1 : Allows only one job that requires the `slot`
    resource (assembly and toxicity prediction steps) to run at a time
-   --cores : Number of CPU cores
-   --use-conda : Runs each rule in its own conda environment

------------------------------------------------------------------------

## Part 7 — Constructing the refined plasmid set

The refined plasmid set is built from the Plasmer-predicted plasmid
contigs (`predPlasmid.fa`, Part 2) that are at least 2 kb long. The
contigs pass through five successive filters; each filter takes the
contigs kept by the previous one, and the contigs that remain after
Filter 5 form the refined plasmid set.

### Filter 1 — Plasmid-likeness scoring (PlasClass)

Scores each contig for plasmid-likeness from its k-mer composition
(k = 3–7), using the pretrained models shipped with PlasClass.
**Contigs with a plasmid score ≥ 0.5 are kept.**

#### Setup

``` bash
git clone https://github.com/Shamir-Lab/PlasClass.git
cd PlasClass

conda create -y -p ./classification-env python=3.8
conda activate ./classification-env

pip install numpy==1.22.0 scipy==1.10.0 joblib==1.2.0 "cython<3"
pip install "setuptools<60" wheel
pip install --no-build-isolation scikit-learn==0.21.3

python setup.py install
```

#### Usage

``` bash
python PlasClass/classify_fasta.py -f input.fa -o output.probs.out -p 8
```

#### Notes

-   -f : Input FASTA file (Plasmer-predicted plasmid contigs ≥ 2 kb)
-   -o : Output file (tab-separated contig name and plasmid score)
-   -p : Number of parallel processes

### Filter 2 — Gene-content classification (viralVerify)

Classifies contigs as plasmid, chromosome or virus from their gene
content (Prodigal gene prediction, `hmmsearch` against informative HMMs
and a Naive Bayes classifier). **Contigs classified as `Virus` or
`Chromosome` are removed; `Plasmid` and all `Uncertain` classes are
kept.**

#### Setup

``` bash
conda create -n viralverify -c conda-forge -c bioconda python=3.10 viralverify hmmer
conda activate viralverify
```

Install the database:

``` bash
wget -O nbc_hmms.hmm.gz https://ndownloader.figshare.com/files/17904323
gunzip nbc_hmms.hmm.gz
```

#### Usage

``` bash
viralverify -f input.fa -o outpath --hmm nbc_hmms.hmm -t 8 -p
```

#### Notes

-   -f : Input FASTA file (contigs kept in Filter 1)
-   -o : Output directory
-   --hmm : HMM database (decompressed)
-   -t : Number of threads
-   -p : Writes predicted plasmid contigs to a separate FASTA file

### Filter 3 — Second k-mer classification (PlasFlow)

Classifies contigs as plasmid or chromosome with a neural network based
on k-mer frequency signatures. **Contigs labelled `plasmid` are kept.**

#### Setup

``` bash
conda create -y -n plasflow -c conda-forge -c bioconda plasflow
conda activate plasflow
conda install -y -c conda-forge "libgfortran=3"
pip install --no-deps --force-reinstall protobuf==3.6.1
```

#### Usage

``` bash
PlasFlow.py --input input.fa --output output.tsv --threshold 0.7
```

#### Notes

-   --input : Input FASTA file (contigs kept in Filter 2)
-   --output : Output table. FASTA files of the classified sequences
    (`_plasmids`, `_chromosomes`, `_unclassified`) are written next to it.
-   --threshold : Minimum class probability for a label; contigs below
    it are `unclassified`

### Filter 4 — rRNA content removal (Infernal)

Detects rRNA genes with Rfam covariance models. **A hit is counted when
it covers ≥ 30% of the model's consensus length. Contigs with at least
one hit are removed.**

#### Setup

``` bash
conda create -y -n infernal -c conda-forge -c bioconda infernal
conda activate infernal
```

Install the Rfam rRNA covariance models:

| Rfam ID | Rfam family | rRNA |
|------------------------------------|------------------------------------|------------------------------------|
| RF00001 | 5S_rRNA | 5S rRNA |
| RF00177 | SSU_rRNA_bacteria | Bacterial 16S rRNA |
| RF02541 | LSU_rRNA_bacteria | Bacterial 23S rRNA |
| RF01959 | SSU_rRNA_archaea | Archaeal 16S rRNA |
| RF02540 | LSU_rRNA_archaea | Archaeal 23S rRNA |

``` bash
RFAM_IDS=(RF00001 RF00177 RF02541 RF01959 RF02540)

for id in "${RFAM_IDS[@]}"; do
    curl -sfL "https://rfam.org/family/${id}/cm"
done > rrna.cm

cmpress -F rrna.cm
```

#### Usage

``` bash
cmscan --noali \
       --cut_ga \
       -g \
       --nohmmonly \
       --rfam \
       --cpu 4 \
       --tblout output.tblout \
       rrna.cm \
       input.fa
```

#### Notes

-   input.fa : Input FASTA file (contigs kept in Filter 3)
-   --tblout : Tabular output file of the hits
-   --cpu : Number of threads
-   --cut_ga : Uses each Rfam family's gathering threshold as the
    reporting cutoff
-   -g : Glocal mode (the full model is aligned within the contig)
-   --nohmmonly : Scores with the covariance models, not only the
    profile-HMM stage
-   --rfam : Rfam-tuned preset for fast searches
-   --noali : Omits alignments from the output

### Filter 5 — Human-contamination removal (minimap2)

Aligns contigs to the human reference genome (GRCh38.p14, RefSeq
GCF_000001405.40). **Contigs with identity ≥ 95% and coverage ≥ 80% are
removed. The remaining contigs form the refined plasmid set.**

#### Setup

``` bash
conda create -y -n minimap2 -c conda-forge -c bioconda minimap2 bedtools samtools
conda activate minimap2
```

Download the reference genome:

``` bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
zcat GCF_000001405.40_GRCh38.p14_genomic.fna.gz > GRCh38.p14.fna
```

#### Usage

``` bash
minimap2 -x asm20 -c --secondary=no -I 1G -t 8 GRCh38.p14.fna input.fa > output.paf
```

#### Notes

-   input.fa : Input FASTA file (contigs kept in Filter 4)
-   -x asm20 : Preset for contig-to-genome alignment (up to ~5%
    divergence)
-   -c : Computes base-level alignments
-   --secondary=no : Reports only primary and supplementary alignments
-   -I 1G : Loads the index in 1 Gbp parts
-   -t : Number of threads

------------------------------------------------------------------------

# References

-   Q. Zhu, S. Gao, B. Xiao, Z. He, and S. Hu, “Plasmer: an Accurate and Sensitive Bacterial Plasmid Prediction Tool Based on Machine Learning of Shared k-mers and Genomic Features”, Microbiol. Spectr., vol. 11, no. 3, p. e0464522, Jun. 2023, doi: 10.1128/spectrum.04645-22. Github: [https://github.com/nekokoe/Plasmer](https://github.com/nekokoe/Plasmer)

-   R. Chikhi et al., “Logan: Planetary-Scale Genome Assembly Surveys Life’s Diversity”, bioRxiv, p. 2024.07.30.605881, Sep. 2025, doi: 10.1101/2024.07.30.605881. Github: [https://github.com/IndexThePlanet/Logan](https://github.com/IndexThePlanet/Logan)

-   W. Shen, B. Sipos, and L. Zhao, “SeqKit2: A Swiss army knife for sequence and alignment processing”, iMeta, vol. 3, no. 3, p. e191, Jun. 2024, doi: 10.1002/imt2.191. Github: [https://github.com/shenwei356/seqkit](https://github.com/shenwei356/seqkit)

-   M. Steinegger and J. Söding, “MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets”, Nat. Biotechnol., vol. 35, no. 11, pp. 1026-1028, Nov. 2017, doi: 10.1038/nbt.3988. Github: [https://github.com/soedinglab/MMseqs2](https://github.com/soedinglab/MMseqs2)

-   D. Hyatt, G.-L. Chen, P. F. LoCascio, M. L. Land, F. W. Larimer, and L. J. Hauser, “Prodigal: prokaryotic gene recognition and translation initiation site identification”, BMC Bioinformatics, vol. 11, p. 119, Mar. 2010, doi: 10.1186/1471-2105-11-119. Github: [https://github.com/hyattpd/prodigal](https://github.com/hyattpd/prodigal)

-   M. G. Durrant and A. S. Bhatt, “Automated Prediction and Annotation of Small Open Reading Frames in Microbial Genomes”, Cell Host Microbe, vol. 29, no. 1, pp. 121-131.e4, Jan. 2021, doi: 10.1016/j.chom.2020.11.002. Github: [https://github.com/bhattlab/SmORFinder](https://github.com/bhattlab/SmORFinder)

-   C. P. Cantalapiedra, A. Hernández-Plaza, I. Letunic, P. Bork, and J. Huerta-Cepas, “eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale”, Mol. Biol. Evol., vol. 38, no. 12, pp. 5825-5829, Dec. 2021, doi: 10.1093/molbev/msab293. Github: [https://github.com/eggnogdb/eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)

-   M. Feldgarden et al., “AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence”, Sci. Rep., vol. 11, no. 1, p. 12728, Jun. 2021, doi: 10.1038/s41598-021-91456-0. Github: [https://github.com/ncbi/amr](https://github.com/ncbi/amr)

-   B. P. Alcock et al., “CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database”, Nucleic Acids Res., vol. 51, no. D1, pp. D690-D699, Jan. 2023, doi: 10.1093/nar/gkac920. Github: [https://github.com/arpcard/rgi](https://github.com/arpcard/rgi)

-   J. Zheng, Q. Ge, Y. Yan, X. Zhang, L. Huang, and Y. Yin, “dbCAN3: automated carbohydrate-active enzyme and substrate annotation”, Nucleic Acids Res., vol. 51, no. W1, pp. W115-W121, Jul. 2023, doi: 10.1093/nar/gkad328. Github: [https://github.com/bcb-unl/run\_dbcan](https://github.com/bcb-unl/run_dbcan)

-   J. Russel, R. Pinilla-Redondo, D. Mayo-Muñoz, S. A. Shah, and S. J. Sørensen, “CRISPRCasTyper: Automated Identification, Annotation, and Classification of CRISPR-Cas Loci”, CRISPR J., vol. 3, no. 6, pp. 462-469, Dec. 2020, doi: 10.1089/crispr.2020.0059. Github: [https://github.com/Russel88/CRISPRCasTyper](https://github.com/Russel88/CRISPRCasTyper)

-   Santos-Júnior CD, Pan S, Zhao XM, Coelho LP. "Macrel: antimicrobial peptide screening in genomes and metagenomes." PeerJ. 2020 Dec 18;8:e10555. doi: 10.7717/peerj.10555. Github: [https://github.com/BigDataBiology/macrel](https://github.com/BigDataBiology/macrel)

-   Erözden AA, Tavşanlı N, Demirel G, Sanli NO, Çalışkan M, Arıkan M. "MetaPepticon: automated prediction of anticancer peptides from microbial genomes and metagenomes." PeerJ. 2026 Mar 27;14:e20990. doi: 10.7717/peerj.20990. Github: [https://github.com/arikanlab/MetaPepticon](https://github.com/arikanlab/MetaPepticon)

-   D. Pellow, I. Mizrahi, and R. Shamir, “PlasClass improves plasmid sequence classification”, PLoS Comput. Biol., vol. 16, no. 4, p. e1007781, Apr. 2020, doi: 10.1371/journal.pcbi.1007781. Github: [https://github.com/Shamir-Lab/PlasClass](https://github.com/Shamir-Lab/PlasClass)

-   D. Antipov, M. Raiko, A. Lapidus, and P. A. Pevzner, “Metaviral SPAdes: assembly of viruses from metagenomic data”, Bioinformatics, vol. 36, no. 14, pp. 4126-4129, Jul. 2020, doi: 10.1093/bioinformatics/btaa490. Github: [https://github.com/ablab/viralVerify](https://github.com/ablab/viralVerify)

-   P. S. Krawczyk, L. Lipinski, and A. Dziembowski, “PlasFlow: predicting plasmid sequences in metagenomic data using genome signatures”, Nucleic Acids Res., vol. 46, no. 6, p. e35, Apr. 2018, doi: 10.1093/nar/gkx1321. Github: [https://github.com/smaegol/PlasFlow](https://github.com/smaegol/PlasFlow)

-   E. P. Nawrocki and S. R. Eddy, “Infernal 1.1: 100-fold faster RNA homology searches”, Bioinformatics, vol. 29, no. 22, pp. 2933-2935, Nov. 2013, doi: 10.1093/bioinformatics/btt509. Github: [https://github.com/EddyRivasLab/infernal](https://github.com/EddyRivasLab/infernal)

-   N. Ontiveros-Palacios et al., “Rfam 15: RNA families database in 2025”, Nucleic Acids Res., vol. 53, no. D1, pp. D258-D267, Jan. 2025, doi: 10.1093/nar/gkae1023. Website: [https://rfam.org](https://rfam.org)

-   H. Li, “Minimap2: pairwise alignment for nucleotide sequences”, Bioinformatics, vol. 34, no. 18, pp. 3094-3100, Sep. 2018, doi: 10.1093/bioinformatics/bty191. Github: [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2)

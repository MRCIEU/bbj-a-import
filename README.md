# Import Biobank Japan GWAS

Data downloaded from here: http://jenger.riken.jp/en/result

The `bbj-a-metadata.csv` is created manually from the HTML tables on that page.

## Define data location

Create a config file called `config.json` which looks like:

```
{
  "datadir": "/path/to/data/dir"
}
```


## Download files

This runs GNU parallel download

```
bash dl.sh
```

## Select phenotypes for analysis

There are 137 QTL and 101 Case-control datasets. Some of these need to be combined because e.g. the autosomes and the x chromosome markers are stored in separate files and need to be combined into a single ID. The 

- Combine autosome and x chromosome files
- (sex stratified) files are tar.gz. They have a separate file per chromosome within the tarball, and have different column names. Perhaps just ignore these files?

Note that some of the files are `tar.gz` files, some are `gz`, some tarballs have pairs of gz files, one has a csv file. This script tries to organise all of it.

```
Rscript organise_phenotypes.r
```

This creates the `/data/dir/ready/input.csv` as well as organising the files into `ready/`

## Run pipeline

At this point we have

1. All the files downloaded and formatted in `/data/dir/ready`
2. A file called `/data/dir/ready/input.csv` which describes the data and specifies the metadata

We can now run the pipeline. Set it up:

```
module add languages/anaconda3/5.2.0-tflow-1.11
git clone --recurse-submodules git@github.com:MRCIEU/igd-hpc-pipeline.git
cd igd-hpc-pipeline/resources/gwas2vcf
python3 -m venv venv
source ./venv/bin/activate
./venv/bin/pip install -r requirements.txt
cd ../..
```



Some manual steps

```
datadir=$(jq -r .datadir ../config.json)

Rscript resources/metadata_to_json.r ${datadir}/ready/input.csv ${datadir}/ready ${datadir}/processed ${datadir}/ready/input_json.csv 8

Rscript resources/setup_directories.r ${datadir}/ready/input_json.csv 8

gwasdir="$(jq -r .datadir ../config.json)/processed"
echo `realpath ${gwasdir}` > gwasdir.txt
p=`pwd`
cd ${gwasdir}
ls --color=none -d * > ${p}/idlist.txt
cd ${p}
head idlist.txt
nid=`cat idlist.txt | wc -l`
echo "${nid} datasets"
```


Now run:

```
module add languages/anaconda3/5.2.0-tflow-1.11
snakemake -prk \
-j 400 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```


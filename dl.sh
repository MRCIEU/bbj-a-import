#!/bin/bash

datadir=$(jq -r .datadir config.json)

echo $datadir
echo $(ls $datadir)

p=$(pwd)
mkdir -p ${datadir}/dl
cd ${datadir}/dl

parallel --jobs 10 curl "http://jenger.riken.jp/{}/ > cc{}.txt.gz" ::: {1..101}
parallel --jobs 10 curl "curl http://jenger.riken.jp/{}analysisresult_qtl_download/ > qtl{}.txt.gz" ::: {1..137}


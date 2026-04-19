#!/usr/bin/env bash
set -euo pipefail

echo "######### source_dir/data/ #########"
path=$(fdog.setup -d ./ --getSourcepath)
ls "$path/data/"

WORKDIR=$(pwd)
DT_DIR="$WORKDIR/dt"

echo "######### TEST FDOG FUNCTIONS #########"

echo "TEST fdog.setup"
fdog.setup -d "$DT_DIR" --woFAS

echo "TEST fdog.checkData"
fdog.checkData -s "$DT_DIR/searchTaxa_dir" \
               -c "$DT_DIR/coreTaxa_dir" \
               -a "$DT_DIR/annotation_dir" \
               --reblast --ignoreAnno

echo "TEST fdog.showTaxa"
fdog.showTaxa

echo "TEST fdog.run"
fdog.run --seqFile infile.fa --jobName test \
         --refspec HUMAN@9606@qfo24_02 \
         --fasOff --group mammalia

echo "TEST fdog.assembly"
fdog.assembly --gene test \
              --refSpec HUMAN@9606@qfo24_02 \
              --augustus \
              --augustusRefSpec human \
              --coregroupPath core_orthologs/ \
              --out test_assembly \
              --fasoff

echo "Prepare seeds"
mkdir -p seeds
for i in 1 2 3; do
    cp "$path/data/infile.fa" "seeds/$i.fa"
done

echo "TEST fdogs.run"
fdogs.run --seqFolder seeds \
          --jobName test_multi \
          --refspec HUMAN@9606@qfo24_02 \
          --fasOff \
          --searchTaxa PARTE@5888@qfo24_02,THAPS@35128@qfo24_02 \
          --hmmScoreType sequence

echo "TEST fdog.addTaxon"
head "$DT_DIR/searchTaxa_dir/HUMAN@9606@qfo24_02/HUMAN@9606@qfo24_02.fa" > hm.fa
fdog.addTaxon -f hm.fa -i 9606 -o ./ -c -a

ls

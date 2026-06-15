#!/usr/bin/env bash
set -euo pipefail

assert_eq() {
    local expected="$1"
    local actual="$2"
    local msg="$3"

    if [ "$expected" != "$actual" ]; then
        echo "ERROR: $msg (expected=$expected, actual=$actual)"
        exit 1
    fi
}

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

lines=$(wc -l < test.phyloprofile)
assert_eq 10 "$lines" "test.phyloprofile line count"

echo "TEST fdog.assembly with miniprot"
fdog.assembly --gene test \
              --refSpec HUMAN@9606@qfo24_02 \
              --augustus \
              --augustusRefSpec human \
              --coregroupPath core_orthologs/ \
              --out test_assembly \
              --fasOff

lines=$(wc -l < test_assembly/test/test_og.fa)
assert_eq 4 "$lines" "test_assembly/test/test_og.fa line count"

echo "TEST fdog.assembly with blast"
fdog.assembly --gene test \
              --refSpec HUMAN@9606@qfo24_02 \
              --augustus \
              --augustusRefSpec human \
              --coregroupPath core_orthologs/ \
              --out test_assembly_blast \
              --fasOff \
              --searchTool blast

lines=$(wc -l < test_assembly_blast/test/test_og.fa)
assert_eq 4 "$lines" "test_assembly_blast/test/test_og.fa line count"

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

lines=$(wc -l < test_multi.phyloprofile)
assert_eq 13 "$lines" "test_multi.phyloprofile line count"

echo "TEST fdog.addTaxon"
head "$DT_DIR/searchTaxa_dir/HUMAN@9606@qfo24_02/HUMAN@9606@qfo24_02.fa" > hm.fa
fdog.addTaxon -f hm.fa -i 9606 -v testAddTaxon -o ./ -c -a

lines=$(wc -l < searchTaxa_dir/HOMSA\@9606\@testAddTaxon/HOMSA\@9606\@testAddTaxon.fa)
assert_eq 10 "$lines" "searchTaxa_dir/HOMSA\@9606\@testAddTaxon/HOMSA\@9606\@testAddTaxon.fa line count"

ls

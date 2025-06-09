# gwas-summary-stats

## Installation

Requires the [Rust programming language](https://www.rust-lang.org/).

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

```bash
cargo install --git https://github.com/GMELab/gwas-summary-stats
```

## Usage

You can use the wrapper James made at `/mnt/nfs/rigenenfs/shared_resources/GWAS_SUMMARY_STATS/GWAS_formatter_Rust_version_wrapper.sh`. It accepts three arguments, the trait name, the google sheets URL, and the output file.

If you would rather do it manually then you can do the following:

```bash
# samtools is annoying now with the new servers, thank you James for figuring this out 
server=$(hostname)
if [[ $server == "rigenep9" || $server == "rigenep10" ]]; then
    samtools="/mnt/nfs/rigenenfs/shared_resources/softwares/samtools-1.22/samtools-1.22_regenep9_10/samtools"
elif [[ $server == "rigenep2"]]; then
    samtools="/mnt/nfs/rigenenfs/shared_resources/softwares/samtools-1.22/samtools-1.22_regenep2/samtools"
elif [[ $server == "rigenep1"]]; then
    samtools="/mnt/nfs/rigenenfs/shared_resources/softwares/samtools-1.22/samtools-1.22_regenep1/samtools"
elif [[ $server == "ristatp15" || $server == "ristatp17" || $server == "ristatp19" || $server == "ristatp21" || $server == "ristatp23" || $server == "ristatp25" ]]; then
    samtools="/mnt/nfs/rigenenfs/shared_resources/softwares/samtools-1.22/samtools-1.22_ristatps/samtools"
fi

# raise the number of files we can open, samtools might error if we don't
ulimit -Sn $(ulimit -Hn)

# you'll likely only want to change the output file and the trait name
# all outputs are placed in the current working directory
output_file="final.txt.gz"
trait_name="example_trait"
gwas-summary-stats \
    --google-sheets-id=1KFuKkrAPYjWs79P4TVlDHg3jP98M1cd_R-aFhyQVBpg \
    --raw-input-dir=/mnt/nfs/rigenenfs/shared_resources/GWAS_SUMMARY_STATS/Raw \
    --liftover=/mnt/nfs/rigenenfs/shared_resources/softwares/LIFTOVER/liftOver \
    --grs-dir=/mnt/nfs/rigenenfs/shared_resources/GWAS_SUMMARY_STATS/formatted \
    --dbsnp-file=/mnt/nfs/rigenenfs/shared_resources/reference_files/DBSNP/DBSNP_155/dbsnp155.hg19_hg38.gmel_variants.info.slim.gnomAD_AF.txt.gz \
    --samtools=${samtools} \
    --fasta-ref=/mnt/nfs/rigenenfs/shared_resources/softwares/Homo_sapiens_assembly38.fasta \
    --liftover-dir=/mnt/nfs/rigenenfs/shared_resources/softwares/LIFTOVER/ \
    --samtools-chunk-size=50000 \
    --output-file=${output_file} \
    --trait-name=${trait_name}
```

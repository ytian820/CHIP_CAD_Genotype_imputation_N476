version 1.0

# Screen one single-sample GRCh37 VCF.gz; keep chromosome named X (not chrX or 23) and biallelic A/C/G/T SNPs.
# Exclude DUPE/FAIL_REF/TRIALLELIC/ZEROED_OUT_ASSAY and all records at duplicated candidate coordinates.
# Preserve sample ID, record order, allele order and GT; remove INFO and other FORMAT fields.
# Input: Terra-localized VCF.gz; no input index needed.
# Outputs: GT-only chrX VCF.gz and tabix index.
# Keep PAR and non-PAR together; no sex inference or ploidy recoding.

# Part 1. Workflow inputs and outputs
workflow ScreenChrX {
    input {
        File vcf
    }

    call ScreenChrXVcf {
        input: vcf = vcf
    }

    output {
        File screened_vcf = ScreenChrXVcf.screened_vcf
        File screened_vcf_index = ScreenChrXVcf.screened_vcf_index
    }
}

# Part 2. Screen one localized VCF
task ScreenChrXVcf {
    input {
        File vcf
    }

    command <<<
        set -euo pipefail

        bcftools --version

        ## Steps 1-4. ChrX -> FILTER -> biallelic SNPs -> GT only.
        bcftools view --no-version -t X -Ou "~{vcf}" \
            | bcftools view --no-version \
                -e 'FILTER~"DUPE" || FILTER~"FAIL_REF" || FILTER~"TRIALLELIC" || FILTER~"ZEROED_OUT_ASSAY"' -Ou \
            | bcftools view --no-version -m2 -M2 -v snps \
                -i 'REF~"^[ACGTacgt]$" && ALT~"^[ACGTacgt]$"' -Ov \
            | awk '/^#/ || toupper($4) != toupper($5)' \
            | bcftools annotate --no-version -x INFO,^FORMAT/GT -Ov -o candidates.vcf

        ## Step 5. Exclude every record at duplicated candidate coordinates.
        ### 1. List every coordinate (CHROM, POS) that occurs more than once.
        awk -F '\t' '!/^#/ { print $1 "\t" $2 }' candidates.vcf \
            | sort | uniq -d > duplicate_sites.tsv

        ### 2. Remove the records at those coordinates (-T ^FILE = exclude listed sites).
        if [ -s duplicate_sites.tsv ]; then
            bcftools view --no-version -T ^duplicate_sites.tsv \
                -Oz -o chrX.screened.vcf.gz candidates.vcf
        else
            bcftools view --no-version \
                -Oz -o chrX.screened.vcf.gz candidates.vcf
        fi
        echo "Duplicated coordinates removed: $(wc -l < duplicate_sites.tsv)"

        ## Step 6. Index and verify the output files.
        bcftools index --tbi chrX.screened.vcf.gz
        test -s chrX.screened.vcf.gz
        test -s chrX.screened.vcf.gz.tbi
        echo "Retained SNPs: $(bcftools index -n chrX.screened.vcf.gz)"
        echo "ChrX_screening: complete"
    >>>

    output {
        File screened_vcf = "chrX.screened.vcf.gz"
        File screened_vcf_index = "chrX.screened.vcf.gz.tbi"
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.0.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}

version 1.0

# Purpose: infer sample sex for manual review; no phenotype input or sample removal.
# Input: merged GRCh37 chrX VCF (476 samples). Outputs: PLINK2 .sexcheck table, SNP list and logs.
# Filters: biallelic A/C/G/T; geno 0.01; MAF 0.05; LD 500kb 1 0.2; no HWE filter.
# X F thresholds: max-female 0.2 (SNPSEX 2), min-male 0.8 (SNPSEX 1); otherwise uncertain.

workflow SexInference {
    input {
        File chrX_vcf
    }

    call InferSex {
        input: chrX_vcf = chrX_vcf
    }

    output {
        File sex_check = InferSex.sex_check
        File pruned_variants = InferSex.pruned_variants
        Array[File] plink_logs = InferSex.plink_logs
    }
}

task InferSex {
    input {
        File chrX_vcf
    }

    command <<<
        set -euo pipefail

        plink2 --version
        echo "GRCh37 non-PAR X; retain all samples; geno 0.01; MAF 0.05; no HWE filter"
        echo "LD pruning: 500kb 1 0.2; max-female-xf=0.2; min-male-xf=0.8"

        # Part 1. Import chrX and split PAR
        plink2 --vcf "~{chrX_vcf}" --double-id --vcf-half-call missing \
            --lax-chrx-import --split-par b37 \
            --snps-only just-acgt --min-alleles 2 --max-alleles 2 \
            --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all \
            --make-pgen --threads 1 --memory 3000 --out chrX

        # Part 2. Filter and prune non-PAR SNPs, then infer sex
        # After splitting PAR, --chr X selects only non-PAR variants.
        plink2 --pfile chrX --chr X --geno 0.01 --maf 0.05 \
            --indep-pairwise 500kb 1 0.2 \
            --threads 1 --memory 3000 --out pruning

        plink2 --pfile chrX --chr X --extract pruning.prune.in \
            --check-sex max-female-xf=0.2 min-male-xf=0.8 cols=fid,xf \
            --threads 1 --memory 3000 --out inferred
    >>>

    output {
        File sex_check = "inferred.sexcheck"
        File pruned_variants = "pruning.prune.in"
        Array[File] plink_logs = ["chrX.log", "pruning.log", "inferred.log"]
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.1.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}

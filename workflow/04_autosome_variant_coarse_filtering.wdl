version 1.0

# Remove autosomal variants with missingness > 5%; retain all samples. No phenotype.
# Input: sorted, GT-only autosomal cohort VCF from 03. Outputs: PGEN/PVAR/PSAM and PLINK log.
# Treat half-missing GTs as missing; preserve sample IDs as both FID and IID.
# No MAF, HWE, LD or additional chromosome filter.

# Part 1. Workflow interface
workflow FilterAutosomes {
    input {
        File vcf
    }

    call FilterMissingness {
        input: vcf = vcf
    }

    output {
        File filtered_pgen = FilterMissingness.filtered_pgen
        File filtered_pvar = FilterMissingness.filtered_pvar
        File filtered_psam = FilterMissingness.filtered_psam
        File plink_log = FilterMissingness.plink_log
    }
}

# Part 2. Filter with PLINK2 and save the retained PGEN dataset
task FilterMissingness {
    input {
        File vcf
    }

    command <<<
        set -euo pipefail
        echo "Variant missingness threshold: 0.05; retain all samples"
        plink2 --version

        # Import and filter variants only; no --mind or sample exclusions.
        plink2 --vcf "~{vcf}" \
            --double-id --vcf-half-call missing \
            --geno 0.05 \
            --make-pgen \
            --threads 1 --memory 3000 \
            --out N476cohort.autosomes.missingness05
    >>>

    output {
        File filtered_pgen = "N476cohort.autosomes.missingness05.pgen"
        File filtered_pvar = "N476cohort.autosomes.missingness05.pvar"
        File filtered_psam = "N476cohort.autosomes.missingness05.psam"
        File plink_log = "N476cohort.autosomes.missingness05.log"
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.0.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}

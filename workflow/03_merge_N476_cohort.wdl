version 1.0

# Merge one GRCh37 sample_set; no phenotype. Inputs: 01/02 VCFs only.
# Outputs: separate multi-sample autosomal and X VCF.gz files with tabix indexes,
# plus the input VCF list of each merge (localized paths; line order = sample order).
# Require sorted inputs, unique sample names and compatible alleles.
# Tabix indexes are rebuilt here from each VCF, so every index matches its own VCF.
# Keep the union of records; absent genotypes remain missing; no new multiallelic records.
# No QC or alignment here; preserve X genotypes and keep PAR/non-PAR together.

# Part 1. Cohort inputs and outputs
workflow MergeCohort {
    input {
        Array[File]+ autosomal_vcfs
        Array[File]+ chrX_vcfs
    }

    # Independent calls: neither consumes the other's outputs.
    call MergeVcfs as MergeAutosomes {
        input:
            vcfs = autosomal_vcfs,
            output_prefix = "N476cohort.autosomes"
    }

    call MergeVcfs as MergeChrX {
        input:
            vcfs = chrX_vcfs,
            output_prefix = "N476cohort.chrX"
    }

    output {
        File merged_autosomal_vcf = MergeAutosomes.merged_vcf
        File merged_autosomal_vcf_index = MergeAutosomes.merged_vcf_index
        File autosomal_input_vcf_list = MergeAutosomes.input_vcf_list
        File merged_chrX_vcf = MergeChrX.merged_vcf
        File merged_chrX_vcf_index = MergeChrX.merged_vcf_index
        File chrX_input_vcf_list = MergeChrX.input_vcf_list
    }
}

# Part 2. Index the localized VCFs and merge
task MergeVcfs {
    input {
        Array[File]+ vcfs
        String output_prefix
    }

    command <<<
        set -euo pipefail

        echo "Input VCF count: ~{length(vcfs)}"
        bcftools --version

        ## Step 1. Build a tabix index beside every localized VCF.
        vcf_list="~{write_lines(vcfs)}"        # one localized VCF path per line
        while read -r vcf; do
            bcftools index --tbi --force "$vcf"
        done < "$vcf_list"
        cp "$vcf_list" "~{output_prefix}.input_vcfs.list"

        ## Step 2. Merge, index and verify the output files.
        bcftools merge --no-version --merge none --threads 1 \
            --file-list "$vcf_list" \
            --output-type z --output "~{output_prefix}.vcf.gz"
        bcftools index --tbi "~{output_prefix}.vcf.gz"
        test -s "~{output_prefix}.vcf.gz"
        test -s "~{output_prefix}.vcf.gz.tbi"
        echo "Cohort_merge: complete"
    >>>

    output {
        File merged_vcf = "~{output_prefix}.vcf.gz"
        File merged_vcf_index = "~{output_prefix}.vcf.gz.tbi"
        File input_vcf_list = "~{output_prefix}.input_vcfs.list"
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.0.0"
        cpu: 2
        memory: "16 GB"
        disks: "local-disk 100 SSD"
    }
}

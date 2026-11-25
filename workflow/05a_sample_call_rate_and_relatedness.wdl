version 1.0

# Report sample missingness and KING-robust kinship; no phenotype.
# Input: autosomal PGEN/PVAR/PSAM from 04, one cohort/sample_set.
# Outputs: SMISS, FID/IID list with F_MISS > 0.02, KING table of all pairs,
# KING table of pairs with kinship >= 0.0884, suggested related-sample exclusion list, PLINK logs.
# Call rate = 1 - F_MISS; the high-missingness and related-sample lists do not remove samples.
# Related-sample list: one member of each pair >= 0.0884 (2nd degree or closer), chosen by PLINK2
# to keep the largest unrelated set; for HWE and unrelated-sample analyses, not for imputation upload.
# Retain all samples and variants; no new PGEN dataset or MAF/HWE/LD filtering.

# Part 1. Workflow interface
workflow SampleCallRateAndKing {
    input {
        File pgen
        File pvar
        File psam
    }

    call ReportSampleQcAndKing {
        input: pgen = pgen, pvar = pvar, psam = psam
    }

    output {
        File sample_missingness = ReportSampleQcAndKing.sample_missingness
        File sample_qc_out_missing = ReportSampleQcAndKing.sample_qc_out_missing
        File king_table = ReportSampleQcAndKing.king_table
        File related_pairs = ReportSampleQcAndKing.related_pairs
        File related_samples_out = ReportSampleQcAndKing.related_samples_out
        Array[File] plink_logs = ReportSampleQcAndKing.plink_logs
    }
}

# Part 2. Sample missingness and relatedness reports
task ReportSampleQcAndKing {
    input {
        File pgen
        File pvar
        File psam
    }

    command <<<
        set -euo pipefail
        echo "Flag sample missingness > 0.02; report KING-robust; retain all samples and variants"
        plink2 --version

        ## Step 1. Report missingness for every input sample.
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --missing sample-only \
            --threads 1 --memory 3000 \
            --out N476cohort.autosomes.sample_qc

        ## Step 2. List samples with missingness > 2% without removing them.
        # .smiss columns: #FID  IID  MISSING_CT  OBS_CT  F_MISS
        awk -F '\t' 'BEGIN { OFS = "\t" }
            NR == 1 {
                if ($1 != "#FID" || $2 != "IID" || $5 != "F_MISS") {
                    print "ERROR: unexpected .smiss header" > "/dev/stderr"
                    exit 1
                }
                print "#FID", "IID"
                next
            }
            $5 > 0.02 { print $1, $2 }
        ' N476cohort.autosomes.sample_qc.smiss > 05a_sample_qc_out_missing.tsv

        ## Step 3. Calculate KING-robust for all input sample pairs.
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --make-king-table \
            --threads 1 --memory 3000 \
            --out N476cohort.autosomes.relatedness

        ## Step 4. List only the pairs with kinship >= 0.0884 (2nd degree or closer).
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --make-king-table --king-table-filter 0.0884 \
            --threads 1 --memory 3000 \
            --out N476cohort.autosomes.related_pairs

        ## Step 5. Suggest one sample to exclude per related pair, without removing anyone.
        # Writes .king.cutoff.in.id (unrelated set) and .king.cutoff.out.id (excluded samples).
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --king-cutoff 0.0884 \
            --threads 1 --memory 3000 \
            --out N476cohort.autosomes.related

        echo "Sample QC reports: complete"
    >>>

    output {
        File sample_missingness = "N476cohort.autosomes.sample_qc.smiss"
        File sample_qc_out_missing = "05a_sample_qc_out_missing.tsv"
        File king_table = "N476cohort.autosomes.relatedness.kin0"
        File related_pairs = "N476cohort.autosomes.related_pairs.kin0"
        File related_samples_out = "N476cohort.autosomes.related.king.cutoff.out.id"
        Array[File] plink_logs = [
            "N476cohort.autosomes.sample_qc.log",
            "N476cohort.autosomes.relatedness.log",
            "N476cohort.autosomes.related_pairs.log",
            "N476cohort.autosomes.related.log"
        ]
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.0.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}

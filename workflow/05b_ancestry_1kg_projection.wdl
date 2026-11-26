version 1.0

# Purpose: build reference PCA axes and project the study cohort onto them.

# Part 1. Workflow interface
workflow AncestryProjection {
    input {
        File study_pgen
        File study_pvar
        File study_psam
        File reference_1kg_pgen
        File reference_1kg_pvar
        File reference_1kg_psam
        File reference_1kg_deg2_related_samples
    }

    call ProjectOnto1kg {
        input:
            study_pgen = study_pgen, study_pvar = study_pvar, study_psam = study_psam,
            reference_1kg_pgen = reference_1kg_pgen, reference_1kg_pvar = reference_1kg_pvar,
            reference_1kg_psam = reference_1kg_psam,
            reference_1kg_deg2_related_samples = reference_1kg_deg2_related_samples
    }

    output {
        File study_pcs = ProjectOnto1kg.study_pcs
        File reference_pcs = ProjectOnto1kg.reference_pcs
        File pca_variants = ProjectOnto1kg.pca_variants
        File allele_weights = ProjectOnto1kg.allele_weights
        File allele_counts = ProjectOnto1kg.allele_counts
        File eigenvalues = ProjectOnto1kg.eigenvalues
        Array[File] plink_logs = ProjectOnto1kg.plink_logs
    }
}

# Part 2. Match SNPs, build the reference PCA and project both datasets
task ProjectOnto1kg {
    input {
        File study_pgen
        File study_pvar
        File study_psam
        File reference_1kg_pgen
        File reference_1kg_pvar
        File reference_1kg_psam
        File reference_1kg_deg2_related_samples
    }

    command <<<
        set -euo pipefail
        plink2 --version

        ## Step 1. Study: autosomal biallelic non-palindromic SNPs, renamed CHROM:POS:allele1:allele2
        ##         (alleles in sorted order), unique IDs only, missingness <= 0.01.
        plink2 --pgen "~{study_pgen}" --pvar "~{study_pvar}" --psam "~{study_psam}" \
            --autosome --snps-only just-acgt --min-alleles 2 --max-alleles 2 --exclude-palindromic-snps \
            --set-all-var-ids '@:#:$1:$2' --rm-dup exclude-all \
            --geno 0.01 --make-pgen --write-snplist \
            --threads 2 --memory 12000 --out study.qc

        ## Step 2. Reference: drop official relatives, apply the same SNP rules, keep SNPs shared
        ##         with the study, then missingness <= 0.01, MAF >= 0.05, outside long-range LD.
        # Long-range high-LD regions, Build 37, 24 regions; treated as 1-based inclusive.
        # Based on Price et al. 2008 and Anderson et al. 2010:
        # https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
        # Write GRCh37 long-range LD regions: chromosome, start, end, region ID.
        cat > long_ld.txt <<'EOF'
        1 48000000 52000000 1
        2 86000000 100500000 2
        2 134500000 138000000 3
        2 183000000 190000000 4
        3 47500000 50000000 5
        3 83500000 87000000 6
        3 89000000 97500000 7
        5 44500000 50500000 8
        5 98000000 100500000 9
        5 129000000 132000000 10
        5 135500000 138500000 11
        6 25000000 35000000 12
        6 57000000 64000000 13
        6 140000000 142500000 14
        7 55000000 66000000 15
        8 7000000 13000000 16
        8 43000000 50000000 17
        8 112000000 115000000 18
        10 37000000 43000000 19
        11 46000000 57000000 20
        11 87500000 90500000 21
        12 33000000 40000000 22
        12 109500000 112000000 23
        20 32000000 34500000 24
        EOF
        # --make-founders clears parental IDs so every retained reference sample counts
        # in the frequency, MAF and LD calculations.
        plink2 --pgen "~{reference_1kg_pgen}" --pvar "~{reference_1kg_pvar}" --psam "~{reference_1kg_psam}" \
            --remove "~{reference_1kg_deg2_related_samples}" --make-founders \
            --autosome --snps-only just-acgt --max-alleles 2 --exclude-palindromic-snps \
            --set-all-var-ids '@:#:$1:$2' --rm-dup exclude-all \
            --extract study.qc.snplist --geno 0.01 --maf 0.05 \
            --exclude range long_ld.txt --make-pgen --sort-vars \
            --threads 2 --memory 12000 --out reference.common

        ## Step 3. LD-prune the reference SNPs.
        plink2 --pfile reference.common --indep-pairwise 500kb 1 0.2 --indep-order 2 \
            --threads 2 --memory 12000 --out reference.ld

        ## Step 4. Reference PCA: 30 PCs, allele weights and allele counts.
        plink2 --pfile reference.common --extract reference.ld.prune.in --freq counts \
            --pca 30 allele-wts vcols=chrom,ref,alt \
            --threads 2 --memory 12000 --out reference.pca

        ## Step 5. Project reference and study samples with the same weights and frequencies.
        # Weight columns 6-25 are PC1-PC20: the first 20 of the 30 reference PCs are projected.
        # Missing calls are skipped (no mean imputation); use the per-sample score averages (PC*_AVG).
        # The score table also reports how many allele observations each sample contributed.
        plink2 --pfile reference.common --extract reference.ld.prune.in \
            --read-freq reference.pca.acount \
            --score reference.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-25 \
            --threads 2 --memory 12000 --out reference.projection
        plink2 --pfile study.qc --extract reference.ld.prune.in \
            --read-freq reference.pca.acount \
            --score reference.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-25 \
            --threads 2 --memory 12000 --out study.projection
        echo "Reference PCA and study projection: complete"
    >>>

    output {
        File study_pcs = "study.projection.sscore"
        File reference_pcs = "reference.projection.sscore"
        File pca_variants = "reference.ld.prune.in"
        File allele_weights = "reference.pca.eigenvec.allele"
        File allele_counts = "reference.pca.acount"
        File eigenvalues = "reference.pca.eigenval"
        Array[File] plink_logs = glob("*.log")
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.0.0"
        cpu: 2
        memory: "16 GB"
        disks: "local-disk 80 SSD"
    }
}

version 1.0

# Purpose: flag autosomal heterozygosity outliers (PLINK2 F) among EUR-assigned study samples.

workflow Heterozygosity {
    input {
        File pgen
        File pvar
        File psam
        File ancestry_assignments
    }

    call FlagHeterozygosity {
        input:
            pgen = pgen, pvar = pvar, psam = psam,
            ancestry_assignments = ancestry_assignments
    }

    output {
        File heterozygosity = FlagHeterozygosity.heterozygosity
        File sample_qc_out_het = FlagHeterozygosity.sample_qc_out_het
        Array[File] plink_logs = FlagHeterozygosity.plink_logs
    }
}

task FlagHeterozygosity {
    input {
        File pgen
        File pvar
        File psam
        File ancestry_assignments
    }

    command <<<
        set -euo pipefail
        plink2 --version

        ## Step 1. Select common autosomal SNPs outside long-range LD, then LD-prune.
        # Long-range high-LD regions, Build 37, 24 regions; treated as 1-based inclusive.
        # Based on Price et al. 2008 and Anderson et al. 2010:
        # https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
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

        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --autosome --snps-only just-acgt --min-alleles 2 --max-alleles 2 \
            --geno 0.01 --maf 0.05 --exclude range long_ld.txt \
            --hwe 1e-6 keep-fewhet \
            --indep-pairwise 500kb 1 0.2 \
            --threads 1 --memory 3000 --out pruning

        ## Step 2. Calculate F on the pruned SNPs only.
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --extract pruning.prune.in --het \
            --threads 1 --memory 3000 --out cohort

        ## Step 3. Flag EUR samples with F outside the EUR mean +/- 3 SD.
        python - <<'PY'
        import pandas as pd

        # .het columns: #FID  IID  O(HOM)  E(HOM)  OBS_CT  F
        samples = pd.read_csv('cohort.het', sep='\t', dtype={'#FID': str, 'IID': str})
        ancestry = pd.read_csv('~{ancestry_assignments}', sep='\t', dtype={'#FID': str, 'IID': str},
                               usecols=['#FID', 'IID', 'predicted_group'])
        samples = samples.merge(ancestry, on=['#FID', 'IID'], how='left', validate='one_to_one')

        is_eur = samples['predicted_group'] == 'EUR'
        mean = samples.loc[is_eur, 'F'].mean()
        sd = samples.loc[is_eur, 'F'].std()
        print('EUR samples:', is_eur.sum(), ' mean F:', round(mean, 5), ' SD:', round(sd, 5))

        samples['het_flag'] = 'not_evaluated'
        samples.loc[is_eur, 'het_flag'] = 'pass'
        samples.loc[is_eur & (samples['F'] < mean - 3 * sd), 'het_flag'] = 'low_F'
        samples.loc[is_eur & (samples['F'] > mean + 3 * sd), 'het_flag'] = 'high_F'
        print('Samples per flag:', samples['het_flag'].value_counts().to_dict())

        samples.to_csv('heterozygosity.tsv', sep='\t', index=False)
        flagged = samples['het_flag'].isin(['low_F', 'high_F'])
        samples.loc[flagged, ['#FID', 'IID']].to_csv('05d_sample_qc_out_het.tsv', sep='\t', index=False)
        PY
    >>>

    output {
        File heterozygosity = "heterozygosity.tsv"
        File sample_qc_out_het = "05d_sample_qc_out_het.tsv"
        Array[File] plink_logs = glob("*.log")
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.1.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}

version 1.0

# Prepare the public 1KG Phase 3 GRCh37 reference for 05b; no phenotype or inputs.
# Full variant set, no INFO annotations, KING-corrected PSAM; no filtering or PCA.
# Outputs: uncompressed PGEN/PVAR/PSAM, official population panel, PLINK2 list of
# second-degree-or-closer relatives to exclude (KING cutoff 0.0884).
# Requires HTTPS access; sources: https://www.cog-genomics.org/plink/2.0/resources

# Part 1. Workflow interface
workflow Prepare1kgReference {
    call Download1kgReference

    output {
        File reference_pgen = Download1kgReference.reference_pgen
        File reference_pvar = Download1kgReference.reference_pvar
        File reference_psam = Download1kgReference.reference_psam
        File population_panel = Download1kgReference.population_panel
        File related_samples_out = Download1kgReference.related_samples_out
    }
}

# Part 2. Download and decompress the reference bundle
task Download1kgReference {
    command <<<
        set -euo pipefail
        ## Step 1. Download the reference files.
        curl -fL --retry 4 --retry-all-errors --retry-delay 10 \
            -o "reference.pgen.zst" \
            "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"

        curl -fL --retry 4 --retry-all-errors --retry-delay 10 \
            -o "reference.pvar.zst" \
            "https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1"

        curl -fL --retry 4 --retry-all-errors --retry-delay 10 \
            -o "reference.psam" \
            "https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1"

        curl -fL --retry 4 --retry-all-errors --retry-delay 10 \
            -o "integrated_call_samples_v3.20130502.ALL.panel" \
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

        curl -fL --retry 4 --retry-all-errors --retry-delay 10 \
            -o "deg2_phase3.king.cutoff.out.id" \
            "https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1"

        ## Step 2. Decompress using PLINK2; all three files share one prefix.
        plink2 --zst-decompress reference.pgen.zst reference.pgen
        plink2 --zst-decompress reference.pvar.zst reference.pvar
    >>>

    output {
        File reference_pgen = "reference.pgen"
        File reference_pvar = "reference.pvar"
        File reference_psam = "reference.psam"
        File population_panel = "integrated_call_samples_v3.20130502.ALL.panel"
        File related_samples_out = "deg2_phase3.king.cutoff.out.id"
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.0.0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 40 SSD"
    }
}

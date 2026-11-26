version 1.0

# Purpose: assign 1KG super-population labels to study samples from the 05b projected PCs.

workflow AncestryAssignment {
    input {
        File reference_pcs
        File study_pcs
    }

    call AssignAncestry {
        input: reference_pcs = reference_pcs, study_pcs = study_pcs
    }

    output {
        File ancestry_assignments = AssignAncestry.ancestry_assignments
    }
}

task AssignAncestry {
    input {
        File reference_pcs
        File study_pcs
    }

    command <<<
        set -euo pipefail
        python - <<'PY'
        import pandas as pd
        from sklearn.ensemble import RandomForestClassifier

        pc_columns = [f'PC{k}_AVG' for k in range(1, 11)]
        min_probability = 0.80
        print('Random forest: PC1-PC10, 500 trees, seed 42, assignment threshold 0.80')

        ## Step 1. Read the PC tables; the reference table already carries the SuperPop label.
        reference = pd.read_csv('~{reference_pcs}', sep='\t', usecols=['SuperPop'] + pc_columns)
        study = pd.read_csv('~{study_pcs}', sep='\t', dtype={'#FID': str, 'IID': str}, usecols=['#FID', 'IID'] + pc_columns)
        print('Reference samples per group:', reference['SuperPop'].value_counts().to_dict())

        ## Step 2. Train on all reference samples and predict the study samples.
        model = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=2)
        model.fit(reference[pc_columns], reference['SuperPop'])
        probabilities = pd.DataFrame(model.predict_proba(study[pc_columns]), columns=model.classes_)

        ## Step 3. Label each study sample; below the threshold it is unassigned.
        result = study[['#FID', 'IID']].copy()
        result['predicted_group'] = probabilities.idxmax(axis=1)
        result['max_probability'] = probabilities.max(axis=1)
        result['assigned_group'] = result['predicted_group'].where(
            result['max_probability'] >= min_probability, 'unassigned')
        result = pd.concat([result, probabilities.add_prefix('prob_')], axis=1)
        result.to_csv('ancestry_assignments.tsv', sep='\t', index=False)
        print('Study samples per group:', result['assigned_group'].value_counts().to_dict())
        PY
    >>>

    output {
        File ancestry_assignments = "ancestry_assignments.tsv"
    }

    runtime {
        docker: "ytian820/chip_cad_genotype:1.1.0"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 10 SSD"
    }
}

SENSITIVITY_INPUTS = f"{WORK}/sample_sensitivity/inputs"
SENSITIVITY_REFITS = f"{RESULTS}/fusarium/sample_sensitivity_models"
SENSITIVITY_RESULTS = f"{RESULTS}/fusarium/sample_sensitivity"


rule sample_sensitivity:
    input:
        SENSITIVITY_RESULTS,


rule complete_study:
    input:
        rules.all.input,
        MACROPHOMINA_RESULTS,
        SENSITIVITY_RESULTS,


rule prepare_sample_sensitivity:
    input:
        data=rules.prepare_expression_inputs.output,
        mapping=rules.summarize_mapping_3dai.output,
        settings="resources/sample_sensitivity.json",
        script=f"{SCRIPTS}/sample_sensitivity.py",
        helper=f"{SCRIPTS}/score_sensitivity.py",
    output:
        directory(SENSITIVITY_INPUTS)
    shell:
        """
        python {input.script:q} prepare \
          --counts {input.data}/counts_primary.tsv \
          --samples {input.data}/samples_primary.tsv \
          --gene-sets {input.data}/functional_gene_sets.tsv \
          --mapping {input.mapping}/sample_mapping.tsv \
          --settings {input.settings:q} --output {output:q}
        """


rule refit_sample_sensitivity:
    input:
        prepared=rules.prepare_sample_sensitivity.output,
        script=f"{SCRIPTS}/refit_sample_sensitivity.R",
    output:
        directory(SENSITIVITY_REFITS)
    shell:
        """
        Rscript {input.script:q} {input.prepared:q} {output:q} {input.prepared}/settings.json
        """


rule summarise_sample_sensitivity:
    input:
        prepared=rules.prepare_sample_sensitivity.output,
        models=rules.refit_sample_sensitivity.output,
        script=f"{SCRIPTS}/sample_sensitivity.py",
        helper=f"{SCRIPTS}/score_sensitivity.py",
    output:
        directory(SENSITIVITY_RESULTS)
    shell:
        """
        python {input.script:q} summarise --prepared {input.prepared:q} \
          --refit {input.models:q} --output {output:q}
        """

rule retrieve_rnaseq:
	input:
		formatted_clinical_summary=formatted_summary
	output:
		formatted_raw_rna_counts=formatted_rnaseq
	conda: '../envs/environment.yaml'
	resources:
		mem_mb=100000
	script:
		'../scripts/retrieve_raw_rnaseq.R'

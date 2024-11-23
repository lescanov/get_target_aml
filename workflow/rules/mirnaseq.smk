rule retrieve_mirnaseq:
	input:
		formatted_clinical_summary=formatted_summary
	output:
		formatted_raw_mirna_counts=formatted_mirnaseq
	conda: '../envs/environment.yaml'
	script:
		'../scripts/retrieve_raw_mirnaseq.R'

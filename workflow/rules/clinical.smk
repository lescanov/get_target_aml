rule format_clinical:
	input:
		clinical_summary=raw_clinical
	output:
		formatted_clinical_summary=formatted_summary
	conda: '../envs/environment.yaml'
	script:
		'../scripts/wrangle_clinical_summary.R'

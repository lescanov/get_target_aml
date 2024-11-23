raw_clinical = config['raw_clinical']
formatted_summary = 'results/clinical/formatted_target_clinical_summary.csv'
formatted_rnaseq='results/rnaseq/formatted_raw_target_rna_counts.csv'
formatted_mirnaseq='results/mirnaseq/formatted_raw_target_mirna_counts.csv'

def get_final_output():
    final_output = [formatted_summary, formatted_rnaseq, formatted_mirnaseq]
    return final_output


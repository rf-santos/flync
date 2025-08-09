import pandas as pd
import sys

def classify_transcripts(cuffcompare_gtf, cpat_orf, cpat_no_orf, cutoff, output_non_coding, output_coding):
    """
    Classifies transcripts as coding or non-coding based on cuffcompare class code and CPAT coding probability.
    """
    # Read cuffcompare gtf
    cuff_gtf = pd.read_csv(cuffcompare_gtf, sep='\t', header=None, comment='#')
    cuff_gtf['class_code'] = cuff_gtf[8].str.extract(r'class_code "(.)";')
    cuff_gtf['transcript_id'] = cuff_gtf[8].str.extract(r'transcript_id "(.*?)";')

    # Read cpat results
    cpat_orf_df = pd.read_csv(cpat_orf, sep='\t')
    cpat_no_orf_df = pd.read_csv(cpat_no_orf, header=None, names=['transcript_id'])

    with open(cutoff) as f:
        cpat_cutoff = float(f.read().strip())

    # Get coding and non-coding transcript ids from cpat
    coding_ids = set(cpat_orf_df[cpat_orf_df['Coding_prob'] > cpat_cutoff]['ID'])
    non_coding_ids = set(cpat_orf_df[cpat_orf_df['Coding_prob'] <= cpat_cutoff]['ID'])
    non_coding_ids.update(set(cpat_no_orf_df['transcript_id']))

    # Classify transcripts
    linc_rna = cuff_gtf[cuff_gtf['class_code'] == 'u']
    intronic_rna = cuff_gtf[cuff_gtf['class_code'] == 'i']
    antisense_rna = cuff_gtf[cuff_gtf['class_code'] == 'x']
    isoforms = cuff_gtf[cuff_gtf['class_code'] == 'j']

    # Filter by coding potential
    new_linc_rna_non_coding = linc_rna[linc_rna['transcript_id'].isin(non_coding_ids)]
    new_intronic_rna_non_coding = intronic_rna[intronic_rna['transcript_id'].isin(non_coding_ids)]
    new_antisense_rna_non_coding = antisense_rna[antisense_rna['transcript_id'].isin(non_coding_ids)]

    new_isoforms_coding = isoforms[isoforms['transcript_id'].isin(coding_ids)]

    # Get micro-ORFs
    micro_orfs_ids = set(cpat_orf_df[(cpat_orf_df['Coding_prob'] > cpat_cutoff) & (cpat_orf_df['ORF'] <= 150)]['ID'])
    new_lncrna_micro_orfs = cuff_gtf[cuff_gtf['transcript_id'].isin(micro_orfs_ids)]

    # Concatenate and write output
    non_coding_gtf = pd.concat([new_linc_rna_non_coding, new_intronic_rna_non_coding, new_antisense_rna_non_coding])
    non_coding_gtf.to_csv(output_non_coding, sep='\t', header=False, index=False)

    coding_gtf = pd.concat([new_isoforms_coding, new_lncrna_micro_orfs])
    coding_gtf.to_csv(output_coding, sep='\t', header=False, index=False)


if __name__ == '__main__':
    cuffcompare_gtf = sys.argv[1]
    cpat_orf = sys.argv[2]
    cpat_no_orf = sys.argv[3]
    cutoff = sys.argv[4]
    output_non_coding = sys.argv[5]
    output_coding = sys.argv[6]
    classify_transcripts(cuffcompare_gtf, cpat_orf, cpat_no_orf, cutoff, output_non_coding, output_coding)

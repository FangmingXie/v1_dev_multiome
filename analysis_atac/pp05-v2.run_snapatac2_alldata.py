import snapatac2 as snap
import os

ddir = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/data/v1_multiome/atac_fragments'

# peak_file = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results_atac/all_l23_peaks_sal.bed'

samples = [
#    'P6a', 'P6b', 'P6c',
#     'P8a', 'P8b', 'P8c',
#     'P10a', 'P10b', 
#     'P12a', 'P12b', 'P12c', 
#     'P14a', 'P14b', 
#     'P17a', 'P17b', 
#     'P21a', 'P21b', 

#     'P12DRa', 'P12DRb', 
#     'P14DRa', 'P14DRb', 
#     'P17DRa', 'P17DRb', 
#     'P21DRa', 'P21DRb',
    ]


for sample in samples:
    fragment_file        = ddir+f'/raw/{sample}/atac_fragments.tsv.gz'
    if not os.path.isfile(fragment_file):
        fragment_file    = ddir+f'/raw/{sample}_reseq/atac_fragments.tsv.gz'
        assert os.path.isfile(fragment_file)

    print(sample, flush=True)
    output_file          = ddir+f'/frag_snap/ATAC_{sample}.h5ad'
    # output_peak_mat_file = ddir+f'/peak_mat_{sample}.h5ad'VDv
    
    print(f'generating {sample} fragment file', flush=True)
    data = snap.pp.import_data(
        fragment_file,
        chrom_sizes=snap.genome.mm10,
        sorted_by_barcode=False,
        file=output_file,  # Optional
    )

    # print(f'generating {sample} peak matrix', flush=True)
    # peak_mat = snap.pp.make_peak_matrix(data, peak_file=peak_file, file=output_peak_mat_file)
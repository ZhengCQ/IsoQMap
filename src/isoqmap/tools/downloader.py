import os
import urllib

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data')
ANNOTATION_URL = 'https://yourserver.org/path/to/transcript_gene_info.tsv.gz'

def get_annotation_file(species, reference_source):
    # Validate species input
    if species.lower() not in ['mouse', 'human']:
        print("The species should be Mouse or Human")
        return None

    # Validate reference_source input
    if reference_source.lower() not in ['gencode', 'ncbi']:
        print("The reference source should be 'gencode' or 'ncbi'")
        return None
    
    # Define the local file path based on species and reference source
    if species.lower() == 'mouse':
        if reference_source.lower() == 'gencode':
            LOCAL_FILENAME = os.path.join(DATA_DIR, 'GRCm39_M36/transcript_gene_info.tsv.gz')
        elif reference_source.lower() == 'ncbi':
            LOCAL_FILENAME = os.path.join(DATA_DIR, 'NCBI_Mm10/transcript_gene_info.tsv.gz')
    elif species.lower() == 'human':
        if reference_source.lower() == 'gencode':
            LOCAL_FILENAME = os.path.join(DATA_DIR, 'GRCh38_V42/transcript_gene_info.tsv.gz')
        elif reference_source.lower() == 'ncbi':
            LOCAL_FILENAME = os.path.join(DATA_DIR, 'NCBI_Hg38/transcript_gene_info.tsv.gz')

    # Ensure the annotation file is available and return its path
    if not os.path.exists(LOCAL_FILENAME):
        os.makedirs(DATA_DIR, exist_ok=True)
        print(f"Downloading annotation data to {LOCAL_FILENAME} ...")
        urllib.request.urlretrieve(ANNOTATION_URL, LOCAL_FILENAME)

    return LOCAL_FILENAME

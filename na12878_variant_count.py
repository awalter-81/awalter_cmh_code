# Imported libraries
import pandas as pd
import gzip

# Same function used in data visualization to extract headers from VCF file
def read_vcf_headers(file_path):
    headers = None
    with gzip.open(file_path, 'rt') as file:
        for line in file:
            if line.startswith('#CHROM'):
                headers = line.strip().lstrip('#').split('\t')
                break
    if headers is None:
        raise ValueError(f"The file {file_path} does not contain a header line starting with '#CHROM'.")
    return headers

# Same function used in data visualization to read the VCF file data and include extracted headers
def read_vcf(file_path, headers):
    with gzip.open(file_path, 'rt') as file:
        return pd.read_csv(file, comment='#', sep='\t', header=None, names=headers)

# Function to identify and filter for variants with an AF threshold of less than 0.01
# Added checks to account for and remove variant rows that do not have AF values or are missing AF altogether
def filter_vcf_by_af(file_path, af_threshold=0.01):
    headers = read_vcf_headers(file_path)
    vcf_data = read_vcf(file_path, headers)
    def has_valid_af(info_str):
        try:
            info_dict = {key: value for key, value in (item.split('=') for item in info_str.split(';') if '=' in item)}
            if 'AF' in info_dict:
                return float(info_dict['AF']) < af_threshold
            return False
        except ValueError:
            return False
    filtered_data = vcf_data[vcf_data['INFO'].apply(has_valid_af)]
    return filtered_data

# VCF files
file_path1 = "NA12878.chr21.slice.vcf.gz"
file_path2 = "gnomad.chr21.slice.vcf.gz"

# Step 1: Filter gnomad.chr21.slice.vcf.gz for AF < 0.01
gnomad_filtered = filter_vcf_by_af(file_path2)

# Extract the POS values from the filtered gnomad data
filtered_pos = gnomad_filtered['POS'].astype(str).tolist()

# Step 2: Filter NA12878.chr21.slice.vcf.gz based on matching POS values
headers1 = read_vcf_headers(file_path1)
na12878_data = read_vcf(file_path1, headers1)
na12878_filtered = na12878_data[na12878_data['POS'].astype(str).isin(filtered_pos)]

# Display the count number for variants < 0.01 in the filtered DataFrame
num_rows = na12878_filtered.shape[0]
print(f"NA12878 variant count < 0.01: {num_rows}")

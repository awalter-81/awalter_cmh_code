# Imported libraries
import pandas as pd
import gzip
import argparse

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

# command line arguments with flags for input files, allele frequency threshold, and filtered output variants
def main():
    parser = argparse.ArgumentParser(description="Filter VCF files and export matching entries to a CSV file.")
    parser.add_argument("-f1", "--file1", type=str, required=True, help="Path to the first VCF gzipped file (e.g., NA12878.chr21.slice.vcf.gz)")
    parser.add_argument("-f2", "--file2", type=str, required=True, help="Path to the second VCF gzipped file (e.g., gnomad.chr21.slice.vcf.gz)")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file name (e.g., filtered_NA12878.csv)")
    parser.add_argument("-af", "--af_threshold", type=float, default=0.01, help="Allele Frequency threshold for filtering (default: 0.01)")
    
    args = parser.parse_args()

    # Step 1: Filter the second VCF file based on AF threshold
    gnomad_filtered = filter_vcf_by_af(args.file2, args.af_threshold)

    # Extract the POS values from the filtered gnomad data
    filtered_pos = gnomad_filtered['POS'].astype(str).tolist()

    # Step 2: Filter the first VCF file based on matching POS values
    headers1 = read_vcf_headers(args.file1)
    na12878_data = read_vcf(args.file1, headers1)
    na12878_filtered = na12878_data[na12878_data['POS'].astype(str).isin(filtered_pos)]

    # Write the filtered data to a new CSV file
    na12878_filtered.to_csv(args.output, index=False)

    # Display the number of rows in the filtered DataFrame (excluding headers)
    num_rows = na12878_filtered.shape[0]
    print(f"Number of rows in the filtered NA12878 data: {num_rows}")
    print(f"Filtered data written to {args.output}")

if __name__ == "__main__":
    main()

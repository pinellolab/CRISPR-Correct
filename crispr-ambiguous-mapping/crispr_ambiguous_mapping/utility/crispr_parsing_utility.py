import re

def test_header_barcode_extraction(fastq_test_header: str, header_barcode_regex: str):
    # Use re.search to find the barcode match in the header
    barcode_match = re.search(header_barcode_regex, fastq_test_header)

    # Check if a barcode match is found
    if barcode_match:
        # Extract the barcode from the match
        barcode = barcode_match.group(1)
        print("Barcode:", barcode)
        return barcode
    else:
        print("No barcode match found.")


def test_header_umi_extraction(fastq_test_header: str, header_umi_regex: str):
    # Use re.search to find the umi match in the header
    umi_match = re.search(header_umi_regex, fastq_test_header)

    # Check if a umi match is found
    if umi_match:
        # Extract the umi from the match
        umi = umi_match.group(1)
        print("UMI:", umi)
        return umi
    else:
        print("No UMI match found.")
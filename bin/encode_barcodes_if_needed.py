#!/usr/bin/env python3

import gzip
import pandas as pd
from pathlib import Path

# Paths
input_path = Path("barcodes.tsv.gz")
output_path = Path("txd_input/barcodes_encoded.tsv.gz")
mapping_path = Path("txd_input/barcode_mapping.tsv")

# Read barcodes
with gzip.open(input_path, 'rt') as f:
    barcodes = [line.strip() for line in f]

# Check for 16bp + unique
all_16bp = all(len(bc) == 16 for bc in barcodes)
all_unique = len(set(barcodes)) == len(barcodes)

if all_16bp and all_unique:
    print("✅ Barcodes are already 16bp and unique. No encoding needed.")
else:
    print("⚠️ Barcodes are not all 16bp or contain duplicates. Encoding...")

    def encode_barcode(index):
        return f"CB{index:014d}"  # 16 characters total

    encoded_barcodes = [encode_barcode(i) for i in range(len(barcodes))]

    # Save encoded barcodes
    with gzip.open(output_path, 'wt') as f:
        for bc in encoded_barcodes:
            f.write(f"{bc}\n")

    # Save mapping
    mapping_df = pd.DataFrame({'original_barcode': barcodes, 'encoded_barcode': encoded_barcodes})
    mapping_df.to_csv(mapping_path, sep='\t', index=False)

    print(f"✅ Encoded {len(barcodes)} barcodes.")

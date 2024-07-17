#!/bin/bash
for f in ./*.barcodes.tsv; do
  # Extract the donor id from the filename and remove extra parts
  don=$(basename "$f" .barcodes.tsv | awk -F'.' '{print $1}')

  # Read the file, add the donor id, and append to the combined file
  awk -v donor_id="$don" 'BEGIN { OFS = "\t" } { print $0, donor_id }' "$f" >> donor_ids.tsv
done

# Add header to the combined file
sed -i '1i\cell\tdonor_id' donor_ids.tsv
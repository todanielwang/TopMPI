#Give me inputs to TopPIC in the format of location of TopPIC, location of database, location of ms2.msalign, and other parameters

#!/bin/bash

# Input file path (provided by user or script argument)
TopPIC="$1"

database="$2"

input_file="$3"

# Extract the base directory
base_dir=$(dirname "$input_file")

# Extract common prefix by trimming the last underscore segment and subsequent text
common_prefix=$(basename "$input_file" | sed -r 's/(_[^_]+)$//')

# New directory based on the modified prefix
new_sub_dir="${base_dir}/${common_prefix}_MSDeplex"

# Create the new directory if it doesn't exist
mkdir -p "$new_sub_dir"

# Define function to copy and rename files
copy_and_rename() {
    local src="$1"
    local dst="$2"
    if [ -f "$src" ]; then
        cp "$src" "$dst"
        echo "Copied and renamed: $(basename "$dst") from $(basename "$src")"
    else
        echo "File not found: $src"
    fi
}

# Generate new file names and copy them
for extension in feature.xml ms1.feature ms1.msalign ms2.feature ms2.msalign; do
    src_file="${base_dir}/${common_prefix}_${extension}"
    dst_file="${new_sub_dir}/A_${extension}"
    copy_and_rename "$src_file" "$dst_file"
done

python3 switchPrecursor.py ${new_sub_dir}/A_ms2.msalign

mv ${new_sub_dir}/A_ms2_modified.msalign ${new_sub_dir}/B_ms2.msalign

for extension in feature.xml ms1.feature ms1.msalign ms2.feature; do
    src_file="${new_sub_dir}/A_${extension}"
    dst_file="${new_sub_dir}/B_${extension}"
    copy_and_rename "$src_file" "$dst_file"
done

$TopPIC $database "${new_sub_dir}/A_ms2.msalign" "${@:4}" -f C57 -d -v 100000 -V 100000 -u 16 -K

$TopPIC $database "${new_sub_dir}/B_ms2.msalign" "${@:4}" -f C57 -d -v 100000 -V 100000 -u 16 -K

python3 noisefilter.py ${new_sub_dir}/

mv ${new_sub_dir}/FirstPrSM_ms2_modified.msalign ${new_sub_dir}/FirstPrSM_ms2.msalign

mv ${new_sub_dir}/SecondPrSM_ms2_modified.msalign ${new_sub_dir}/SecondPrSM_ms2.msalign

for extension in feature.xml ms1.feature ms1.msalign ms2.feature; do
    src_file="${new_sub_dir}/SecondPrSM_${extension}"
    dst_file="${new_sub_dir}/SecondPrSM_${extension}"
    copy_and_rename "$src_file" "$dst_file"
done

$TopPIC $database "${new_sub_dir}/SecondPrSM_ms2.msalign" "${@:4}" -f C57 -d -t FDR -T FDR -u 16


#Give me inputs to TopMPI in the format of location of TopPIC, location of database, location of ms2.msalign, TopPIC parameters, and TopMPI parameters

#!/bin/bash

# Ensure at least 3 positional arguments are provided
if [[ "$#" -lt 3 ]]; then
    echo "Usage: $0 <TopPIC> <fasta file> <ms2.msalign> [--TopPIC-flag=<flags>] [--TopMPI-flag=<flags>]"
    exit 1
fi

# Read the first three positional arguments
TopPIC="$1"
database="$2"
input_file="$3"
shift 3 # Shift the positional arguments off the list

# Initialize optional flags for Program 1 and Program 2
TopPIC_flags=""
TopMPI_flags=""

# Parse remaining optional arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --prog1-flag=*) # Optional flags for Program 1
            TopPIC_flags+=" ${1#*=}"
            ;;
        --prog2-flag=*) # Optional flags for Program 2
            TopMPI_flags+=" ${1#*=}"
            ;;
        *) # Catch all for unrecognized arguments
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
    shift
done

# # Display the parsed arguments (optional)
# echo "Positional Argument 1: $pos1"
# echo "Positional Argument 2: $pos2"
# echo "Positional Argument 3: $pos3"

# # Run Program 1 with its optional flags
# echo "Running Program 1 with flags: $prog1_flags"
# program1 "$pos1" "$pos2" "$pos3" $prog1_flags

# # Run Program 2 with its optional flags
# echo "Running Program 2 with flags: $prog2_flags"
# program2 "$pos1" "$pos2" "$pos3" $prog2_flags


# Extract the base directory
base_dir=$(dirname "$input_file")

# Extract common prefix by trimming the last underscore segment and subsequent text
common_prefix=$(basename "$input_file" | sed -r 's/(_[^_]+)$//')

# New directory based on the modified prefix
new_sub_dir="${base_dir}/${common_prefix}_TopMPI"

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

# Process Program 1 flags
if [[ "$prog1_flags" == *"-t FDR"* ]]; then
    echo "Ignoring '-t FDR' for Program 1."
    prog1_flags=$(echo "$prog1_flags" | sed 's/-t FDR//g')
fi

if [[ "$prog1_flags" == *"-v 0.01"* ]]; then
    echo "Ignoring '-v 0.01' for Program 1."
    prog1_flags=$(echo "$prog1_flags" | sed 's/-v 0.01//g')
fi

if [[ "$prog1_flags" == *"-d"* ]]; then
    echo "Adding '-K' because '-d' is provided for Program 1."
    prog1_flags+=" -K"
fi

# Trim any extra spaces from prog1_flags
prog1_flags=$(echo "$prog1_flags" | xargs)

copy_and_rename "${base_dir}/${common_prefix}_ms2.msalign" "${new_sub_dir}/First_ms2.msalign"

# Generate new file names and copy them
if [[ "$prog1_flags" != *"-x"* ]]; then
    for extension in feature.xml ms1.feature ms1.msalign ms2.feature; do
        src_file="${base_dir}/${common_prefix}_${extension}"
        dst_file="${new_sub_dir}/First_${extension}"
        copy_and_rename "$src_file" "$dst_file"
    done
fi

python3 switchPrecursor.py ${new_sub_dir}/First_ms2.msalign

mv ${new_sub_dir}/First_ms2_modified.msalign ${new_sub_dir}/Second_ms2.msalign

if [[ "$prog1_flags" != *"-x"* ]]; then
    for extension in feature.xml ms1.feature ms1.msalign ms2.feature; do
        src_file="${new_sub_dir}/First_${extension}"
        dst_file="${new_sub_dir}/Second_${extension}"
        copy_and_rename "$src_file" "$dst_file"
    done
fi

if [[ "$prog1_flags" == *"-d"* ]]; then
    echo "Adding '-K' because '-d' is provided for Program 1."
    TopPIC_flags+=" -K"
fi

temprunflags="${TopPIC_flags}"

if [[ "$prog1_flags" == *"-t FDR"* || "$prog1_flags" == *"-T FDR"* ]]; then
    echo "Ignoring '-t FDR' and '-T FDR' for Program 1."
    temprunflags=$(echo "$temprunflags" | sed 's/-t FDR//g; s/-T FDR//g')
fi

if [[ "$prog1_flags" == *"-v "* || "$prog1_flags" == *"-V "* ]]; then
    echo "Ignoring '-v <value>' and '-V <value>' for Program 1."
    temprunflags=$(echo "$temprunflags" | sed 's/-v [^ ]*//g; s/-V [^ ]*//g')
fi

temprunflags=$(echo "$temprunflags" | xargs)

$TopPIC $database "${new_sub_dir}/First_ms2.msalign" -v 100000 -V 100000 $temprunflags

$TopPIC $database "${new_sub_dir}/Second_ms2.msalign" -v 100000 -V 100000 $temprunflags

python3 checkAndRemovePeaks.py ${new_sub_dir}/ $TopMPI_flags

mv ${new_sub_dir}/Primary_ms2_modified.msalign ${new_sub_dir}/Primary_ms2.msalign

mv ${new_sub_dir}/Secondary_ms2_modified.msalign ${new_sub_dir}/Secondary_ms2.msalign

for extension in feature.xml ms1.feature ms1.msalign ms2.feature; do
    src_file="${new_sub_dir}/Primary_${extension}"
    dst_file="${new_sub_dir}/Secondary_${extension}"
    copy_and_rename "$src_file" "$dst_file"
done

$TopPIC $database "${new_sub_dir}/Secondary_ms2.msalign" $TopPIC_flags
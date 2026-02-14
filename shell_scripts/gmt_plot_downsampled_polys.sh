#!/usr/bin/env bash


# Initialize variables
INPUT_FILE=""
CPT_RNG=""
PROJECTION="15c"  # Default projection: Mercator with 15cm width


# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -cpt_rng)
            if [[ -z "$2" ]]; then
                echo "Error: -cpt_rng requires a value (format: min/max/interval)" >&2
                echo "Usage: $0 input_file [-cpt_rng min/max/interval] [-project 4i]" >&2
                exit 1
            fi
            CPT_RNG="$2"
            shift 2
            ;;
        -project)
            if [[ -z "$2" ]]; then
                echo "Error: -project requires a projection specification (e.g., 15c, 12c, 15c)" >&2
                echo "Usage: $0 input_file [-cpt_rng min/max/interval] [-project 4i]" >&2
                exit 1
            fi
            PROJECTION="$2"
            shift 2
            ;;
        *)
            if [[ -z "$INPUT_FILE" ]]; then
                INPUT_FILE="$1"
                shift
            else
                echo "Error: Unexpected argument $1" >&2
                echo "Usage: $0 input_file [-cpt_rng min/max/interval] [-project 4i]" >&2
		echo " Example: "
		echo "   gmt_plot_downsampled_polys.sh 20250307_20250412_SPO_algcl_xs64_i32_r4_r_ll_meter_quadds_ll.dat.xypoly.gmt -cpt_rng -0.2/0.2/0.01 -project 4i"
                exit 1
            fi
            ;;
    esac
done


# Validate input file
if [[ -z "$INPUT_FILE" ]]; then
    echo "Usage: $0 input_file [-cpt_rng min/max/interval] [-project projection_spec]" >&2
    echo " Example: "
    echo "     > gmt_plot_downsampled_polys.sh 20250307_20250412_SPO_algcl_xs64_i32_r4_r_ll_meter_quadds_ll.dat.xypoly.gmt -cpt_rng -0.2/0.2/0.01 -project 4i"
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist!" >&2
    exit 1
fi


# Generate output filename: remove any existing .png extension first, then add .png
FILENAME=$(basename "$INPUT_FILE")
BASE_NAME="${FILENAME%.*}"
OUTPUT_IMAGE="${BASE_NAME}"


# Extract raw values using line number rule (every 6th line, starting at 1)
echo "Diagnosing raw_value distribution..."
RAW_VALUES=$(awk '
    NR % 6 == 1 {
        line_num = NR
        if (match($0, /(-?[0-9]+\.?[0-9]*)/, arr)) {
            if (arr[1] ~ /^-?[0-9]+(\.[0-9]+)?$/) {
                print arr[1]
            } else {
                print "Warning: Invalid numeric value in line " line_num > "/dev/stderr"
                empty_count++
                print ""
            }
        } else {
            print "Warning: No valid raw_value in line " line_num > "/dev/stderr"
            empty_count++
            print ""
        }
    }
    END {
        if (empty_count > 0) print empty_count > "/dev/stderr"
    }
' "$INPUT_FILE")

# Statistics processing
EMPTY_COUNT=$(awk 'END {print empty_count+0}' <(echo "$RAW_VALUES" 2>&1 | grep -oP '\d+' | tail -1))
TOTAL_COUNT=$(echo "$RAW_VALUES" | wc -l)
VALID_VALUES=$(echo "$RAW_VALUES" | awk '$0 != ""')
VALID_COUNT=$(echo "$VALID_VALUES" | wc -l)
ZERO_COUNT=$(echo "$VALID_VALUES" | awk '$1 == 0 {count++} END {print count+0}')
NONZERO_COUNT=$((VALID_COUNT - ZERO_COUNT))

# Handle min/max values
if [ $VALID_COUNT -gt 0 ]; then
    MIN_RAW=$(echo "$VALID_VALUES" | sort -n | head -1)
    MAX_RAW=$(echo "$VALID_VALUES" | sort -n | tail -1)
else
    MIN_RAW="N/A (no valid values)"
    MAX_RAW="N/A (no valid values)"
fi

#
echo "Raw value statistics:"
echo " - Total polygons: $TOTAL_COUNT"
echo " - Valid values: $VALID_COUNT"
echo " - Zero values: $ZERO_COUNT"
echo " - Non-zero values: $NONZERO_COUNT"
echo " - Min value: $MIN_RAW"
echo " - Max value: $MAX_RAW"

#
if [ $VALID_COUNT -eq 0 ]; then
    echo "Error: No valid raw_values found in file - cannot generate color palette" >&2
    exit 1
fi

if [ $NONZERO_COUNT -eq 0 ]; then
    echo "Warning: All valid raw_values are zero - colors will be uniform" >&2
fi

if [ $EMPTY_COUNT -gt 0 ]; then
    echo "Error: Failed to extract $EMPTY_COUNT raw_values from file (see above for line numbers)" >&2
fi


# Get region from input file with higher precision
# Use -I0 to preserve original precision (no rounding)
RAW_REGION=$(gmt info -I0 "$INPUT_FILE" | awk '{print substr($1, 2)}' | sed 's/[^0-9/.+-]//g')
IFS='/' read -r WEST EAST SOUTH NORTH <<< "$RAW_REGION"


# Longitude correction function (preserve decimals)
fix_extreme_longitude() {
    local lon=$1
    lon=$(echo "$lon" | sed 's/[^0-9.+-]//g')
    
    if [ -z "$lon" ]; then
        echo "Warning: Empty longitude, using default" >&2
        echo "-10.0"
        return
    fi
    
    # Use bc for floating-point arithmetic to preserve decimals
    if (( $(echo "$lon > 180" | bc -l) )); then
        lon=$(echo "$lon - 360" | bc -l)
    elif (( $(echo "$lon < -180" | bc -l) )); then
        lon=$(echo "$lon + 360" | bc -l)
    fi
    
    if (( $(echo "$lon > 180 || $lon < -180" | bc -l) )); then
        echo "Warning: Invalid longitude, using default" >&2
        lon="-10.0"
    fi
    echo "$lon"
}


# Fix coordinates with floating-point handling to preserve decimals
WEST=$(fix_extreme_longitude "$WEST")
EAST=$(fix_extreme_longitude "$EAST")

# Use bc to limit latitude to [-90, 90] while preserving decimals
SOUTH=$(echo "scale=6; if ($SOUTH < -90) -90 else if ($SOUTH > 90) 90 else $SOUTH" | bc -l)
NORTH=$(echo "scale=6; if ($NORTH < -90) -90 else if ($NORTH > 90) 90 else $NORTH" | bc -l)

REGION="$WEST/$EAST/$SOUTH/$NORTH"
echo "Using valid region: $REGION"
echo "Using projection: $PROJECTION"


# Determine color palette range
if [[ -n "$CPT_RNG" ]]; then
    if [[ "$CPT_RNG" != *\/*\/* ]]; then
        echo "Error: Invalid -cpt_rng format. Use 'min/max/interval' (e.g., -0.03/0.03/0.005)" >&2
        exit 1
    fi
    echo "Creating color palette with user-specified range: $CPT_RNG"
    CPT_ARGS="$CPT_RNG"
else
    echo "Creating color palette with data-derived range: $MIN_RAW/$MAX_RAW/0.001"
    CPT_ARGS="$MIN_RAW/$MAX_RAW/0.001"
fi



# Plot map with data processing
echo "Plotting map..."
gmt begin "$OUTPUT_IMAGE" png,pdf


   gmt set MAP_FRAME_TYPE plain
   gmt set FONT_ANNOT_PRIMARY 10p
   gmt set FONT_LABEL 12p


   gmt basemap -R$REGION -JM$PROJECTION -Baf -BWsNe

   gmt makecpt -Cjet -T"$CPT_ARGS" -D -Z
   #
   # Generate GMT format data (every 6 lines as a group, first line as -Z value)
   gmt psxy $INPUT_FILE -C -G+z -W0.1p,white -L


   # Add color bar
   gmt colorbar -DJBC+o0/0.75c+w${PROJECTION}/0.5c+ml+h -Baf -Bx+l"LOS deformation(m)" 


gmt end


# Cleanup

echo "Plotting completed. Output: $OUTPUT_IMAGE.png and $OUTPUT_IMAGE.pdf"

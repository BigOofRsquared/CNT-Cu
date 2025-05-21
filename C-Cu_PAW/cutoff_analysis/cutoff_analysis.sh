if [ $# -lt 1 ]; then
  echo "specify filename"
  exit 1
fi

virgin_name="$1"
virgin_inputfile="${virgin_name}.in"

# energy cutoff that will be used #
cutoff_list=(10 20 30 40 50 60 70 80)
# # # # # # # # # # # # # # # # # #

for cutoff in "${cutoff_list[@]}"; do
    rhocutoff=$((cutoff * 8))
    name="${virgin_name}_ecut_${cutoff}_Ry"
    inputfile="${name}.in"
    outputfile="${name}.out"

    sed \
        -e "s|ECUTWFC|${cutoff}|g" \
        -e "s|ECUTRHO|${rhocutoff}|g" \
        -e "s|OUTDIR|${name}|g" \
        "$virgin_inputfile" > "$inputfile"

    echo "# # # # CUTOFF = ${cutoff} Ry # # # #"

    pw.x < "$inputfile" > "$outputfile"

done

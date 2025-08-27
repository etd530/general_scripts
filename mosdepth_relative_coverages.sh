for file in *txt; do
sample=$(echo $file | cut -f1 -d'.');
cat $file | grep -v chrom | grep -v total | grep -v mitochondrion | grep -v scf_ | awk -v sample=$sample '{
            lines[NR] = $0
            col4[NR] = $4
            sum += $4
        }
        END {
            mean = sum / NR
            for (i = 1; i <= NR; i++) {
                split(lines[i], fields)
                ratio = col4[i] / mean
                print sample "\t" fields[1] "\t" col4[i] "\t" ratio
            }
        }';
done
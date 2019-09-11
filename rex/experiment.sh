#!/usr/bin/env bash

./intersection.py ../data/blast.1.tsv ../data/blast.2.tsv | ./eutil.py efetch - > library.fna

glsearch36 -T 4 -m 8C ../data/query.1.fna library.fna > glsearch.1.tsv
glsearch36 -T 4 -m 8C ../data/query.2.fna library.fna > glsearch.2.tsv

./extract.R glsearch.1.tsv glsearch.2.tsv > range.tsv 2> /dev/null

samtools faidx -r range.tsv library.fna | sed 's/:.*//g' > rex.fna

cut -f 1 -d : range.tsv | ./eutil.py esummary - -params "retmode=json" | \
	tee docsum.json | \
	jq -r ".result | del(.uids) | map([.accessionversion, .subtype, .subname] | @tsv) | .[]" | \
	./subtype.R - 2> /dev/null | \
	./lubridate.R - collection_date > subtype.tsv 2> /dev/null

field="$(head -n 1 subtype.tsv | tr '\t' '\n' | cat -n | grep collection_date | tr -d ' ' | cut -f 1)"
cut -f "1,$field" subtype.tsv | awk 'NR > 1 && $2 != NA { print $1; }'

cut -f "1,$field" subtype.tsv | awk 'NR > 1 && $2 != "NA" { printf "/^>/ s/%s/%s_%s/g\n", $1, $1, $2; }' > rec.sed
cut -f "1,$field" subtype.tsv | awk 'NR > 1 && $2 != "NA" { print $1; }' | \
	xargs samtools faidx rex.fna | sed -f rec.sed > rec.fna

./clusters.R rec.fna 2> /dev/null | \
while read -r ele; do \
	msa="${ele/.fna/.msa.fna}"
	log="${ele/.fna/.log}"
	echo "$ele -> ${ele/.fna/.msa}"
	mafft --auto --adjustdirection --thread -1 "$ele" 2> "$log" | sed 's/^>_R_/>/g' > "$msa" && \
	iqtree -s "$msa" -pre "${ele/.fna/.phy}" -m TESTONLY
done

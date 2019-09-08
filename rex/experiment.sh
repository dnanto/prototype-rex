#!/usr/bin/env bash

./intersection.py ../data/blast.1.tsv ../data/blast.2.tsv | ./eutil.py efetch - -params > library.fna

glsearch36 -T 4 -m 8C ../data/query.1.fna library.fna > glsearch.1.tsv
glsearch36 -T 4 -m 8C ../data/query.2.fna library.fna > glsearch.2.tsv

./extract.py glsearch.1.tsv glsearch.2.tsv > range.tsv

samtools faidx -r range.tsv library.fna | sed 's/:.*//g' > rec.fna

cut -f 1 -d : range.tsv | ./eutil.py esummary - -params "retmode=json" | \
	tee docsum.json | \
	jq -r ".result | del(.uids) | map([.accessionversion, .subtype, .subname] | @tsv) | .[]" | \
	./subtype.R - 2> /dev/null | \
	./lubridate.R - collection_date > subtype.tsv 2> /dev/null

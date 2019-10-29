#!/usr/bin/env bash

cd "$(dirname "$0")" || exit

bgzip -c -d ops.fna.gz | \
	makeblastdb \
		-parse_seqids -hash_index \
		-blastdb_version 5 -dbtype nucl \
		-title ops -out ops -logfile ops.log -taxid_map ops.ssv

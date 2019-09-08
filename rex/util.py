#!/usr/bin/env python3

from collections import OrderedDict


def batchify(entries, size=10):
	batch = []
	for i, e in enumerate(entries, start=1):
		batch.append(e)
		if i % size == 0:
			yield batch
			batch = []

	if batch:
		yield batch


def parse_outfmt7(file):
	fields = []
	for line in map(str.strip, file):
		if line.startswith("# Fields: "):
			fields = line[10:].split(", ")
		elif line and not line.startswith("#"):
			yield OrderedDict(zip(fields, line.split("\t")))

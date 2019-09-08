#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from functools import reduce

import pandas as pd

from rex.util import parse_outfmt7


def parse_args(argv):
	parser = ArgumentParser(description="extract", formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("hits", nargs="+", help="the glsearch hits files")
	args = parser.parse_args(argv)
	return args


def main(argv):
	args = parse_args(argv[1:])

	fields = ["subject id", "s. start", "s. end", "% identity"]

	df_list = []
	for path in args.hits:
		with open(path) as file:
			df = pd.DataFrame(parse_outfmt7(file))
			for field in fields[1:]:
				df[field] = pd.to_numeric(df[field])
			df = df.loc[df.groupby(fields[0])[fields[-1]].idxmax()]
			df_list.append(df[fields[:-1]])

	df = reduce(lambda x, y: pd.merge(x, y, on=fields[0]), df_list)
	df["min"] = df.min(1)
	df["max"] = df.max(1)

	rangeify = lambda row: "%s:%d-%d" % (row[fields[0]], row["min"], row["max"])
	print(*df.apply(rangeify, axis=1).values, sep="\n")

	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))

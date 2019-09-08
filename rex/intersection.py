#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from operator import itemgetter

from rex.util import parse_outfmt7


def parse_args(argv):
	parser = ArgumentParser(description="intersection", formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("hits", nargs="+", help="the blast hits files")
	parser.add_argument("-saccver", default="subject acc.ver", help="the subject field key")
	args = parser.parse_args(argv)
	return args


def main(argv):
	args = parse_args(argv[1:])

	accs = set()
	getter = itemgetter(args.saccver)
	for path in args.hits:
		with open(path) as file:
			if accs:
				accs &= set(map(getter, parse_outfmt7(file)))
			else:
				accs = set(map(getter, parse_outfmt7(file)))

	print(*accs, sep="\n")

	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))

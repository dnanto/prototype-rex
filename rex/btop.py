import sys

# https://www.ncbi.nlm.nih.gov/books/NBK279682/#_cookbook_Traceback_operations_BTOP_
# BTOP operations consist of 
# 1.) a number with a count of matching letters, 
# 2.) two letters showing a mismatch (e.g., “AG” means A was replaced by G), or 
# 3.) a dash (“-“) and a letter showing a gap. 

def mutations(btop):
	pos, digis, chars = 0, [], []
	for ele in btop:
		if ele.isdigit():
			digis.append(ele)
		else:
			chars.append(ele)
		if chars and len(chars) % 2 == 0:
			if digis:
				pos += int("".join(digis))
				digis = []
			pos += 1
			yield (*chars, pos)
			chars = []

btop = sys.argv[1]

print(*mutations(btop), sep="\n")

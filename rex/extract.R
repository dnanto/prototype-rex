#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)

read_outfmt7 <- function(path)
{
	lines <- read_lines(path)
	fields <-
		str_split_fixed(lines[grep("^# Fields: ", lines)], " ", 3)[3] %>%
		str_split(", ", simplify = T)
	read_tsv(lines, col_names = fields, comment = "#")
}

process_hits <- function(df)
{
	select(df, `subject id`, `s. start`, `s. end`, `% identity`) %>%
		group_by(`subject id`) %>%
		top_n(1, `% identity`) %>%
		ungroup() %>%
		select(-`% identity`)
}

lapply(args, read_outfmt7) %>%
	lapply(process_hits) %>%
	Reduce(function(...) merge(..., by = "subject id"), .) %>%
	mutate(
		start = apply(.[tail(names(.), -1)], 1, min),
		end = apply(.[tail(names(.), -1)], 1, max)
	) %>%
	mutate(range = paste(`subject id`, ":", start, "-", end, sep = "")) %>%
	pull(range) %>%
	cat(sep = "\n")

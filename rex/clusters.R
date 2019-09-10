#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)

path.ani <- args[1]
path.fna <- args[2]

read_dist <- function(path, sep = "\t")
{
	lines <- read_lines(path)
	n <- as.integer(lines[1])
	tokens <- str_split_fixed(lines, sep, 2)
	names <- tail(tokens[ , 1], -1)
	tail(tokens[ , 2], -1) %>%
		str_split_fixed("\t", n) %>%
		apply(2, as.numeric) %>%
		as.data.frame(row.names = names) %>%
		set_names(names)
}

coor <- read_dist(path.ani) %>% as.dist() %>% cmdscale()
devnull <- capture.output(res <- NbClust::NbClust(coor, method = "kmeans"))
k <- max(res$Best.partition)

# as.data.frame(coor) %>%
# 	mutate(class = as.factor(res$Best.partition)) %>%
# 	ggplot(aes(V1, V2)) +
# 	geom_point(aes(color = class), alpha = 0.5) +
# 	theme_minimal() +
# 	theme(legend.position = "bottom")

rec <- ape::read.FASTA(path.fna)
width <- nchar(as.character(k))
batches <- split(rec, res$Best.partition)
for (k in names(batches))
{
	path <- sprintf("k-%0*d.fna", width, as.integer(k))
	ape::write.FASTA(batches[[k]], path)
	cat(path, fill = T)
}

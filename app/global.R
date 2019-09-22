library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(tidyverse)

roots <- c(home = "../data")

datatable <- function(df)
{
	DT::datatable(
		df,
		rownames = FALSE,
		style = "bootstrap",
		class = "table-bordered table-striped table-hover responsive",
		filter = list(position = "top")
	)
}

read_header <- function(path)
{
	read_lines(path) %>%
		.[grep("^>", .)] %>%
		str_remove("^>") %>%
		str_split_fixed(" ", 2) %>%
		as.data.frame() %>%
		set_names(c("id", "definition"))
}

list_blastdb <- function(path)
{
	system(paste("blastdbcmd", "-list", path, "-list_outfmt", "'%f    %p    %t    %d    %l    %n    %U'"), intern = T) %>%
		enframe(name = NULL) %>%
		separate(value, c("path", "type", "title", "update", "bases", "sequences", "bytes"), "    ", convert = T)
}

read_dist <- function(path)
{
	lines <- read_lines(path)
	n <- as.integer(lines[1])
	str_split_fixed(tail(lines, -1), "\t", n) %>%
		as_tibble() %>%
		column_to_rownames("V1") %>%
		add_column(" ") %>%
		set_names(rownames(.))
}

NbClust <- function(...)
{
	pdf(file = NULL)
	null <- capture.output(suppressWarnings(result <- NbClust::NbClust(...)))
	null <- capture.output(dev.off())
	result
}

mash <- function(path, k = 21)
{
	cmd <- paste("mash", "triangle", "-k", k, path)
	read_dist(system(cmd, intern = T))
}

db <- str_split(Sys.getenv("BLASTDB"), ":", simplify = T) %>% apply(2, list_blastdb) %>% bind_rows()


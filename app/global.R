library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(tidyverse)

roots <- c(home = "~")

np <- parallel::detectCores()

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

rec_info <- function(path)
{
	ape::read.FASTA(path) %>%
		sapply(length) %>%
		enframe(value = "length")
}

list_blastdb <- function(path)
{
	system(paste("blastdbcmd", "-list", path, "-list_outfmt", "'%f    %p    %t    %d    %l    %n    %U'"), intern = T) %>%
		enframe(name = NULL) %>%
		separate(value, c("path", "type", "title", "update", "bases", "sequences", "bytes"), "    ", convert = T)
}

parse_outfmt7 <- function(lines)
{
	str_split_fixed(lines[grep("^# Fields: ", lines)], " ", 3)[3] %>%
		str_split(", ", simplify = T) %>%
		read_tsv(lines, comment = "#", col_names = .)
}

read_outfmt7 <- function(path)
{
	lines <- read_lines(path) %>% head(-1)
	enc <- startsWith(lines, "#") %>% rle()
	enc$values <- ceiling(seq_along(enc$values) / 2)
	split(lines, inverse.rle(enc)) %>%
		lapply(parse_outfmt7)
}

contextify <- function(df)
{
	group_by(df, `subject acc.ver`) %>%
		summarise(
			`s. start` = min(`s. start`, `s. end`), `s. end` = max(`s. start`, `s. end`),
			`query length` = head(`query length`, 1), `subject length` = head(`subject length`, 1)
		) %>%
		ungroup() %>%
		mutate(
			start = ifelse(`s. start` > `s. end`, `s. end`, `s. start`) - `query length`,
			end = ifelse(`s. start` > `s. end`, `s. start`, `s. end`) + `query length`
		) %>%
		mutate(
			start = ifelse(start < 0, 1, start),
			end = ifelse(end > `subject length`, `subject length`, end)
		)
}

blast <- function(qry, db, out, np = 1)
{
	cmd <- paste(
		Sys.which("blastn"),
		"-task", "megablast",
		"-query", qry,
		"-db", db,
		"-outfmt", "'7 std qcovhsp qlen slen sstrand btop'",
		"-out", out,
		"-num_threads", np
	)
	system(cmd)
}

blastdbcmd <- function(batch, db, out = "-", ...)
{
	cmd <- paste(Sys.which("blastdbcmd"), "-db", db, "-entry_batch", "-", "-out", out, ...)
	system(cmd, input = batch, intern = out == "-")
}

glsearch <- function(qry, db, out, np = 1)
{
	cmd <- paste(Sys.which("glsearch36"), "-T", np, "-m", "8CB", qry, db, ">", out)
	system(cmd)
}

decode_btop <- function(btop)
{
	matches <- as.integer(str_extract_all(btop, "\\d+", simplify = T))
	matches <- if (str_starts(btop, "\\d", negate = T)) c(0, matches) else matches
	pos <- qpos <- spos <- 0
	n <- 0; m <- 0; muts <- list();
	for (ele in Filter(nchar, str_split(btop, "\\d+", simplify = T)))
	{
		val <- matches[n<-n+1]
		pos <- pos + val
		qpos <- qpos + val
		spos <- spos + val
		for (i in seq(1, nchar(ele), 2))
		{
			tok <- substr(ele, i, i + 1)
			qpos <- qpos + str_starts(tok, "-", negate = T)
			spos <- spos + str_ends(tok, "-", negate = T)
			muts[[m<-m+1]] <- list(pos = pos<-pos+1, qpos = qpos, spos = spos, mut = tok)
		}
	}
	muts
}

call_mut <- function(mut)
{
	case_when(
		mut == "AC" ~ "trv", mut == "CA" ~ "trv",
		mut == "AG" ~ "trs", mut == "GA" ~ "trs",
		mut == "AT" ~ "trv", mut == "TA" ~ "trv",
		mut == "CG" ~ "trv", mut == "GC" ~ "trv",
		mut == "CT" ~ "trs", mut == "TC" ~ "trs",
		mut == "GT" ~ "trv", mut == "TG" ~ "trv",
		endsWith(mut, "-") ~ "ins", startsWith(mut, "-") ~ "del"
	)
}

plot_hits <- function(hits, snp)
{
	ggplot(hits, aes(x = `q. start`, y = `subject acc.ver`)) +
		geom_segment(aes(xend = `q. end`, yend = `subject acc.ver`, alpha = `% qidentity`), size = 2) +
		geom_point(data = snp, aes(color = call, shape = call, size = call)) +
		scale_color_manual(values = c(trv = "magenta", trs = "cyan", ins = "orange", del = "blue")) +
		scale_shape_manual(values = c(trv = "|", trs = "|", ins = "+", del = "x")) +
		scale_size_manual(values = c(trv = 4, trs = 4, ins = 8, del = 6)) +
		facet_wrap(~ `query acc.ver`, nrow = 1, scales = "free_x") +
		theme_minimal() +
		theme(
			legend.position = "bottom",
			text = element_text(family = "mono"),
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
		)
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

db <- flatten_chr(str_split(Sys.getenv("BLASTDB"), ":")) %>% lapply(list_blastdb) %>% bind_rows()

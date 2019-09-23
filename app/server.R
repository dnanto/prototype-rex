function(input, output, session) {

	# query

	## the database table

	output$db <- DT::renderDT(datatable(db))

	## the query file
	shinyFileChoose(input, "import_qry", roots = roots)
	path_qry <- eventReactive(input$import_qry, req(pull(parseFilePaths(roots, input$import_qry), datapath)))
	output$path_qry <- eventReactive(path_qry(), path_qry())

	## the query table
	output$query <- DT::renderDT(datatable(rec_info(path_qry())))

	bdb <- eventReactive(input$db, req(input$db))

	observeEvent(input$run, {
		req(path <- path_qry(), db <- bdb())

		withProgress({
			hits.1 <- file.path(dirname(path), "blast.tsv")
			hits.2 <- file.path(dirname(path), "glsearch.tsv")
			lib <- file.path(dirname(path), "library.fna")
			ext <- file.path(dirname(path), "extract.fna")

			incProgress(1/4, "blastn...")

			blast(path, db, hits.1, np)

			read_outfmt7(hits.1) %>%
				lapply(contextify) %>%
				Reduce(function(...) merge(..., by = "subject acc.ver"), .) %>%
				mutate(
					start = apply(.[names(.)[grep("^start", names(.))]], 1, min),
					end = apply(.[names(.)[grep("^end", names(.))]], 1, max)
				) %>%
				with(paste(`subject acc.ver`, " ", start, "-", end, sep = "")) %>%
				blastdbcmd(db, lib)

			incProgress(1/4, "glsearch...")

			glsearch(path, lib, hits.2, np)

			incProgress(1/4, "results...")

			read_outfmt7(hits.2) %>%
				lapply(select, `subject id`, `s. start`, `s. end`) %>%
				lapply(distinct, `subject id`, .keep_all = T) %>%
				Reduce(function(...) merge(..., by = "subject id"), .) %>%
				mutate(
					start = apply(.[names(.)[grep("^s. start", names(.))]], 1, min),
					end = apply(.[names(.)[grep("^s. end", names(.))]], 1, max)
				) %>%
				separate(`subject id`, c("subject acc.ver", "start.1", "end.1"), sep = "[-: ]", convert = T) %>%
				mutate(offset = start.1 - start) %>%
				mutate(start = start + offset, end = end + offset) %>%
				with(paste(`subject acc.ver`, " ", start, "-", end, sep = "")) %>%
				blastdbcmd(db, ext)

			hits1 <-
				read_outfmt7(hits.1) %>%
				bind_rows() %>%
				mutate(`% qidentity` = `% identity` * `% query coverage per hsp` / 100)

			snps1 <-
				apply(hits1, 1, function(row) {
					val <- decode_btop(row["BTOP"])
					if (!is_empty(val))
						bind_rows(val) %>%
						mutate(
							`q. start` = qpos + as.integer(row["q. start"]) - 1,
							`query acc.ver` = row["query acc.ver"],
							`subject acc.ver` = row["subject acc.ver"]
						)
				}) %>%
				bind_rows() %>%
				mutate(call = call_mut(mut))

			output$plt.hits.1 <- renderPlot(plot_hits(hits1, snps1))
			output$tbl.hits.1 <- DT::renderDT(datatable(hits1))
			output$tbl.snps.1 <- DT::renderDT(datatable(snps1))

			hits2 <-
				read_outfmt7(hits.2) %>%
				bind_rows() %>%
				distinct(`query id`, `subject id`, .keep_all = T) %>%
				mutate(
					qstart = ifelse(`q. start` < `q. end`, `q. start`, `q. end`),
					qend = ifelse(`q. start` < `q. end`, `q. end`, `q. start`)
				) %>%
				mutate(`q. start` = qstart, `q. end` = qend) %>%
				mutate(`% qidentity` = `% identity`) %>%
				rename(`query acc.ver` = `query id`, `subject acc.ver` = `subject id`)

			snps2 <-
				apply(hits2, 1, function(row) {
					val <- decode_btop(row["BTOP"])
					if (!is_empty(val))
						bind_rows(val) %>%
						mutate(
							`q. start` = qpos + as.integer(row["q. start"]) - 1,
							`query acc.ver` = row["query acc.ver"],
							`subject acc.ver` = row["subject acc.ver"]
						)
				}) %>%
				bind_rows() %>%
				mutate(call = call_mut(mut))

			output$plt.hits.2 <- renderPlot(plot_hits(hits2, snps2))
			output$tbl.hits.2 <- DT::renderDT(datatable(hits2))
			output$tbl.snps.2 <- DT::renderDT(datatable(snps2))

			output$extract <- DT::renderDT(datatable(rec_info(ext)))

			incProgress(1/4, "complete!")
		})
	})
}

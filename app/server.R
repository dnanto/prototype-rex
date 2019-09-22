function(input, output, session) {

	# query

	## the database table

	output$db <- DT::renderDT(datatable(db))

	## the query file
	shinyFileChoose(input, "import_qry", roots = roots)
	path_qry <- eventReactive(input$import_qry, req(pull(parseFilePaths(roots, input$import_qry), datapath)))
	output$path_qry <- eventReactive(path_qry(), path_qry())

	## the query table
	output$query <- DT::renderDT(datatable(read_header(path_qry())))

	bdb <- eventReactive(input$db, req(input$db))

	hits.1 <- eventReactive(input$run, {
		req(path <- path_qry(), db <- bdb())
		out <- file.path(dirname(path), "blast.tsv")
		blast(path, db, out, np)
		bind_rows(read_outfmt7(out)) %>%
			mutate(`% qidentity` = `% identity` * `% query coverage per hsp` / 100)
	})

	snps.1 <- eventReactive(hits.1(), {
		req(hits.1()) %>%
			apply(1, function(row) {
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
	})

	output$plt.hits.1 <- renderPlot(plot_hits(hits.1(), snps.1()))

}

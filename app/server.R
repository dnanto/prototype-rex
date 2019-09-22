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

}

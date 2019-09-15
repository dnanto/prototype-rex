shinyServer(function(input, output, session) {

	#### get the path to the FASTA file
	shinyFileChoose(input, "import_fas", roots = roots)
	path_fas <- eventReactive(input$import_fas, req(pull(parseFilePaths(roots, input$import_fas), datapath)))
	output$path_fas <- eventReactive(path_fas(), path_fas())

	#### output the taxa table
	output$taxa <- DT::renderDT(datatable(read_header(path_fas())))

	#### run
	observeEvent(input$run, {
		coor <- mash(path_fas(), k = input$k) %>% as.dist() %>% cmdscale() %>% scale()
		res <- NbClust(coor, method = "kmeans")
		output$taxa <- DT::renderDT(datatable(mutate(read_header(path_fas()), class = res$Best.partition)))
		output$plot <- renderPlot(
			as.data.frame(coor) %>%
				setNames(c("dim.1", "dim.2")) %>%
				mutate(class = as.factor(res$Best.partition)) %>%
				ggplot(aes(dim.1, dim.2)) +
				geom_point(aes(color = class), alpha = 0.5) +
				theme_minimal() +
				theme(legend.position = "bottom")
		)
	})

})

fluidPage(
	wellPanel(actionButton("run", "run", width = "100%")),
	tabsetPanel(
		tabPanel(
			"alignment",
			wellPanel(
				selectizeInput("db", "db", db$title)
			),
			DT::DTOutput("db"),
			wellPanel(
				shinyFilesButton("import_qry", "Import FASTA", "Please select a FASTA file...", multiple = F, icon = icon("file-import")),
				verbatimTextOutput("path_qry", placeholder = T)
			),
			DT::DTOutput("query")
		),
		tabPanel(
			"hits.1",
			plotOutput("plt.hits.1", height = 600),
			DT::DTOutput("tbl.hits.1"),
			DT::DTOutput("tbl.snps.1")
		),
		tabPanel(
			"hits.2",
			plotOutput("plt.hits.2", height = 600),
			DT::DTOutput("tbl.hits.2"),
			DT::DTOutput("tbl.snps.2")
		),
		tabPanel(
			"extract",
			DT::DTOutput("extract")
		)
	)
)

fluidPage(
	tabsetPanel(
		tabPanel(
			"database",
			wellPanel(
				selectizeInput("db", "db", db$title)
			),
			DT::DTOutput("db")
		),
		tabPanel(
			"query",
			wellPanel(
				shinyFilesButton("import_qry", "Import FASTA", "Please select a FASTA file...", multiple = F, icon = icon("file-import")),
				verbatimTextOutput("path_qry", placeholder = T)
			),
			DT::DTOutput("query")
		),
		tabPanel(
			"blast",
			wellPanel(
				actionButton("run_blast", "run_blast")
			),
			DT::DTOutput("hits_blast")
		)
	)
)

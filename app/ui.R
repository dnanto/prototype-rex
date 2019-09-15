fluidPage(
	wellPanel(
		shinyFilesButton("import_fas", "Import FASTA", "Please select a FASTA file...", multiple = F, icon = icon("file-import")),
		verbatimTextOutput("path_fas", placeholder = T),
		numericInput("k", "k", 5, min = 1, width = "100%"),
		actionButton("run", "run", width = "100%")
	),
	fluidRow(column(4, DT::DTOutput("taxa")), column(8, plotOutput("plot")))
)

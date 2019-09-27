navbarPage(
	"rex",
	tabPanel(
		"alignment",
		tabsetPanel(
			tabPanel(
				"params",
				wellPanel(actionButton("run", "run", width = "100%")),
				wellPanel(
					shinyFilesButton(
						"import_qry", "Import FASTA", "Please select a FASTA file...",
						multiple = F, icon = icon("file-import")
					),
					verbatimTextOutput("path_qry", placeholder = T)
				),
				DT::DTOutput("query"),
				wellPanel(
					selectizeInput("db", "db", db$title)
				),
				DT::DTOutput("db")
			),
			tabPanel(
				"blast",
				tabsetPanel(
					tabPanel("hits", DT::DTOutput("tbl.hits.1")),
					tabPanel("snp", DT::DTOutput("tbl.snps.1")),
					tabPanel("plot", plotOutput("plt.hits.1", height = 1000))
				)
			),
			tabPanel(
				"glsearch",
				tabsetPanel(
					tabPanel("hits", DT::DTOutput("tbl.hits.2")),
					tabPanel("snp", DT::DTOutput("tbl.snps.2")),
					tabPanel("plot", plotOutput("plt.hits.2", height = 1000))
				)
			),
			tabPanel(
				"extract",
				DT::DTOutput("extract")
			)
		)
	)
)

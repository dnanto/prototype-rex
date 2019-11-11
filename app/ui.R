navbarPage(
	"rex",
	tabPanel(
		"extract",
		tabsetPanel(
			tabPanel(
				"params",
				wellPanel(actionButton("run_extract", "extract", width = "100%")),
				wellPanel(
					shinyFilesButton(
						"import_qry", "Import FASTA", "Please select a FASTA file...",
						multiple = F, icon = icon("file-import")
					),
					verbatimTextOutput("path_qry", placeholder = T),
					DT::DTOutput("query")
				),
				wellPanel(
					selectizeInput("db", "db", db$title),
					DT::DTOutput("db")
				)
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
				"records",
				DT::DTOutput("result_ext")
			)
		)
	),
	tabPanel(
		"cluster",
		tabsetPanel(
			tabPanel(
				"params",
				wellPanel(actionButton("run_cluster", "cluster", width = "100%")),
				wellPanel(
					shinyFilesButton(
						"import_ext", "Import FASTA", "Please select a FASTA file...",
						multiple = F, icon = icon("file-import")
					),
					verbatimTextOutput("path_ext", placeholder = T),
					DT::DTOutput("extract")
				)
			),
			tabPanel(
				"clusters",
				plotOutput("plt.hclust", height = 750),
				DT::DTOutput("clusters")
			)
		)
	),
	tabPanel(
		"chrono",
		tabsetPanel(
			tabPanel(
				"params",
				wellPanel(actionButton("run_chrono", "chrono", width = "100%")),
				wellPanel(
					shinyFilesButton(
						"import_cls", "Import FASTA", "Please select a FASTA file...",
						multiple = F, icon = icon("file-import")
					),
					verbatimTextOutput("path_cls", placeholder = T),
					DT::DTOutput("cluster")
				)
			),
			tabPanel(
				"metadata",
				DT::DTOutput("metadata")
			),
			tabPanel(
				"signal",
				plotOutput("signal")
			),
			tabPanel(
				"MCMC",
				plotOutput("chronogram"),
				plotOutput("trace"),
				verbatimTextOutput("ESS")
			)
		)
	)
)

function(input, output, session) {

	shinyFileChoose(input, "import_qry", roots = roots)
	path_qry <- eventReactive(input$import_qry, req(pull(parseFilePaths(roots, input$import_qry), datapath)))
	output$path_qry <- eventReactive(path_qry(), path_qry())
	output$query <- DT::renderDT(datatable(rec_info(path_qry())))

	bdb <- eventReactive(input$db, req(input$db))
	output$db <- DT::renderDT(datatable(db))

	observeEvent(input$run_extract, {
		req(path <- path_qry(), db <- bdb())
		root <- dirname(path)
		withProgress({
			hits.1 <- file.path(root, "hits-1.tsv")
			hits.2 <- file.path(root, "hits-2.tsv")
			lib <- file.path(root, "lib.fna")
			ext <- file.path(root, "ext.fna")

			incProgress(1/4, "blastn...")

			blast(path, db, hits.1, np)

			entry <-
				read_outfmt7(hits.1) %>%
				lapply(contextify) %>%
				Reduce(function(...) merge(..., by = "subject acc.ver"), .) %>%
				mutate(
					start = apply(.[names(.)[grep("^start|end", names(.))]], 1, min),
					end = apply(.[names(.)[grep("^start|end", names(.))]], 1, max)
				)

			with(entry, paste(`subject acc.ver`, " ", start, "-", end, sep = "")) %>%
				blastdbcmd(db, "-", "-outfmt", "'%a    %t    %s'") %>%
				enframe(name = NULL) %>%
				separate(value, c("subject acc.ver", "title", "seq"), "    ") %>%
				merge(entry, by = "subject acc.ver") %>%
				mutate(accession = paste(`subject acc.ver`, ":", start, "-", end, sep = "")) %>%
				apply(1, function(row) paste(">", row["accession"], " ", row["title"], "\n", row["seq"], sep = "")) %>%
				write_lines(lib)

			incProgress(1/4, "glsearch...")

			glsearch(path, lib, hits.2, np)

			incProgress(1/4, "results...")

			entry <-
				read_outfmt7(hits.2) %>%
				lapply(select, `subject id`, `s. start`, `s. end`) %>%
				lapply(distinct, `subject id`, .keep_all = T) %>%
				Reduce(function(...) merge(..., by = "subject id"), .) %>%
				mutate(
					start = apply(.[names(.)[grep("start", names(.))]], 1, min),
					end = apply(.[names(.)[grep("end", names(.))]], 1, max)
				) %>%
				separate(`subject id`, c("subject acc.ver", "start.1", "end.1"), sep = "[-: ]", convert = T) %>%
				mutate(offset = start.1 - start) %>%
				mutate(start = start + offset, end = end + offset)

			rec <-
				with(entry, paste(`subject acc.ver`, " ", start, "-", end, sep = "")) %>%
				blastdbcmd(db, "-", "-outfmt", "'%a    %t    %s'") %>%
				enframe(name = NULL) %>%
				separate(value, c("subject acc.ver", "title", "seq"), "    ", convert = T) %>%
				merge(entry, by = "subject acc.ver") %>%
				mutate(
					accession = paste(`subject acc.ver`, " ", start, "-", end, sep = ""),
					length = end - start + 1
				)

			apply(rec, 1, function(row) paste(">", row["accession"], " ", row["title"], "\n", row["seq"], sep = "")) %>%
				write_lines(ext)

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
				mutate(
					`q. start` = qstart,
					`q. end` = qend,
					`% qidentity` = `% identity`,
					`query acc.ver` = `query id`,
					`subject acc.ver` = `subject id`
				)

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

			output$result_ext <- DT::renderDT(datatable(select(rec, `subject acc.ver`, start, end, title, length)))

			incProgress(1/4, "complete!")
		})
	})


	shinyFileChoose(input, "import_ext", roots = roots)
	path_ext <- eventReactive(input$import_ext, req(pull(parseFilePaths(roots, input$import_ext), datapath)))
	output$path_ext <- eventReactive(path_ext(), path_ext())
	output$extract <- DT::renderDT(datatable(rec_info(path_ext())))

	observeEvent(input$run_cluster, {
		req(path <- path_ext())
		root <- dirname(path)
		withProgress({
			incProgress(1/5, "mafft...")
			rec <- ape::read.FASTA(path)
			#msa <- ips::mafft(rec, method = "retree 1", option = c("--adjustdirection"), exec = Sys.which("mafft"))
			msa <- ape::read.FASTA("../data/msa.fna")
			names(msa) <- str_remove(names(msa), "^_R_")

			incProgress(1/5, "cluster...")
			dst <- ape::dist.dna(msa, model = "raw", pairwise.deletion = F)
			devnull <- capture.output(res <- NbClust(diss = dst, distance = NULL, method = "ward.D2", index = "silhouette", max.nc = 30))
			recoding <- arrange(enframe(table(res$Best.partition)), desc(value)) %>% with(setNames(seq_along(name), name))
			bp <- recode(res$Best.partition, !!!recoding)

			output$plt.hclust <- renderPlot(
				ape::as.phylo(hclust(dst)) %>%
					groupOTU(split(names(rec), bp)) %>%
					ggtree(aes(color = group), size = 2) +
					geom_tiplab()
			)

			output$clusters <- DT::renderDT(datatable(mutate(rec_info(path_ext()), cluster = bp)))

			incProgress(1/5, "batch...")
			nk <- max(res$Best.partition)
			width <- nchar(as.character(nk))
			batches <- split(rec, bp)

			for (k in 1:nk)
			{
				incProgress(0, sprintf("batch %0*d/%d...", width, k, nk))
				# xrec <- batches[[k]]
				path <- file.path(root, sprintf("k-%0*d.rec.fna", width, k))
				ape::write.dna(batches[[k]], path, format = "fasta", colsep = "")
			}

			incProgress(1/5, "complete!")
		})
	})


	shinyFileChoose(input, "import_cls", roots = roots)
	path_cls <- eventReactive(input$import_cls, req(pull(parseFilePaths(roots, input$import_cls), datapath)))
	output$path_cls <- eventReactive(path_cls(), path_cls())
	output$cluster <- DT::renderDT(datatable(rec_info(path_cls())))

	observeEvent(input$run_chrono, {
		req(path <- path_cls())
		root <- dirname(path)
		withProgress({
			incProgress(1/5, "metadata...")
			rec <- ape::read.FASTA(path)
			lab <- setNames(as.data.frame(str_split_fixed(names(rec), " ", 2)), c("acc", "title"))
			tag <- first(strsplit(basename(path), "\\.")[[1]])
			out <- file.path(root, paste(tag, "obj.json", sep = "."))
			# system(paste("../rex/eutil.py", "esummary", "-", "-params", "retmode=json", ">", out), input = lab$acc)
			meta <-
				jqr::jq(file(out), ".result | del(.uids) | map([.accessionversion, .subtype, .subname]) | .[]") %>%
				textConnection() %>%
				jsonlite::stream_in() %>%
				set_names(c("acc", "subtype", "subname")) %>%
				bind_cols(
					bind_rows(
						apply(., 1, function(row) {
							key <- str_split(row["subtype"], "\\|")[[1]]
							val <- str_split(row["subname"], "\\|")[[1]]
							data.frame(as.list(setNames(val, key)), stringsAsFactors = F)
						})
					)
				) %>%
				mutate_at("collection_date", lubridate::parse_date_time, orders = c("dbY", "Ymd", "bY", "Y"), quiet = T) %>%
				filter(complete.cases(collection_date)) %>%
				select(acc, collection_date, everything(), -subname, -subtype)

			output$metadata <- DT::renderDT(datatable(meta))


			incProgress(1/5, "mafft...")
			rec <- rec[meta$acc]
			names(rec) <- with(meta, paste(acc, strftime(collection_date, format = "%Y-%m-%d"), sep = "_"))
			msa <- ips::mafft(rec, method = "retree 1", option = c("--adjustdirection"), exec = Sys.which("mafft"))
			rownames(msa) <- str_remove(rownames(msa), "^_R_")
			out <- file.path(root, paste(tag, "msa.fna", sep = "."))
			ape::write.dna(msa, out, format = "fasta", colsep = "")


			incProgress(1/5, "ML...")
			pre <- file.path(root, paste(tag, "phy", sep = "."))
			system(paste("iqtree", "-pre", pre, "-s", out, "-m", "TESTONLY"))

			incProgress(1/5, "rtt...")
			phy <- file.path(root, paste(tag, "phy.treefile", sep = ".")) %>% ape::read.tree()
			tip.date <- phy$tip.label %>% strsplit("_") %>% sapply(last) %>% lubridate::ymd()
			phy <- initRoot(phy, as.numeric(tip.date))
			output$signal <- renderPlot(roottotip(phy, as.numeric(tip.date)))

			incProgress(1/5, "mcmc...")
			result <- bactdate(phy, as.numeric(tip.date), model = "relaxedgamma", nbIts = 1000000, thin = 1000)
			output$chronogram <- renderPlot(plot(result, "treeCI"))
			output$trace <- renderPlot(plot(result, "trace"))
			output$ESS <- renderText(coda::effectiveSize(coda::as.mcmc(result)))
		})
	})

}

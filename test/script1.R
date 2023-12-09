        dat <- phys_data_ready[[1]]
        attribute_group <- dat$Attribute_group |>
            {\(y) y[!is.na(y)]}() |>
            unique()
        attribute_type <- dat$Attribute_type |>
            {\(y) y[!is.na(y)]}() |>
            unique()
        attribute_nms <- dat$Attribute |>
            {\(y) y[!is.na(y)]}() |>
            unique()
        
        dat_n_tax <- length(unique(dat$NCBI_ID))
        node_list <- split(dat, factor(dat$NCBI_ID))
        
        not_in_ncbi_tree <- all(!names(node_list) %in% ncbi_nodes)
        
        if (not_in_ncbi_tree) {
            msg <- paste0(
                'Not enough data for propagation for ', phys_name, '. Rank:', rank_var, '.',
                ' Quitting R script. Finishing at ', Sys.time(), '.'
            )
            log_print(msg, blank_after = TRUE)
            quit(save = 'no')
        }
        
        # ncbi_tree$Do(function(node) {
        #     if (node$name %in% names(node_list))
        #         node$attribute_tbl <- node_list[[node$name]]
        # })
        # ncbi_tree$Do(
        #     function(node) {
        #         taxPool(
        #             node = node,
        #             grp = attribute_group,
        #             typ = attribute_type)
        #     },
        #     traversal = 'post-order'
        # )
        # ncbi_tree$Do(inh1, traversal = 'pre-order')
        # new_dat <- ncbi_tree$Get(
        #     'attribute_tbl', filterFun = function(node) {
        #         grepl('^[gst]__', node$name)
        #     }
        # ) |>
        #     discard(~ all(is.na(.x))) |>
        #     bind_rows() |>
        #     arrange(NCBI_ID, Attribute) |>
        #     filter(!NCBI_ID %in% dat$NCBI_ID) |>
        #     bind_rows(dat)
        # 
        # new_taxids <- new_dat |> 
        #     pull(taxid) |> 
        #     unique() |>
        #     {\(y) y[!is.na(y)]}()
        # per <- mean(tip_data$taxid %in% new_taxids) * 100
        # if (per < 1) {
        #     return(per)
        # }
        
        per_annotated_tips <- round(mean(unique(dat$NCBI_ID) %in% tree$tip.label) * 100)
        
        tip_data_annotated <- left_join(
            x = tip_data,
            y = select(dat, taxid, Attribute, Score),
            by = 'taxid', # joining by taxid helps to join cases like 561 and g__561
            relationship = 'many-to-many'
        )
        
        annotated_tips <- tip_data_annotated |>
            select(tip_label, Attribute, Score) |>
            filter(!is.na(Attribute)) |>
            pivot_wider(
                names_from = 'Attribute', values_from = 'Score', values_fill = 0
            ) |>
            tibble::column_to_rownames(var = 'tip_label') |>
            as.matrix()
        
        annotated_gn_tips <- grep('g__', rownames(annotated_tips), value = TRUE)
        remove_gn_tips <- gn_tips[!gn_tips %in% annotated_gn_tips]
        
        tip_data <- tip_data |> 
            filter(!tip_label %in% remove_gn_tips)
        
        if (attribute_type %in% c('binary', 'multistate-union')) {
            no_annotated_tips <- tip_data |>
                filter(!tip_label %in% rownames(annotated_tips)) |>
                select(tip_label) |>
                mutate(
                    Attribute = factor(
                        attribute_nms[[1]], levels = attribute_nms
                    )
                ) |>
                complete(tip_label, Attribute) |>
                mutate(Score = ifelse(grepl('--FALSE$', Attribute), 1, 0)) |>
                pivot_wider(names_from = 'Attribute', values_from = 'Score') |>
                tibble::column_to_rownames(var = 'tip_label') |>
                as.matrix() |>
                {\(y) y[,colnames(annotated_tips)]}()
            
        } else if (attribute_type == 'multistate-intersection') {
            no_annotated_tip_names <- tip_data |>
                filter(!tip_label %in% rownames(annotated_tips)) |>
                pull(tip_label)
            fill_value <- 1 / length(attribute_nms)
            vct <- rep(
                fill_value,
                length(no_annotated_tip_names) * length(attribute_nms)
            )
            no_annotated_tips <- matrix(
                data = vct,
                nrow = length(no_annotated_tip_names),
                ncol = length(attribute_nms),
                dimnames = list(no_annotated_tip_names, attribute_nms)
            )
        }
       
        tree <- ape::drop.tip(phy = tree, tip = remove_gn_tips) ## After removing this, node numbers change
        input_matrix <- rbind(annotated_tips, no_annotated_tips)
        input_matrix <- input_matrix[tree$tip.label,]
        
        models <- c(ER = 'ER', ARD = 'ARD', SYM = 'SYM')
        fittedModels <- map(models, ~ {
                message('Fitting ', .x)
                start_time <- Sys.time()
                message(start_time)
                if (.x == 'ER') {
                        r <- fitMk(
                                tree = tree, x = input_matrix, model = .x,
                                pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE
                        )
                } else {
                        r <- fitMk.parallel(
                                tree = tree, x = input_matrix, model = .x,
                                pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE,
                                ncores = n_cores
                        )
                }
                end_time <- Sys.time()
                elapsed_time <- difftime(end_time, start_time)
                message('Took ', elapsed_time)
                return(r)
        })
        
        bestModel <- names(sort(aic.w(aic = map_dbl(fittedModels, AIC)), decreasing = TRUE))[1]
        asr <- ancr(fittedModels[[bestModel]], tips = TRUE)
        res <- asr$ace |> 
            as.data.frame() |> 
            mutate(label = c(tree$tip.label, tree$node.label)) |> # Use this instead of rownames_to_colmn because labels are lost in res$ace
            filter(!grepl('^n\\d+$', label)) |> 
            filter(label != 'NA') |> 
            filter(!label %in% rownames(annotated_tips)) # this also remove all remaining genus tips # this correspond to original annotations so the will be in the sources
        res_tips_df <- res |>
            as.data.frame() |>
            filter(!grepl('^\\d+(\\+\\d+)*$', label)) |> 
            rename(tip_label = label)
        res_nodes_df <- res |>
            as.data.frame() |>
            filter(grepl('^\\d+(\\+\\d+)*$', label)) |> 
            rename(node_label = label)
        
        ## Get annotations for tips and nodes
        new_tips_data <- tip_data |>
            filter(tip_label %in% unique(res_tips_df$tip_label)) |>
            select(tip_label, taxid, Taxon_name, Rank) |>
            group_by(taxid) |>
            slice_max(order_by = tip_label, n = 1) |>
            ungroup() |>
            left_join(res_tips_df, by = 'tip_label') |>
            mutate(
                NCBI_ID = case_when(
                    Rank == 'species' ~ paste0('s__', taxid),
                    Rank == 'strain' ~ paste0('t__', taxid)
                )
            ) |>
            filter(Rank %in% c('species', 'strain')) |>
            mutate(Evidence = 'asr') |>
            relocate(NCBI_ID, taxid, Taxon_name, Rank, Evidence) |>
            pivot_longer(
                cols = 7:last_col(), names_to = 'Attribute', values_to = 'Score'
            ) |>
            mutate(
                Attribute_source = NA,
                Confidence_in_curation = NA,
                Attribute_group = attribute_group,
                Attribute_type = attribute_type,
                Frequency = case_when(
                    Score == 1 ~ 'always',
                    Score > 0.9 ~ 'usually',
                    Score >= 0.5 ~ 'sometimes',
                    Score > 0 & Score < 0.5 ~ 'rarely',
                    Score == 0 ~ 'never'
                )
            ) |>
            select(-tip_label)
        
        new_nodes_data <- node_data |>
            filter(node_label %in% unique(res_nodes_df$node_label)) |>
            select(node_label, taxid, Taxon_name, Rank) |>
            left_join(res_nodes_df, by = 'node_label') |>
            mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
            mutate(
                NCBI_ID = case_when(
                    Rank == 'kingdom' ~ paste0('k__', taxid),
                    Rank == 'phylum' ~ paste0('p__', taxid),
                    Rank == 'class' ~ paste0('c__', taxid),
                    Rank == 'order' ~ paste0('o__', taxid),
                    Rank == 'family' ~ paste0('f__', taxid),
                    Rank == 'genus' ~ paste0('g__', taxid),
                    Rank == 'species' ~ paste0('s__', taxid),
                    Rank == 'strain' ~ paste0('t__', taxid)
                )
            ) |>
            filter(
                Rank %in% c(
                    'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                    'species', 'strain'
                )
            ) |>
            mutate(Evidence = 'asr') |>
            relocate(NCBI_ID, taxid, Taxon_name, Rank, Evidence) |>
            pivot_longer(
                cols = 7:last_col(), names_to = 'Attribute', values_to = 'Score'
            ) |>
            mutate(
                Attribute_source = NA,
                Confidence_in_curation = NA,
                Attribute_group = attribute_group,
                Attribute_type = attribute_type,
                Frequency = case_when(
                    Score == 1 ~ 'always',
                    Score > 0.9 ~ 'usually',
                    Score >= 0.5 ~ 'sometimes',
                    Score > 0 & Score < 0.5 ~ 'rarely',
                    Score == 0 ~ 'never'
                )
            ) |>
            select(-node_label) |>
            filter(!NCBI_ID %in% dat$NCBI_ID)
        
        new_taxa_for_ncbi_tree <- bind_rows(new_tips_data, new_nodes_data) |>
            relocate(NCBI_ID, Rank, Attribute, Score, Evidence)
        
        new_taxa_for_ncbi_tree_list <- split(
            new_taxa_for_ncbi_tree, factor(new_taxa_for_ncbi_tree$NCBI_ID)
        )
        
        ## Add annotations from sources
        ncbi_tree$Do(function(node) {
            if (node$name %in% names(node_list))
                node$attribute_tbl <- node_list[[node$name]]
        })
       
        ## Add annotatiosn from ASR 
        ncbi_tree$Do(function(node) {
            cond1 <- node$name %in% names(new_taxa_for_ncbi_tree_list)
            cond2 <- is.null(node$attribute_tbl) || all(is.na(node$attribute_tbl))
            if (cond1 && cond2) {
                node$attribute_tbl <- new_taxa_for_ncbi_tree_list[[node$name]]
            }
        })
        
        ncbi_tree$Do(function(nd) inh1(node = nd, adjF = 0.1, evidence_label = 'inh2'), traversal = 'pre-order')
        
        result <- ncbi_tree$Get(
            attribute = 'attribute_tbl', simplify = FALSE,
            filterFun = function(node) {
                node$name != 'ArcBac' && !is.null(node$attribute_tbl)
                # grepl('^gst__', node$name)
                
            }
        ) |>
            bind_rows() |>
            discard(~ all(is.na(.x)))
        
        min_thr <- 1 / length(unique(dat$Attribute))
        
        final_result <- bind_rows(
            dat, # contains source and tax
            new_taxa_for_ncbi_tree, # only contains asr (not in new_dat)
            filter(result, Evidence == 'inh2')
        ) |>
            filter(Score > min_thr)
        ncbi_tree$Do(cleanNode)

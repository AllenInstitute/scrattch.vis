#' Combine data for plotting from the included Mouse V1 data
#' Current scrattch no longer contains this dataset.
# get_internal_data <- function(genes,group_by,clusters) {
#   data <- scrattch::v1_data
#   all.anno <- scrattch::v1_anno
#   
#   cluster_order <- data.frame(clusters=clusters) %>%
#     mutate(cluster_x=1:n())
#   
#   data <- data %>%
#     filter(gene %in% genes)
#   
#   genes <- sub("-",".",genes)
#   
#   # Reformat the retrieved data
#   row.names(data) <- sub("-",".",data[,1])
#   data <- data %>% 
#     select(-1) %>% 
#     t() %>% 
#     as.data.frame()
#   
#   data <- data %>%
#     mutate(sample_id=row.names(data)) %>%
#     select(one_of(c("sample_id",genes)))
#   
#   genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
#   names(data)[grepl("^[0-9]",names(data))] <- paste0("X",names(data)[grepl("^[0-9]",names(data))])
#   
#   # Filter and order the rows
#   data <- left_join(data,all.anno,by="sample_id") %>%
#     rename_("plot_id" = paste0(group_by,"_id"),
#             "plot_label" = paste0(group_by,"_label"),
#             "plot_color" = paste0(group_by,"_color")) %>%
#     filter(plot_id %in% clusters) %>%
#     left_join(cluster_order,by=c("plot_id"="clusters")) %>%
#     arrange(cluster_x) %>%
#     mutate(xpos=1:n()) %>%
#     select(-plot_id) %>%
#     rename_("plot_id"="cluster_x")
#   
#   return(data)
# }


#' Get data from a SQLite3 database for scrattch plotting
#' SQLite3 is retired from use. Won't be used in scrattch release (at least initially).
#' 
# get_db_data <- function(data_source,genes,group_by,group_ids) {
#   library(dplyr)
#   library(DBI)
#   library(RSQLite)
#   
#   if(file.exists(data_source)) {
#     con <- dbConnect(RSQLite::SQLite(),data_source)
#     
#     # Get annotations
#     get <- "SELECT * FROM anno;"
#     res <- dbSendQuery(con,get)
#     all.anno <- dbFetch(res,n=-1)
#     dbClearResult(res)
#     
#     # Rename based on group_by, and select annotations in the group_ids
#     anno <- all.anno %>%
#       rename_("plot_id" = paste0(group_by,"_id"),
#               "plot_label" = paste0(group_by,"_label"),
#               "plot_color" = paste0(group_by,"_color")) %>%
#       filter(plot_id %in% group_ids)
#     
#     # Get data for genes
#     #  getgenes <- chr_to_sql(toupper(genes))
#     #  get <- paste("SELECT * FROM data WHERE upper(gene) IN ",getgenes,sep="")
#     getgenes <- chr_to_sql(genes)
#     get <- paste("SELECT * FROM data WHERE gene IN ",getgenes,sep="")
#     res <- dbSendQuery(con,get)
#     data <- dbFetch(res,n=-1)
#     dbClearResult(res)
#     
#     dbDisconnect(con)
#     
#     genes <- sub("-",".",genes)
#     
#     # transpose genes table for joining to annotations
#     row.names(data) <- sub("-",".",data[,1])
#     data <- data %>% 
#       select(-1) %>% 
#       t() %>% 
#       as.data.frame()
#     names(data) <- toupper(names(data))
#     
#     # Filter data to only samples in the filtered annotations
#     data <- data %>%
#       mutate(sample_id=row.names(data)) %>%
#       filter(sample_id %in% anno$sample_id) %>%
#       select(one_of(c("sample_id",toupper(genes))))
#     names(data) <- c("sample_id",genes)
#     
#     # join data and annotations
#     data <- left_join(data,anno,by="sample_id")
#     
#     # Set cluster order based on the input group_ids
#     cluster_order <- data.frame(plot_id = group_ids) %>%
#       filter(group_ids %in% unique(data$plot_id)) %>%
#       mutate(cluster_x=1:n())
#     
#     # sort data by input cluster order
#     data <- data %>%
#       left_join(cluster_order,by="plot_id") %>%
#       arrange(cluster_x) %>%
#       mutate(xpos=1:n()) %>%
#       select(-plot_id) %>%
#       rename_("plot_id" = "cluster_x")
#     
#     return(data)
#   } else {
#     cat("Database file not found!\n")
#   }
# }


#' Format data provided in list format for scrattch plots
#' 
#' Currently only compatible with data from feather_to_list()
# get_list_data <- function(data_list, genes, group_by, group_ids) {
#   
#   library(dplyr)
#   
#   # Read annotations and convert factors
#   anno <- data_list$anno %>%
#     mutate_if(is.factor, as.character)
#   
#   # If an _id column was a factor, it's now a character. Convert to numeric for sorting.
#   id_cols <- names(anno)[grepl("_id$", names(anno)) & names(anno) != "sample_id"]
#   anno[id_cols] <- lapply(anno[id_cols], as.numeric)
#   
#   # Check the provided genes against the column names in data_file
#   data_names <- names(data_list$data)
#   
#   if(sum(genes %in% data_names) != length(genes)) {
#     # Report if names don't match after ignorning case
#     not_found <- genes[!toupper(genes) %in% toupper(data_names)]
#     
#     warning(paste(paste0(not_found, collapse = ", "), "not found in data table!"))
#     
#     # Update genes to use names as formatted in data
#     genes <- data_names[toupper(data_names) %in% toupper(genes)]
#   }
#   
#   gene_data <- data_list$data[,c("sample_id",genes)]
#   
#   # Change - to . in column names and genes
#   colnames(gene_data) <- gsub("-",".",colnames(gene.data))
#   genes <- gsub("-",".",genes)
#   
#   # rename the _id, _label, and _color for the group_by values for use in plotting
#   all_anno <- anno %>%
#     rename_("plot_id" = paste0(group_by,"_id"),
#             "plot_label" = paste0(group_by,"_label"),
#             "plot_color" = paste0(group_by,"_color"))
#   
#   # use the group_ids to retain the order provided by the group_ids argument
#   cluster_order <- data.frame(group_ids = group_ids) %>%
#     mutate(cluster_x = 1:n())
#   
#   # Filter and order the rows
#   data <- left_join(all_anno, gene.data, by = "sample_id") %>%
#     filter(plot_id %in% group_ids) %>%
#     left_join(cluster_order, by = c("plot_id" = "group_ids")) %>%
#     arrange(cluster_x) %>%
#     mutate(xpos = 1:n()) %>%
#     select(-plot_id) %>%
#     rename_("plot_id" = "cluster_x")
#   
#   return(data)
# }


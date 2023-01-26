#' Inject LSI object from ArchR into Seurat Object
#' Adds the lsi to the input Seurat object as a DimReducObject
#' under the specified assay
#' @param so The Seruat Object
#' @param proj The ArchR object
#' @param assay The assay to add the LSI 
#' @param reduction.name The name of the new LSI reduction
#' @param reduction.key The key for the reduction dimension names
#' @return The original seurat object with injected LSI
inject_lsi_seurat <- function(so, proj, assay = "ATAC", reduction.name = "lsi_atac", reduction.key="lsit_"){
    lsi <- getReducedDims(proj)
    atac_meta <- ArchR::getCellColData(proj)
    atac_meta <- atac_meta[, c("barcodes", setdiff(colnames(atac_meta), colnames(so@meta.data)))]
    
    # clean up lsi rownames to just the cell barcode part
    rownames(lsi) <- gsub(".*#", "" , rownames(lsi))
    rownames(atac_meta) <- gsub(".*#", "" , rownames(atac_meta))
    
    if(!all(colnames(so) %in% rownames(lsi))){
        common_cells <- intersect(colnames(so), rownames(lsi))
        message(sprintf("More cells in seurat object [%s] than reduction matrix [%s]. Filtering to [%s] common cells", 
                        ncol(so), 
                        length(rownames(lsi)), 
                        length(common_cells)))
        so <- so[, common_cells]
    }

    # create a new assay. note that 
    if(!assay %in% Seurat::Assays(so)){
        message(sprintf("Creating new assay '%s' to store lsi", assay))
        dummy_tile <- matrix(data = 0, nrow = 2, ncol = ncol(so))
        colnames(dummy_tile) <- rownames(so@meta.data)
        rownames(dummy_tile) <- LETTERS[1:2]
        so[[assay]] <- Seurat::CreateAssayObject(counts = dummy_tile)
    }

    # Add lsi as dim reduction
    message("Adding LSI to seurat Object")
    so[[reduction.name]] <- CreateDimReducObject(embeddings = lsi, key = reduction.key, assay = assay)
    
    # Add metadata
    message("Adding ATAC metadata to seurat Object")
    so <- AddMetaData(so, as.data.frame(atac_meta))
    
    return(so)
}




get_umap_df <- function(so, reduction_name=c('umap','pca'), keep_all = FALSE, reduction_index = c(1:2)){
    meta_umap <- so@meta.data
    if(!keep_all){
        if(is.list(reduction_index)){
            assertthat::assert_that(all(reduction_name) %in% names(reduction_index))
        }
    }
    for(rname in reduction_name){
        reduction_df <- so@reductions[[rname]]@cell.embeddings
        if(!keep_all){
            if(is.list(reduction_index)){
                ikeep <- reduction_index[[rname]]
            } else { 
                ikeep <- reduction_index
            }
            reduction_df <- reduction_df[,ikeep]
        }
        assertthat::assert_that(all(rownames(meta_umap) == rownames(reduction_df)))
        meta_umap <- cbind(meta_umap, reduction_df)
    }
    return(meta_umap)
}

create_seurat <- function(h5_file_list, keep_assays = c("ADT")){
    
    so_mat_list <- lapply(h5_file_list, read_h5_dgCMatrix,  target = "ADT", feature_names = "id", sample_names = "barcodes", index1 = TRUE)
    mat <- tryCatch(dplyr::bind_cols(so_mat_list), 
                    error = function(e){
                        print("")
                    }
    )
    
    adt <- read_h5_dgCMatrix(h5_file, "ADT", feature_names = "id")
    colnames(adt) <- colnames(so)
    so[["ADT"]] <- CreateAssayObject(adt)
    
    seurat_labels <- read.csv(seurat_label_file)
    PedVsSnr_labels <- read.csv(PedVsSnr_label_file)
    
    cell_ids <- rownames(so@meta.data)
    so@meta.data <- left_join(so@meta.data, seurat_labels)
    so@meta.data <- left_join(so@meta.data, PedVsSnr_labels)
    rownames(so@meta.data) <- cell_ids
    
    so
}

# Keep data and general metadata only, remove dataset-specific analysis results
strip_seurat <- function(so, keep_assays = c("RNA", "ADT"), rm_pattern =c("_res.[0-9.]+$", "[.]weight$")){
    DefaultAssay(so) <- keep_assays[1]
    so <- DietSeurat(so, assays=keep_assays)
    
    for (pat in rm_pattern){
        ipat <- grep(pat, names(so@meta.data))
        if(length(ipat)>0){
            message(sprintf("Removing column %s from metadata\n", names(so@meta.data)[ipat]))
            so@meta.data[,ipat] <- NULL
        }
    }
    return(so)
}
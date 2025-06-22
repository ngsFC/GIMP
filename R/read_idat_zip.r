#' Check Available minfi Functions
#' @export
check_minfi_functions <- function() {
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("minfi package not installed")
  }
  
  library(minfi)
  
  cat("minfi Package Information\n")
  cat("========================\n")
  cat("Version:", as.character(packageVersion("minfi")), "\n\n")
  
  # Get all functions
  all_functions <- ls(envir = asNamespace("minfi"))
  
  # Look for detection-related functions
  cat("Detection-related functions:\n")
  detection_funcs <- grep("etection|etP", all_functions, value = TRUE, ignore.case = TRUE)
  if (length(detection_funcs) > 0) {
    for (func in detection_funcs) {
      cat("  -", func, "\n")
    }
  } else {
    cat("  No detection functions found\n")
  }
  
  # Look for preprocessing functions
  cat("\nPreprocessing functions:\n")
  preprocess_funcs <- grep("preprocess", all_functions, value = TRUE, ignore.case = TRUE)
  for (func in head(preprocess_funcs, 10)) {
    cat("  -", func, "\n")
  }
  
  # Look for get functions
  cat("\nGetter functions:\n")
  get_funcs <- grep("^get", all_functions, value = TRUE, ignore.case = TRUE)
  for (func in head(get_funcs, 10)) {
    cat("  -", func, "\n")
  }
  
  # Test specific functions we need
  cat("\nTesting required functions:\n")
  
  required_funcs <- c("read.metharray.exp", "getBeta", "getDetectionP", "detectionP", 
                      "preprocessQuantile", "preprocessSWAN", "preprocessFunnorm", "preprocessNoob")
  
  for (func in required_funcs) {
    if (exists(func, envir = asNamespace("minfi"))) {
      cat("  ✅", func, "- available\n")
    } else {
      cat("  ❌", func, "- NOT available\n")
    }
  }
  
  return(all_functions)
}

#' Alternative Detection P-value Calculation
#' @param rgSet RGChannelSet object
#' @export
calculate_detection_pvalues <- function(rgSet) {
  
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("minfi package required")
  }
  
  library(minfi)
  
  # Try different approaches
  tryCatch({
    # Method 1: Standard function
    if (exists("detectionP", envir = asNamespace("minfi"))) {
      return(detectionP(rgSet))
    }
  }, error = function(e) {
    cat("Method 1 failed:", e$message, "\n")
  })
  
  tryCatch({
    # Method 2: Alternative name
    if (exists("getDetectionP", envir = asNamespace("minfi"))) {
      return(getDetectionP(rgSet))
    }
  }, error = function(e) {
    cat("Method 2 failed:", e$message, "\n")
  })
  
  tryCatch({
    # Method 3: Manual calculation
    cat("Attempting manual detection p-value calculation...\n")
    
    # Get control probes
    if (exists("getProbeInfo", envir = asNamespace("minfi"))) {
      ctrl_probes <- getProbeInfo(rgSet, type = "Control")
      
      # Simple approach: use background correction as proxy for detection
      # This is a simplified version and may not be as accurate
      green <- getGreen(rgSet)
      red <- getRed(rgSet)
      
      # Calculate background levels (very simplified)
      green_bg <- apply(green, 2, function(x) quantile(x, 0.05, na.rm = TRUE))
      red_bg <- apply(red, 2, function(x) quantile(x, 0.05, na.rm = TRUE))
      
      # Simple detection: signal significantly above background
      detection_matrix <- matrix(0.01, nrow = nrow(green), ncol = ncol(green))
      rownames(detection_matrix) <- rownames(green)
      colnames(detection_matrix) <- colnames(green)
      
      # Mark low-intensity probes as failed
      for (i in 1:ncol(green)) {
        low_green <- green[, i] < (green_bg[i] * 2)
        low_red <- red[, i] < (red_bg[i] * 2)
        detection_matrix[low_green | low_red, i] <- 0.1  # High p-value = failed
      }
      
      cat("Manual detection p-values calculated (simplified approach)\n")
      return(detection_matrix)
    }
  }, error = function(e) {
    cat("Method 3 failed:", e$message, "\n")
  })
  
  # If all methods fail, return NULL
  cat("All detection p-value methods failed. Proceeding without detection QC.\n")
  return(NULL)
}

#' Read IDAT Files from ZIP Archive (Robust version with proper indexing)
#'
#' @param zip_file Path to ZIP file containing IDAT files and sample sheet
#' @param sample_sheet_name Name of the sample sheet file in the ZIP (default: "samplesheet.csv")
#' @param array_type Array type for annotation ("450k", "EPIC", "EPICv2")
#' @param temp_dir Temporary directory for extraction (default: creates temporary directory)
#' @param normalize_method Normalization method for minfi ("quantile", "SWAN", "funnorm", "noob")
#' @param detection_pval P-value threshold for detection (default: 0.01)
#' @param remove_failed_samples Remove samples with >10% failed probes (default: TRUE)
#' @return A list containing beta_matrix, sample_info, qc_metrics, failed_samples
#' @export

read_idat_zip <- function(zip_file, 
                          sample_sheet_name = "samplesheet.csv",
                          array_type = c("EPIC", "450k", "EPICv2"),
                          temp_dir = NULL,
                          normalize_method = c("quantile", "SWAN", "funnorm", "noob"),
                          detection_pval = 0.01,
                          remove_failed_samples = TRUE) {
  
  # Check and load required packages
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("minfi package is required. Install with: BiocManager::install('minfi')")
  }
  
  # Load minfi
  library(minfi)
  
  # Validate inputs
  array_type <- match.arg(array_type)
  normalize_method <- match.arg(normalize_method)
  
  if (!file.exists(zip_file)) {
    stop("ZIP file not found: ", zip_file)
  }
  
  # Create temporary directory if not provided
  if (is.null(temp_dir)) {
    temp_dir <- tempdir()
    cleanup_temp <- TRUE
  } else {
    cleanup_temp <- FALSE
  }
  
  extract_dir <- file.path(temp_dir, paste0("idat_extraction_", Sys.getpid()))
  
  # Clean up any existing extraction directory
  if (dir.exists(extract_dir)) {
    unlink(extract_dir, recursive = TRUE)
  }
  
  # Initialize variables
  sample_sheet <- NULL
  rgSet <- NULL
  original_sample_count <- 0
  
  tryCatch({
    cat("=== IDAT PROCESSING STARTED ===\n")
    cat("Extracting ZIP file to:", extract_dir, "\n")
    
    # Extract ZIP file
    extracted_files <- unzip(zip_file, exdir = extract_dir)
    cat("Extracted", length(extracted_files), "files\n")
    
    # Find sample sheet
    sample_sheet_path <- file.path(extract_dir, sample_sheet_name)
    
    # Try alternative sample sheet names if not found
    if (!file.exists(sample_sheet_path)) {
      possible_sheets <- c("SampleSheet.csv", "sample_sheet.csv", "samples.csv", "metadata.csv")
      found_sheet <- FALSE
      
      for (sheet in possible_sheets) {
        alt_path <- file.path(extract_dir, sheet)
        if (file.exists(alt_path)) {
          sample_sheet_path <- alt_path
          cat("Using sample sheet:", sheet, "\n")
          found_sheet <- TRUE
          break
        }
      }
      
      if (!found_sheet) {
        # Look for CSV files in extracted directory
        csv_files <- list.files(extract_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
        if (length(csv_files) > 0) {
          sample_sheet_path <- csv_files[1]
          cat("Using first CSV file found as sample sheet:", basename(sample_sheet_path), "\n")
        } else {
          stop("No sample sheet found. Expected: ", sample_sheet_name, 
               " or one of: SampleSheet.csv, sample_sheet.csv, samples.csv")
        }
      }
    }
    
    # Read sample sheet
    cat("Reading sample sheet from:", sample_sheet_path, "\n")
    sample_sheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)
    original_sample_count <- nrow(sample_sheet)
    
    cat("Sample sheet dimensions:", nrow(sample_sheet), "rows,", ncol(sample_sheet), "columns\n")
    cat("Sample sheet columns:", paste(colnames(sample_sheet), collapse = ", "), "\n")
    
    # Validate sample sheet columns
    required_cols <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position")
    missing_cols <- setdiff(required_cols, colnames(sample_sheet))
    
    if (length(missing_cols) > 0) {
      # Try alternative column names
      alt_names <- list(
        Sample_Name = c("Sample_ID", "SampleID", "sample_name", "sample_id"),
        Sentrix_ID = c("Slide", "Array", "sentrix_id", "chip_id"),
        Sentrix_Position = c("Array_Position", "Position", "sentrix_position", "array_position")
      )
      
      for (col in missing_cols) {
        found <- FALSE
        for (alt in alt_names[[col]]) {
          if (alt %in% colnames(sample_sheet)) {
            colnames(sample_sheet)[colnames(sample_sheet) == alt] <- col
            found <- TRUE
            cat("Mapped column", alt, "to", col, "\n")
            break
          }
        }
        if (!found) {
          stop("Required column '", col, "' not found in sample sheet.\nAvailable columns: ", 
               paste(colnames(sample_sheet), collapse = ", "),
               "\nRequired columns: ", paste(required_cols, collapse = ", "))
        }
      }
    }
    
    # Clean up sample sheet data and ensure character types
    sample_sheet$Sample_Name <- as.character(sample_sheet$Sample_Name)
    sample_sheet$Sentrix_ID <- as.character(sample_sheet$Sentrix_ID)
    sample_sheet$Sentrix_Position <- as.character(sample_sheet$Sentrix_Position)
    
    # Remove any empty rows - be very careful with indexing here
    valid_rows <- !is.na(sample_sheet$Sample_Name) & 
      sample_sheet$Sample_Name != "" &
      !is.na(sample_sheet$Sentrix_ID) & 
      sample_sheet$Sentrix_ID != "" &
      !is.na(sample_sheet$Sentrix_Position) & 
      sample_sheet$Sentrix_Position != ""
    
    cat("Valid rows before filtering:", sum(valid_rows), "out of", length(valid_rows), "\n")
    
    sample_sheet <- sample_sheet[valid_rows, , drop = FALSE]
    
    cat("Valid samples in sheet after filtering:", nrow(sample_sheet), "\n")
    
    if (nrow(sample_sheet) == 0) {
      stop("No valid samples found in sample sheet after filtering")
    }
    
    # Get all IDAT files in the directory
    all_idat_files <- list.files(extract_dir, pattern = "\\.idat$", full.names = TRUE, recursive = TRUE)
    
    if (length(all_idat_files) == 0) {
      stop("No IDAT files found in the ZIP archive")
    }
    
    cat("Found", length(all_idat_files), "IDAT files in ZIP\n")
    
    # Print some example IDAT file names for debugging
    idat_basenames <- unique(gsub("_(Red|Grn)\\.idat$", "", basename(all_idat_files)))
    cat("Example IDAT identifiers (first 5):\n")
    cat(paste(head(idat_basenames, 5), collapse = "\n"), "\n")
    cat("Total unique IDAT basenames:", length(idat_basenames), "\n")
    
    # Create Basename column for minfi - try multiple strategies
    cat("Attempting to match IDAT files with sample sheet entries...\n")
    
    # Strategy 1: Standard format (Sentrix_ID_Sentrix_Position)
    sample_sheet$Basename <- file.path(dirname(all_idat_files[1]), 
                                       paste0(sample_sheet$Sentrix_ID, "_", sample_sheet$Sentrix_Position))
    
    # Check if files exist with standard format
    idat_files_red <- paste0(sample_sheet$Basename, "_Red.idat")
    idat_files_grn <- paste0(sample_sheet$Basename, "_Grn.idat")
    missing_files <- !(file.exists(idat_files_red) & file.exists(idat_files_grn))
    
    cat("Standard matching: ", sum(!missing_files), "out of", length(missing_files), "samples matched\n")
    
    if (any(missing_files)) {
      cat("Trying flexible matching for", sum(missing_files), "unmatched samples...\n")
      
      # Strategy 2: Flexible matching
      for (i in which(missing_files)) {
        # Try different patterns
        patterns_to_try <- c(
          paste0(sample_sheet$Sentrix_ID[i], "_", sample_sheet$Sentrix_Position[i]),
          paste0(sample_sheet$Sample_Name[i], "_", sample_sheet$Sentrix_ID[i], "_", sample_sheet$Sentrix_Position[i]),
          sample_sheet$Sample_Name[i],
          gsub("GSM", "", sample_sheet$Sample_Name[i])
        )
        
        found_match <- FALSE
        for (pattern in patterns_to_try) {
          matching_files <- grep(pattern, idat_basenames, value = TRUE)
          
          if (length(matching_files) == 1) {
            # Found unique match
            sample_sheet$Basename[i] <- file.path(dirname(all_idat_files[1]), matching_files[1])
            found_match <- TRUE
            cat("  Matched sample", sample_sheet$Sample_Name[i], "with pattern:", matching_files[1], "\n")
            break
          } else if (length(matching_files) > 1) {
            # Multiple matches - try to find exact match
            exact_match <- grep(paste0("^", pattern, "$"), idat_basenames, value = TRUE)
            if (length(exact_match) == 1) {
              sample_sheet$Basename[i] <- file.path(dirname(all_idat_files[1]), exact_match[1])
              found_match <- TRUE
              cat("  Matched sample", sample_sheet$Sample_Name[i], "with exact pattern:", exact_match[1], "\n")
              break
            }
          }
        }
        
        if (!found_match) {
          cat("  ⚠️  No match found for sample:", sample_sheet$Sample_Name[i], 
              "(", sample_sheet$Sentrix_ID[i], "_", sample_sheet$Sentrix_Position[i], ")\n")
        }
      }
    }
    
    # Final validation of IDAT files
    idat_files_red <- paste0(sample_sheet$Basename, "_Red.idat")
    idat_files_grn <- paste0(sample_sheet$Basename, "_Grn.idat")
    final_missing_files <- !(file.exists(idat_files_red) & file.exists(idat_files_grn))
    
    if (any(final_missing_files)) {
      missing_samples <- sample_sheet$Sample_Name[final_missing_files]
      
      cat("\n❌ FINAL IDAT FILE MATCHING FAILED\n")
      cat("==================================\n")
      cat("Could not find IDAT files for", length(missing_samples), "samples:\n")
      
      for (i in which(final_missing_files)[1:min(5, length(which(final_missing_files)))]) {
        cat("\nSample:", sample_sheet$Sample_Name[i], "\n")
        cat("  Expected pattern:", basename(sample_sheet$Basename[i]), "\n")
        cat("  Looking for Red:", basename(paste0(sample_sheet$Basename[i], "_Red.idat")), "\n")
        cat("  Looking for Grn:", basename(paste0(sample_sheet$Basename[i], "_Grn.idat")), "\n")
      }
      
      if (length(missing_samples) > 5) {
        cat("... and", length(missing_samples) - 5, "more samples\n")
      }
      
      cat("\nAvailable IDAT patterns (first 10):\n")
      cat(paste(head(idat_basenames, 10), collapse = "\n"), "\n")
      
      stop("Missing IDAT files for samples: ", paste(head(missing_samples, 10), collapse = ", "),
           if(length(missing_samples) > 10) paste(" ... and", length(missing_samples) - 10, "more") else "")
    }
    
    # Filter sample sheet to only include samples with IDAT files
    valid_sample_indices <- !final_missing_files
    sample_sheet_filtered <- sample_sheet[valid_sample_indices, , drop = FALSE]
    
    cat("Processing", nrow(sample_sheet_filtered), "samples with valid IDAT files\n")
    
    if (nrow(sample_sheet_filtered) == 0) {
      stop("No samples with valid IDAT files found")
    }
    
    # Read methylation data with minfi
    cat("Reading methylation data with minfi...\n")
    cat("Sample sheet for minfi has", nrow(sample_sheet_filtered), "rows\n")
    
    # Ensure rownames are set properly for minfi
    rownames(sample_sheet_filtered) <- sample_sheet_filtered$Sample_Name
    
    rgSet <- read.metharray.exp(targets = sample_sheet_filtered, verbose = TRUE)
    
    cat("RGChannelSet dimensions:", dim(rgSet), "\n")
    cat("RGChannelSet sample names:", paste(head(colnames(rgSet)), collapse = ", "), "\n")
    
    # Quality control with proper dimension handling
    cat("Performing quality control...\n")
    
    # Initialize QC variables
    detP <- NULL
    failed_per_sample <- rep(0, ncol(rgSet))
    names(failed_per_sample) <- colnames(rgSet)
    
    # Try to get detection p-values with various methods
    tryCatch({
      # Try different detection functions
      if (exists("detectionP", envir = asNamespace("minfi"))) {
        detP <- detectionP(rgSet)
        cat("  Using detectionP function\n")
      } else {
        cat("  detectionP function not available, using simplified QC\n")
        # Use signal intensity as proxy for quality
        if (exists("getGreen", envir = asNamespace("minfi")) && exists("getRed", envir = asNamespace("minfi"))) {
          green_intensities <- getGreen(rgSet)
          red_intensities <- getRed(rgSet)
          
          # Calculate median intensities per sample
          green_medians <- apply(green_intensities, 2, median, na.rm = TRUE)
          red_medians <- apply(red_intensities, 2, median, na.rm = TRUE)
          
          # Flag samples with very low intensity
          low_intensity_threshold <- 1000
          low_intensity_samples <- (green_medians < low_intensity_threshold) | (red_medians < low_intensity_threshold)
          failed_per_sample[low_intensity_samples] <- 0.15
          
          cat("  Identified", sum(low_intensity_samples), "samples with low signal intensity\n")
        }
      }
    }, error = function(e) {
      cat("  Warning: Could not perform detailed QC:", e$message, "\n")
    })
    
    # Calculate failed probe proportions if we have detection p-values
    if (!is.null(detP)) {
      cat("  Calculating failed probe proportions...\n")
      cat("  Detection p-value matrix dimensions:", dim(detP), "\n")
      cat("  RGChannelSet dimensions:", dim(rgSet), "\n")
      
      # Ensure dimensions match
      if (ncol(detP) == ncol(rgSet)) {
        failed_per_sample <- colMeans(detP > detection_pval, na.rm = TRUE)
        names(failed_per_sample) <- colnames(rgSet)
      } else {
        cat("  Warning: Dimension mismatch between detP and rgSet\n")
      }
    }
    
    # Identify failed samples
    failed_samples <- names(failed_per_sample)[failed_per_sample > 0.1]
    
    if (length(failed_samples) > 0) {
      cat("Samples with potential quality issues:", length(failed_samples), "\n")
      cat("Failed samples:", paste(head(failed_samples, 5), collapse = ", "), 
          if(length(failed_samples) > 5) "..." else "", "\n")
      
      if (remove_failed_samples && length(failed_samples) < ncol(rgSet)) {
        cat("Removing", length(failed_samples), "potentially problematic samples\n")
        
        # Get indices of samples to keep
        keep_sample_names <- setdiff(colnames(rgSet), failed_samples)
        keep_indices <- colnames(rgSet) %in% keep_sample_names
        
        cat("Keeping", sum(keep_indices), "samples out of", ncol(rgSet), "\n")
        
        # Filter rgSet
        rgSet <- rgSet[, keep_indices]
        
        # Filter detection p-values if available
        if (!is.null(detP)) {
          detP <- detP[, keep_indices, drop = FALSE]
        }
        
        # Filter sample sheet to match
        sample_sheet_filtered <- sample_sheet_filtered[sample_sheet_filtered$Sample_Name %in% keep_sample_names, , drop = FALSE]
        
        cat("After QC filtering: ", nrow(sample_sheet_filtered), "samples remaining\n")
      }
    }
    
    # Normalization
    cat("Normalizing data using", normalize_method, "method...\n")
    
    mSet <- switch(normalize_method,
                   "quantile" = preprocessQuantile(rgSet),
                   "SWAN" = preprocessSWAN(rgSet),
                   "funnorm" = preprocessFunnorm(rgSet),
                   "noob" = preprocessNoob(rgSet)
    )
    
    cat("Normalized data dimensions:", dim(mSet), "\n")
    
    # Get beta values
    cat("Extracting beta values...\n")
    beta_matrix <- getBeta(mSet)
    
    cat("Beta matrix dimensions:", dim(beta_matrix), "\n")
    
    # Remove probes with high detection p-values if available
    if (!is.null(detP)) {
      cat("Filtering probes based on detection p-values...\n")
      
      # Ensure dimensions match before filtering
      if (nrow(detP) == nrow(beta_matrix) && ncol(detP) == ncol(beta_matrix)) {
        failed_probes <- rowMeans(detP > detection_pval, na.rm = TRUE) > 0.1
        cat("Removing", sum(failed_probes), "probes with high detection p-values\n")
        beta_matrix <- beta_matrix[!failed_probes, , drop = FALSE]
      } else {
        cat("Dimension mismatch - skipping probe filtering\n")
      }
    } else {
      # Alternative probe filtering
      cat("Filtering probes based on variance...\n")
      probe_vars <- apply(beta_matrix, 1, var, na.rm = TRUE)
      low_var_probes <- is.na(probe_vars) | probe_vars < 0.001
      cat("Removing", sum(low_var_probes), "probes with low variance\n")
      beta_matrix <- beta_matrix[!low_var_probes, , drop = FALSE]
    }
    
    cat("Final beta matrix dimensions:", dim(beta_matrix), "\n")
    
    # Ensure column names match sample names
    colnames(beta_matrix) <- sample_sheet_filtered$Sample_Name
    
    # Prepare sample information for GIMP
    sample_info_gimp <- sample_sheet_filtered$Sample_Name
    names(sample_info_gimp) <- sample_sheet_filtered$Sample_Name
    
    # Add group information if available
    if ("Sample_Group" %in% colnames(sample_sheet_filtered)) {
      sample_info_gimp <- sample_sheet_filtered$Sample_Group
      names(sample_info_gimp) <- sample_sheet_filtered$Sample_Name
    } else if ("Group" %in% colnames(sample_sheet_filtered)) {
      sample_info_gimp <- sample_sheet_filtered$Group
      names(sample_info_gimp) <- sample_sheet_filtered$Sample_Name
    }
    
    # Quality control metrics
    qc_metrics <- list(
      original_samples = original_sample_count,
      total_samples = nrow(sample_sheet_filtered),
      total_probes = nrow(beta_matrix),
      failed_samples = failed_samples,
      failed_sample_proportion = failed_per_sample[sample_sheet_filtered$Sample_Name],
      array_type = array_type,
      normalization = normalize_method,
      detection_qc_available = !is.null(detP)
    )
    
    cat("✅ Successfully processed", ncol(beta_matrix), "samples and", nrow(beta_matrix), "probes\n")
    cat("=== IDAT PROCESSING COMPLETED ===\n")
    
    # Cleanup
    if (cleanup_temp && dir.exists(extract_dir)) {
      unlink(extract_dir, recursive = TRUE)
    }
    
    return(list(
      beta_matrix = beta_matrix,
      sample_info = sample_info_gimp,
      sample_sheet = sample_sheet_filtered,
      qc_metrics = qc_metrics,
      failed_samples = failed_samples
    ))
    
  }, error = function(e) {
    # Cleanup on error
    if (cleanup_temp && dir.exists(extract_dir)) {
      unlink(extract_dir, recursive = TRUE)
    }
    
    # Detailed error information
    cat("\n=== ERROR DEBUGGING INFO ===\n")
    cat("Error message:", e$message, "\n")
    cat("Original sample count:", original_sample_count, "\n")
    
    if (!is.null(sample_sheet)) {
      cat("Sample sheet loaded: YES (", nrow(sample_sheet), "rows )\n")
      cat("Sample sheet columns:", paste(colnames(sample_sheet), collapse = ", "), "\n")
    } else {
      cat("Sample sheet loaded: NO\n")
    }
    
    if (exists("rgSet") && !is.null(rgSet)) {
      cat("RGChannelSet created: YES (", paste(dim(rgSet), collapse = "x"), ")\n")
    } else {
      cat("RGChannelSet created: NO\n")
    }
    
    if (exists("all_idat_files")) {
      cat("IDAT files found:", length(all_idat_files), "\n")
    } else {
      cat("IDAT files found: NO\n")
    }
    
    cat("==========================\n")
    
    stop("Error processing IDAT files: ", e$message)
  })
}
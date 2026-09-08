#!/usr/bin/env Rscript 

suppressPackageStartupMessages({
  library(fs)
  library(stringr)
  library(purrr)
  library(arrow)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir        <- if (length(args) >= 1) fs::path_expand(args[[1]]) else getwd()
db_hint_arg     <- if (length(args) >= 2) fs::path_expand(args[[2]]) else fs::path(base_dir, "data", "reference-fasta")
fns_script_arg  <- if (length(args) >= 3) fs::path_expand(args[[3]]) else fs::path(base_dir, "workflow", "estimate_fdp.R")
level_arg       <- if (length(args) >= 4) tolower(args[[4]]) else "peptide"

if (!level_arg %in% c("peptide", "protein")) {
  stop("Arg #4 ('level') must be either 'peptide' or 'protein'.")
}

old_wd <- getwd()
setwd(base_dir)
on.exit(setwd(old_wd), add = TRUE)

resolve_functions_script <- function(base_dir, candidate) {
  if (fs::file_exists(candidate)) return(candidate)
  cand2 <- fs::path(base_dir, "workflow", "estimate_fdp.R")
  if (fs::file_exists(cand2)) return(cand2)
  found <- fs::dir_ls(base_dir, recurse = TRUE, type = "file", glob = "**/estimate_fdp.R")
  if (length(found) >= 1) return(found[[1]])
  stop("Could not find estimate_fdp.R.")
}

fns_script <- resolve_functions_script(base_dir, fns_script_arg)
message("Sourcing functions from: ", fns_script)
source(fns_script, chdir = TRUE)

resolve_scan_root <- function(base_dir, level = "peptide") {
  exact_candidates <- if (level == "peptide") {
    c(
      fs::path(base_dir, "output", "diann", "peptide_entrapment", "accurate-design"),
      fs::path(base_dir, "output", "diann", "peptide_entrapment", "accurate_design")
    )
  } else {
    c(
      fs::path(base_dir, "output", "diann", "protein_entrapment", "accurate-design"),
      fs::path(base_dir, "output", "diann", "protein_entrapment", "accurate_design")
    )
  }

  exact_existing <- exact_candidates[fs::dir_exists(exact_candidates)]
  if (length(exact_existing) >= 1) return(exact_existing[[1]])

  fallback_root <- if (level == "peptide") {
    fs::path(base_dir, "output", "diann", "peptide_entrapment")
  } else {
    fs::path(base_dir, "output", "diann", "protein_entrapment")
  }
  if (fs::dir_exists(fallback_root)) return(fallback_root)

  broader_root <- fs::path(base_dir, "output", "diann")
  if (fs::dir_exists(broader_root)) return(broader_root)

  stop("Could not find expected DIANN output directory.")
}

sanitize_tag <- function(x) {
  x |>
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

classify_run_dir <- function(run_dir, level = "peptide") {
  run_name <- fs::path_file(run_dir)
  run_name_lower <- stringr::str_to_lower(run_name)

  method <- dplyr::case_when(
    stringr::str_detect(run_name_lower, "library-free")   ~ "library_free",
    stringr::str_detect(run_name_lower, "library-carafe") ~ "library_carafe",
    stringr::str_detect(run_name_lower, "library-prosit") ~ "library_prosit",
    stringr::str_detect(run_name_lower, "library-ms2pip") ~ "library_ms2pip",
    TRUE                                                  ~ "unknown_method"
  )

  raws <- dplyr::case_when(
    stringr::str_detect(run_name_lower, "one[-_]?raw")   ~ "one_raw",
    stringr::str_detect(run_name_lower, "four[-_]?raw")  ~ "four_raws",
    TRUE                                                 ~ "unknown_raws"
  )

  entrapment <- dplyr::case_when(
    stringr::str_detect(run_name_lower, "random[-_]?shuffling") ~ "random_shuffling",
    stringr::str_detect(run_name_lower, "foreign[-_]?sp")       ~ "foreign_celegans",
    stringr::str_detect(run_name_lower, "c[-_]?elegans")        ~ "foreign_celegans",
    stringr::str_detect(run_name_lower, "celegans")             ~ "foreign_celegans",
    stringr::str_detect(run_name_lower, "foreign")              ~ "foreign_celegans",
    TRUE                                                        ~ "unknown_entrapment"
  )

  organism <- dplyr::case_when(
    stringr::str_detect(run_name_lower, "mouse") ~ "mouse",
    stringr::str_detect(run_name_lower, "human") ~ "human",
    TRUE                                         ~ "unknown_organism"
  )

  fdr_tag <- stringr::str_extract(run_name_lower, "[0-9]+fdr")
  if (is.na(fdr_tag)) fdr_tag <- "unknown_fdr"

  run_tag <- sanitize_tag(paste(method, raws, organism, fdr_tag, entrapment, sep = "-"))
  prefix  <- paste0("diann_", level, "_", run_tag)

  list(
    dir         = run_dir,
    run_name    = run_name,
    method      = method,
    raws        = raws,
    organism    = organism,
    entrapment  = entrapment,
    fdr_tag     = fdr_tag,
    level       = level,
    report      = fs::path(run_dir, "report.tsv"),
    prefix      = prefix,
    fdp_csv     = if (level == "peptide") fs::path(run_dir, paste0(prefix, "-diann_fdp_precursor.csv")) else fs::path(run_dir, paste0(prefix, "-diann_fdp_protein.csv")),
    pdf         = fs::path(run_dir, paste0("fdr_fdp_", run_tag, ".pdf")),
    db_log      = if (level == "peptide") fs::path(run_dir, paste0(prefix, "-selected_entrapment_db.txt")) else NA_character_
  )
}

resolve_peptide_db <- function(info, db_hint, base_dir) {
  expected_name <- dplyr::case_when(
    info$method == "library_ms2pip" && info$entrapment == "random_shuffling" ~
      "UP000000589_10090_mouse_20241210_entrapment_diann_pep_7_30.txt",
    info$method == "library_ms2pip" && info$entrapment == "foreign_celegans" ~
      "UP000000589_10090_mouse_20241210_foreignsp_diann_pep_7_30.txt",
    info$method != "library_ms2pip" && info$entrapment == "random_shuffling" ~
      "UP000000589_10090_mouse_20241210_diann_entrapment_pep.txt",
    info$method != "library_ms2pip" && info$entrapment == "foreign_celegans" ~
      "UP000000589_foreign_species_C_elegans_entrapment_diann_pep.txt",
    TRUE ~ NA_character_
  )

  if (is.na(expected_name)) {
    warning("Could not determine peptide DB for run: ", info$run_name)
    return(NA_character_)
  }

  search_dirs <- c(
    if (info$method == "library_ms2pip") fs::path(base_dir, "data", "reference-fasta", "ms2pip") else character(),
    if (info$method == "library_ms2pip") fs::path(db_hint, "ms2pip") else character(),
    if (fs::dir_exists(db_hint)) db_hint else character(),
    if (fs::file_exists(db_hint)) fs::path_dir(db_hint) else character(),
    fs::path(base_dir, "data", "reference-fasta"),
    fs::path(base_dir, "data", "databases", "diann_peptide_entrapment_database"),
    fs::path(base_dir, "data", "databases"),
    base_dir
  ) |>
    unique()

  search_dirs <- search_dirs[fs::dir_exists(search_dirs)]
  exact_hits <- fs::path(search_dirs, expected_name)
  exact_hits <- exact_hits[fs::file_exists(exact_hits)]

  if (length(exact_hits) >= 1) {
    message("Selected peptide DB for ", info$run_name, ": ", exact_hits[[1]])
    return(exact_hits[[1]])
  }

  found <- tryCatch(
    fs::dir_ls(base_dir, recurse = TRUE, type = "file", glob = paste0("**/", expected_name)),
    error = function(e) character()
  )

  if (length(found) >= 1) {
    message("Selected peptide DB for ", info$run_name, " by recursive fallback: ", found[[1]])
    return(found[[1]])
  }

  warning(
    "Could not find expected peptide DB for run '", info$run_name, "'. Expected file: ", expected_name
  )
  NA_character_
}

infer_r_from_peptide_db <- function(pep_file) {
  if (is.na(pep_file) || !fs::file_exists(pep_file)) return(NA_real_)
  x <- tryCatch(readr::read_tsv(pep_file, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (is.null(x) || !all(c("peptide_type", "decoy") %in% names(x))) return(NA_real_)
  x <- x %>% dplyr::filter(!stringr::str_detect(.data$decoy, "^Yes"))
  n_t <- sum(x$peptide_type == "target", na.rm = TRUE)
  n_p <- sum(x$peptide_type == "p_target", na.rm = TRUE)
  if (n_t <= 0 || n_p <= 0) return(NA_real_)
  n_p / n_t
}

save_plot_safely <- function(p0, pdf_file, run_name) {
  ok <- FALSE
  if (capabilities("cairo")) {
    try({
      grDevices::cairo_pdf(pdf_file, width = 4.5, height = 4.5, family = "sans")
      print(p0 + ggplot2::theme(text = ggplot2::element_text(family = "sans")))
      grDevices::dev.off()
      ok <- TRUE
    }, silent = TRUE)
  }
  if (!ok) {
    try({
      grDevices::pdf(pdf_file, width = 4.5, height = 4.5, family = "Helvetica")
      print(p0 + ggplot2::theme(text = ggplot2::element_text(family = "Helvetica")))
      grDevices::dev.off()
      ok <- TRUE
    }, silent = TRUE)
  }
  if (!ok) {
    png_file <- fs::path_ext_set(pdf_file, "png")
    warning("PDF plotting failed for ", run_name, ". Saving PNG instead: ", png_file)
    grDevices::png(png_file, width = 1200, height = 1200, res = 200)
    print(p0 + ggplot2::theme(text = ggplot2::element_text(family = "sans")))
    grDevices::dev.off()
  }
}


make_fdp_fdr_plot <- function(fdp_csv,
                              level = "peptide",
                              fdr_max = NULL,
                              add_numbers = FALSE,
                              fixed_fdr_max = TRUE) {
  # Prefer the user's existing plotting helper if available. It already matches
  # the desired visual style used in this project.
  if (exists("plot_fdp_fdr_modified", mode = "function")) {
    p_try <- tryCatch(
      plot_fdp_fdr_modified(
        fdp_csv,
        fdr_max = fdr_max,
        add_numbers = add_numbers,
        fixed_fdr_max = fixed_fdr_max
      ),
      error = function(e) NULL
    )
    if (!is.null(p_try)) return(p_try)
  }

  # Fallback plotter that works for both peptide and protein FDRBench CSVs.
  x <- readr::read_csv(fdp_csv, show_col_types = FALSE)

  fdr_col <- dplyr::case_when(
    "q_value" %in% names(x) ~ "q_value",
    "fdr" %in% names(x) ~ "fdr",
    "FDR" %in% names(x) ~ "FDR",
    "fdr_threshold" %in% names(x) ~ "fdr_threshold",
    TRUE ~ NA_character_
  )

  if (is.na(fdr_col)) {
    stop("Could not find an FDR/q-value column in: ", fdp_csv)
  }

  method_cols <- intersect(
    c("combined_fdp", "paired_fdp", "lower_bound_fdp"),
    names(x)
  )

  if (length(method_cols) == 0) {
    stop("Could not find FDP columns in: ", fdp_csv)
  }

  if (!is.null(fdr_max)) {
    x <- x[x[[fdr_col]] <= fdr_max, , drop = FALSE]
  }

  long <- dplyr::bind_rows(lapply(method_cols, function(mc) {
    data.frame(
      fdr = x[[fdr_col]],
      fdp = x[[mc]],
      method = dplyr::case_when(
        mc == "combined_fdp" ~ "Combined method",
        mc == "paired_fdp" ~ "Paired method",
        mc == "lower_bound_fdp" ~ "Lower bound",
        TRUE ~ mc
      ),
      stringsAsFactors = FALSE
    )
  }))

  # Some FDRBench outputs store FDP as fractions. Plot in percent.
  long$FDR_percent <- long$fdr * 100
  long$FDP_percent <- long$fdp * 100

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$FDR_percent, y = .data$FDP_percent, color = .data$method)
  ) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::labs(
      x = "FDR threshold",
      y = "Estimated FDP",
      color = "Method",
      title = paste0("FDP vs FDR: ", level, " level")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  p
}

add_protein_fdp_labels <- function(p0, fdp_csv, fdr_limit = 0.01) {
  x <- readr::read_csv(fdp_csv, show_col_types = FALSE)

  fdr_col <- dplyr::case_when(
    "q_value" %in% names(x) ~ "q_value",
    "fdr" %in% names(x) ~ "fdr",
    "FDR" %in% names(x) ~ "FDR",
    "fdr_threshold" %in% names(x) ~ "fdr_threshold",
    TRUE ~ NA_character_
  )

  if (is.na(fdr_col)) {
    warning("Could not find FDR/q-value column for plot labels: ", fdp_csv)
    return(p0)
  }

  x <- x |>
    dplyr::filter(.data[[fdr_col]] <= fdr_limit) |>
    dplyr::arrange(dplyr::desc(.data[[fdr_col]]))

  if (nrow(x) == 0) {
    warning("No FDP rows found at or below FDR ", fdr_limit, " in: ", fdp_csv)
    return(p0)
  }

  row1 <- x[1, , drop = FALSE]

  id_col <- intersect(
    c(
      "total_discoveries",
      "n_discoveries",
      "discoveries",
      "num_discoveries",
      "total_ids",
      "n_ids",
      "ids"
    ),
    names(row1)
  )

  total_ids <- if (length(id_col) >= 1) {
    as.character(row1[[id_col[1]]])
  } else {
    "NA"
  }

  fmt_fdp <- function(v) {
    if (length(v) == 0 || is.na(v)) return("NA")
    sprintf("%.2f%%", as.numeric(v) * 100)
  }

  combined_fdp <- if ("combined_fdp" %in% names(row1)) fmt_fdp(row1$combined_fdp) else "NA"
  paired_fdp   <- if ("paired_fdp" %in% names(row1)) fmt_fdp(row1$paired_fdp) else "NA"
  lower_fdp    <- if ("lower_bound_fdp" %in% names(row1)) fmt_fdp(row1$lower_bound_fdp) else "NA"

  label_text <- paste0(
    "Total discoveries: ", total_ids, "\n",
    "Combined method: ", combined_fdp, "\n",
    "Paired method: ", paired_fdp, "\n",
    "Lower bound: ", lower_fdp
  )

  p0 +
    ggplot2::annotate(
      "text",
      x = 0.00025,
      y = 0.00975,
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      lineheight = 1.05
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, fdr_limit),
      ylim = c(0, fdr_limit),
      clip = "off"
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

# plot_fdp_if_available <- function(info, level_label) {
#   if (fs::file_exists(info$fdp_csv)) {
#     message("Saving ", level_label, " FDP vs FDR plot: ", info$pdf)
#     p0 <- make_fdp_fdr_plot(
#       info$fdp_csv,
#       level = level_label,
#       fdr_max = NULL,
#       add_numbers = TRUE,
#       fixed_fdr_max = TRUE
#     )

#     # Change x- and y-lim to [0,1]
#     if (level_label == "protein") {
#       p0 <- p0 +
#         ggplot2::coord_cartesian(
#           xlim = c(0, 0.05),
#           ylim = c(0, 0.05)
#         )
#     }
#     save_plot_safely(p0, info$pdf, info$run_name)
#   } else {
#     warning("Expected FDP CSV not found for ", level_label, ": ", info$fdp_csv)
#   }
# }

plot_fdp_if_available <- function(info, level_label) {
  if (fs::file_exists(info$fdp_csv)) {
    message("Saving ", level_label, " FDP vs FDR plot: ", info$pdf)

    if (level_label %in% c("peptide", "protein")) {

      # For protein entrapment:
      # - calculate/report FDP values at 1% FDR
      # - display plot from 0% to 5%
      summary_fdr <- 0.01
      plot_limit  <- 0.05

      p0 <- make_fdp_fdr_plot(
        info$fdp_csv,
        level = level_label,
        fdr_max = plot_limit,
        add_numbers = TRUE,
        fixed_fdr_max = FALSE
      )

      p0 <- p0 +
        ggplot2::coord_cartesian(
          xlim = c(0, plot_limit),
          ylim = c(0, plot_limit),
          clip = "off"
        )

      save_plot_safely(p0, info$pdf, info$run_name)

    } else {

      p0 <- make_fdp_fdr_plot(
        info$fdp_csv,
        level = level_label,
        fdr_max = NULL,
        add_numbers = TRUE,
        fixed_fdr_max = TRUE
      )

      save_plot_safely(p0, info$pdf, info$run_name)
    }

  } else {
    warning("Expected FDP CSV not found for ", level_label, ": ", info$fdp_csv)
  }
}

process_one <- function(parquet_path, db_hint, base_dir, level = "peptide") {
  run_dir <- fs::path_dir(parquet_path)
  info <- classify_run_dir(run_dir, level = level)

  message("\n=== Processing: ", parquet_path, " ===")
  message("Run folder: ", info$run_name)
  message("Detected method: ", info$method,
          " | entrapment: ", info$entrapment,
          " | level: ", info$level)

  dt <- arrow::read_parquet(parquet_path)
  readr::write_tsv(dt, info$report)

  if (info$level == "peptide") {
    pep_db_file <- resolve_peptide_db(info, db_hint = db_hint, base_dir = base_dir)

    if (!is.na(pep_db_file)) writeLines(pep_db_file, info$db_log)

    if (is.na(pep_db_file) || !fs::file_exists(pep_db_file)) {
      warning("Peptide entrapment DB not found; skipping: ", info$run_name)
      return(invisible(FALSE))
    }

    message("Using peptide entrapment DB: ", pep_db_file)

    r_arg <- NULL
    pep_file_for_fdp <- pep_db_file

    if (info$entrapment == "foreign_celegans") {
      r_arg <- infer_r_from_peptide_db(pep_db_file)
      if (is.na(r_arg)) r_arg <- 0.7

      # For foreign-species peptide entrapment, FDREval's combined method can
      # fail if the real peptide table does not contain every DIANN-reported
      # peptide. This is especially common for library-free runs. Passing this
      # one-line pointer file makes FDREval fall back to Protein.Group-based
      # target/entrapment classification while still using the correct r value.
      #
      # The real selected DB path remains recorded in info$db_log.
      foreign_pointer_file <- fs::path(info$dir, paste0(info$prefix, "-foreign_combined_fallback_pointer.txt"))
      writeLines(pep_db_file, foreign_pointer_file)
      pep_file_for_fdp <- foreign_pointer_file

      message("Foreign-species entrapment: using combined method with r=", signif(r_arg, 4))
      message("Foreign-species FDP classification mode: Protein.Group fallback via pointer file: ", foreign_pointer_file)
    }

    out <- tryCatch(
      run_diann_fdp_analysis(
        report_file      = info$report,
        level            = "peptide",
        pep_file         = pep_file_for_fdp,
        prefix           = info$prefix,
        k_fold           = 1,
        out_dir          = info$dir,
        r                = r_arg,
        entrapment_type  = info$entrapment
      ),
      error = function(e) {
        if (stringr::str_detect(conditionMessage(e), "unused argument")) {
          warning("run_diann_fdp_analysis() does not accept entrapment_type; retrying without it.")
          run_diann_fdp_analysis(
            report_file = info$report,
            level       = "peptide",
            pep_file    = pep_file_for_fdp,
            prefix      = info$prefix,
            k_fold      = 1,
            out_dir     = info$dir,
            r           = r_arg
          )
        } else {
          stop(e)
        }
      }
    )

    plot_fdp_if_available(info, "peptide")
  } else {
    # Use the foreign-species protein ratio when appropriate
    # Random shuffling runs retain the one-fold paired method
    r_arg <- if (info$entrapment == "foreign_celegans") 0.7 else NULL 

    if (is.null(r_arg)){
      message("Random shuffling protein entrapment: using paired method (no r value needed).")
    } else {
      message("Foreign-species protein entrapment: using combined method with r=", signif(r_arg, 4))
    }

    run_diann_fdp_analysis(
      report_file = info$report,
      level       = "protein",
      pep_file    = NULL,
      prefix      = info$prefix,
      k_fold      = 1,
      out_dir     = info$dir,
      r           = r_arg,
      entrapment_type = info$entrapment
    )

    plot_fdp_if_available(info, "protein")
  }

  invisible(TRUE)
}

scan_root <- resolve_scan_root(base_dir, level = level_arg)
message("Scanning for DIANN parquet results under: ", scan_root)

all_parquets <- fs::dir_ls(scan_root, recurse = TRUE, type = "file", glob = "**/result.parquet") |>
  purrr::discard(~ stringr::str_detect(.x, "result-first-pass\\.parquet$")) |>
  unique()

if (length(all_parquets) == 0L) {
  message("No result.parquet files found under: ", scan_root)
  quit(status = 0)
}

message("Found ", length(all_parquets), " parquet file(s).")
purrr::walk(all_parquets, ~ process_one(.x, db_hint_arg, base_dir = base_dir, level = level_arg))
message("\nAll done.")

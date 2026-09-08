#!/usr/bin/env Rscript

# Portions of this script were adapted from Noble-Lab/FDRBench:
# https://github.com/Noble-Lab/FDRBench
# Licensed under the Apache License 2.0.
# Modifications were made for the DIA-NN analyses in this repository.

options(repr.matrix.max.cols=150, repr.matrix.max.rows=200)

library(readr)
library(ggplot2)
library(DBI)
library(ggpubr)
library(tidyr)
library(dplyr)
library(stringr)
library(here)
library(tidyverse)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

get_fdrbench_jar <- function() {
  jar <- Sys.getenv("FDRBENCH_JAR", unset = "")

  if (!nzchar(jar)) {
    stop(
      "Set FDRBENCH_JAR to the location of the FDRBench .jar file before running this analysis."
    )
  }

  jar <- path.expand(jar)
  if (!file.exists(jar)) {
    stop("FDRBench .jar not found: ", jar)
  }

  normalizePath(jar)
}

# -----------------------------------------------------------------------------
# Helper: infer entrapment type from peptide-pair file name
# -----------------------------------------------------------------------------
detect_entrapment_type_from_pep_file <- function(pep_file) {
  if (is.null(pep_file) || is.na(pep_file) || !nzchar(pep_file)) {
    return("none")
  }

  pep_name <- tolower(basename(pep_file))

  if (grepl("random|shuff", pep_name)) {
    return("random_shuffling")
  }

  if (grepl("celegans|c_elegans|foreign", pep_name)) {
    return("foreign_celegans")
  }

  return("unknown")
}

# -----------------------------------------------------------------------------
# Main function
# report_file: DIANN report.tsv file
# pep_file: peptide-level pair / entrapment file for FDRBench
# -----------------------------------------------------------------------------
run_diann_fdp_analysis <- function(report_file = "",
                                   level = "protein",
                                   pep_file = NULL,
                                   prefix = "test",
                                   k_fold = 1,
                                   pick_one_protein_method = "first",
                                   out_dir = NULL,
                                   r = NULL,
                                   entrapment_type = "auto") {

  a <- read_tsv(report_file, show_col_types = FALSE)
  n_run <- a %>% select(Run) %>% distinct() %>% nrow()

  if (level == "protein") {
    # Apply the foreign species entrapment ratio at protein level
    if (entrapment_type == "auto" && !is.null(pep_file)) {
      entrapment_type <- detect_entrapment_type_from_pep_file(pep_file)
    }

    if (entrapment_type == "foreign_celegans" && is.null(r)) {
      r <- 0.7
    }

    cat("Protein-level entrapment type:", entrapment_type, "\n")
    cat("Protein-level entrapment ratio (r):", r, "\n")
    
    if (n_run >= 2) {
      cat("Multiple runs in the report file:", n_run, "\n")
      b <- a %>% select(`Protein.Group`, `Lib.PG.Q.Value`) %>% distinct()
      b$q_value <- b$Lib.PG.Q.Value
    } else {
      cat("Single run in the report file\n")
      b <- a %>% select(`Protein.Group`, `PG.Q.Value`) %>% distinct()
      b$q_value <- b$PG.Q.Value
    }

    b$protein <- b$Protein.Group
    set.seed(2024)
    b$score <- sample(x = seq_len(nrow(b)), size = nrow(b), replace = FALSE)
    b <- b %>% arrange(q_value, score) %>% mutate(score = row_number())

    cat("The number of proteins:", nrow(b), "\n")

    if (is.null(out_dir)) {
      out_dir <- dirname(report_file)
    } else {
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
      }
    }

    out_file <- paste(out_dir, "/", prefix, "-fdp_protein_input.tsv", sep = "")
    write_tsv(b, out_file)

    fdp_file <- paste(out_dir, "/", prefix, "-diann_fdp_protein.csv", sep = "")

    if (!is.null(r)) {
      cmd <- paste(
        "java -jar", shQuote(get_fdrbench_jar()),
        "-i", shQuote(out_file),
        "-level protein",
        "-o", shQuote(fdp_file),
        "-score 'score:0'",
        "-r", r,
        "-pick", pick_one_protein_method
      )
    } else {
      cmd <- paste(
        "java -jar", shQuote(get_fdrbench_jar()),
        "-i", shQuote(out_file),
        "-level protein",
        "-o", shQuote(fdp_file),
        "-score 'score:0'",
        "-fold", k_fold,
        "-pick", pick_one_protein_method
      )
    }

    cat("Running ", cmd, "\n")
    out <- system(cmd, intern = TRUE)
    cat(paste(out, collapse = "\n"), "\n")
    return(fdp_file)
  }

  # else if (level == "peptide" || level == "precursor") {
  #   if (n_run >= 2) {
  #     cat("Multiple runs in the report file:", n_run, "\n")
  #     b <- a %>%
  #       select(`Run`, `Stripped.Sequence`, `Modified.Sequence`, `Precursor.Charge`,
  #              `Lib.Q.Value`, `PEP`, `Protein.Group`) %>%
  #       rename(
  #         q_value     = `Lib.Q.Value`,
  #         run         = `Run`,
  #         peptide     = `Stripped.Sequence`,
  #         mod_peptide = `Modified.Sequence`,
  #         charge      = `Precursor.Charge`,
  #         protein     = `Protein.Group`
  #       )

  #     # keep the top precursor
  #     b <- b %>%
  #       group_by(peptide, mod_peptide, charge) %>%
  #       arrange(q_value, PEP) %>%
  #       filter(row_number() == 1) %>%
  #       ungroup()
  #   } else {
  #     cat("Single run in the report file\n")
  #     b <- a %>%
  #       select(`Run`, `Stripped.Sequence`, `Modified.Sequence`, `Precursor.Charge`,
  #              `Q.Value`, `PEP`, `Protein.Group`) %>%
  #       distinct() %>%
  #       rename(
  #         q_value     = `Q.Value`,
  #         run         = `Run`,
  #         peptide     = `Stripped.Sequence`,
  #         mod_peptide = `Modified.Sequence`,
  #         charge      = `Precursor.Charge`,
  #         protein     = `Protein.Group`
  #       )
  #   }

  #   set.seed(2024)
  #   b <- b %>% arrange(q_value, PEP) %>% mutate(score = row_number())

  #   cat("The number of peptides:", nrow(b), "\n")

  #   if (is.null(out_dir)) {
  #     out_dir <- dirname(report_file)
  #   } else {
  #     if (!dir.exists(out_dir)) {
  #       dir.create(out_dir, recursive = TRUE)
  #     }
  #   }

  #   out_file <- paste(out_dir, "/", prefix, "-fdp_precursor_input.tsv", sep = "")
  #   message("Rows before precursor deduplication: ", nrow(a))

  #   dup_summary <- a %>%
  #     dplyr::count(peptide, charge, name = "n") %>%
  #     dplyr::filter(n > 1)

  #   message("Duplicated peptide+charge keys before deduplication: ", nrow(dup_summary))

  #   if (nrow(dup_summary) > 0) {
  #     message("Maximum duplicate count for peptide+charge: ", max(dup_summary$n))
  #   }

  #   a <- a %>%
  #     dplyr::mutate(
  #       peptide = as.character(peptide),
  #       charge = as.integer(charge),
  #       q_value = as.numeric(q_value),
  #       PEP = as.numeric(PEP),
  #       score = as.numeric(score),
  #       q_value_sort = ifelse(is.na(q_value), Inf, q_value),
  #       PEP_sort = ifelse(is.na(PEP), Inf, PEP),
  #       score_sort = ifelse(is.na(score), Inf, score)
  #     ) %>%
  #     dplyr::arrange(q_value_sort, PEP_sort, score_sort) %>%
  #     dplyr::group_by(peptide, charge) %>%
  #     dplyr::slice(1) %>%
  #     dplyr::ungroup() %>%
  #     dplyr::select(-q_value_sort, -PEP_sort, -score_sort)

  #   message("Rows after precursor deduplication: ", nrow(a))
  #   write_tsv(b, out_file)

  #   fdp_file <- paste(out_dir, "/", prefix, "-diann_fdp_precursor.csv", sep = "")

  #   if (!is.null(pep_file)) {
  #     if (entrapment_type == "auto") {
  #       entrapment_type <- detect_entrapment_type_from_pep_file(pep_file)
  #     }

  #     cat("Detected entrapment type:", entrapment_type, "\n")
  #     cat("Using peptide entrapment file:", pep_file, "\n")

  #     # IMPORTANT:
  #     # - random shuffling works with the paired k-fold route (-fold ...)
  #     # - foreign C. elegans appears to fail in FDRBench's paired peptide lookup
  #     #   path, so route it through the -r mode instead of -fold
  #     if (entrapment_type == "foreign_celegans") {
  #       if (is.null(r)) {
  #         r <- 0.7
  #       }

  #       cmd <- paste(
  #         "java -jar /path/to/fdrbench.jar",
  #         "-i", out_file,
  #         "-pep", pep_file,
  #         "-level precursor",
  #         "-o", fdp_file,
  #         "-r", r,
  #         "-score 'score:0'"
  #       )
  #     } else {
  #       if (!is.null(r)) {
  #         cmd <- paste(
  #           "java -jar /path/to/fdrbench.jar",
  #           "-i", out_file,
  #           "-pep", pep_file,
  #           "-level precursor",
  #           "-o", fdp_file,
  #           "-r", r,
  #           "-score 'score:0'"
  #         )
  #       } else {
  #         cmd <- paste(
  #           "java -jar /path/to/fdrbench.jar",
  #           "-i", out_file,
  #           "-fold", k_fold,
  #           "-pep", pep_file,
  #           "-level precursor",
  #           "-o", fdp_file,
  #           "-score 'score:0'"
  #         )
  #       }
  #     }

  #     cat("Running ", cmd, "\n")
  #     out <- system(cmd, intern = TRUE)
  #     cat(paste(out, collapse = "\n"), "\n")
  #     return(fdp_file)
  #   } else {
  #     cat("No paired peptide file\n")
  #   }
  # }
  else if (level == "peptide" || level == "precursor") {
  if (n_run >= 2) {
    cat("Multiple runs in the report file:", n_run, "\n")

    b <- a %>%
      select(
        `Run`,
        `Stripped.Sequence`,
        `Modified.Sequence`,
        `Precursor.Charge`,
        `Lib.Q.Value`,
        `PEP`,
        `Protein.Group`
      ) %>%
      rename(
        q_value     = `Lib.Q.Value`,
        run         = `Run`,
        peptide     = `Stripped.Sequence`,
        mod_peptide = `Modified.Sequence`,
        charge      = `Precursor.Charge`,
        protein     = `Protein.Group`
      )
  } else {
    cat("Single run in the report file\n")

    b <- a %>%
      select(
        `Run`,
        `Stripped.Sequence`,
        `Modified.Sequence`,
        `Precursor.Charge`,
        `Q.Value`,
        `PEP`,
        `Protein.Group`
      ) %>%
      rename(
        q_value     = `Q.Value`,
        run         = `Run`,
        peptide     = `Stripped.Sequence`,
        mod_peptide = `Modified.Sequence`,
        charge      = `Precursor.Charge`,
        protein     = `Protein.Group`
      )
  }

  # Basic cleanup
  b <- b %>%
    mutate(
      peptide = as.character(peptide),
      mod_peptide = as.character(mod_peptide),
      charge = as.integer(charge),
      q_value = as.numeric(q_value),
      PEP = as.numeric(PEP),
      protein = as.character(protein)
    ) %>%
    filter(
      !is.na(peptide),
      !is.na(charge),
      !is.na(q_value)
    )

  # First remove exact duplicate rows
  b <- b %>%
    distinct(
      run,
      peptide,
      mod_peptide,
      charge,
      q_value,
      PEP,
      protein,
      .keep_all = TRUE
    )

  # Assign score before final precursor-level deduplication
  # Lower score is better for FDRBench because you pass -score 'score:0'
  b <- b %>%
    arrange(q_value, PEP) %>%
    mutate(score = row_number())

  cat("The number of precursor rows before peptide+charge deduplication:", nrow(b), "\n")

  dup_summary <- b %>%
    count(peptide, charge, name = "n") %>%
    filter(n > 1)

  message("Duplicated peptide+charge keys before deduplication: ", nrow(dup_summary))

  if (nrow(dup_summary) > 0) {
    message("Maximum duplicate count for peptide+charge: ", max(dup_summary$n))
  }

  # Important:
  # FDRBench paired/k-fold precursor mode uses stripped peptide + charge.
  # Therefore, keep only one row per peptide + charge.
  b <- b %>%
    mutate(
      q_value_sort = ifelse(is.na(q_value), Inf, q_value),
      PEP_sort = ifelse(is.na(PEP), Inf, PEP),
      score_sort = ifelse(is.na(score), Inf, score)
    ) %>%
    arrange(q_value_sort, PEP_sort, score_sort) %>%
    group_by(peptide, charge) %>%
    slice(1) %>%
    ungroup() %>%
    select(-q_value_sort, -PEP_sort, -score_sort)

  # Re-rank after deduplication so score is continuous and unique
  b <- b %>%
    arrange(q_value, PEP) %>%
    mutate(score = row_number())

  message("Rows after precursor deduplication: ", nrow(b))

  if (is.null(out_dir)) {
    out_dir <- dirname(report_file)
  } else {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
  }

  out_file <- paste(out_dir, "/", prefix, "-fdp_precursor_input.tsv", sep = "")
  write_tsv(b, out_file)

  fdp_file <- paste(out_dir, "/", prefix, "-diann_fdp_precursor.csv", sep = "")

  if (!is.null(pep_file)) {
    if (entrapment_type == "auto") {
      entrapment_type <- detect_entrapment_type_from_pep_file(pep_file)
    }

    cat("Detected entrapment type:", entrapment_type, "\n")
    cat("Using peptide entrapment file:", pep_file, "\n")

    if (entrapment_type == "foreign_celegans") {
      if (is.null(r)) {
        r <- 0.7
      }

      cmd <- paste(
        "java -jar", shQuote(get_fdrbench_jar()),
        "-i", shQuote(out_file),
        "-pep", shQuote(pep_file),
        "-level precursor",
        "-o", shQuote(fdp_file),
        "-r", r,
        "-score 'score:0'"
      )
    } else {
      if (!is.null(r)) {
        cmd <- paste(
          "java -jar", shQuote(get_fdrbench_jar()),
          "-i", shQuote(out_file),
          "-pep", shQuote(pep_file),
          "-level precursor",
          "-o", shQuote(fdp_file),
          "-r", r,
          "-score 'score:0'"
        )
      } else {
        cmd <- paste(
          "java -jar", shQuote(get_fdrbench_jar()),
          "-i", shQuote(out_file),
          "-fold", k_fold,
          "-pep", shQuote(pep_file),
          "-level precursor",
          "-o", shQuote(fdp_file),
          "-score 'score:0'"
        )
      }
    }

    cat("Running ", cmd, "\n")
    out <- system(cmd, intern = TRUE)
    cat(paste(out, collapse = "\n"), "\n")
    return(fdp_file)
  } else {
    cat("No paired peptide file\n")
    return(out_file)
   }
  }
}

# report_file: report.tsv file, pep_file: peptide level pair file
# run_diann_fdp_analysis=function(report_file="",level="protein",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
#   a <- read_tsv(report_file)
#   n_run <- a %>% select(Run) %>% distinct() %>% nrow()
#   if(level=="protein"){
#     if(n_run>=2){
#       ## multiple runs
#       cat("Multiple runs in the report file:",n_run,"\n")
#       b <- a %>% select(`Protein.Group`,`Lib.PG.Q.Value`) %>% distinct()
#       b$q_value <- b$Lib.PG.Q.Value
#     }else{
#       cat("Single run in the report file\n")
#       b <- a %>% select(`Protein.Group`,`PG.Q.Value`) %>% distinct()
#       b$q_value <- b$PG.Q.Value
#     }
#     #b$protein <- sapply(b$Protein.Group,function(x){ y<-  str_split(x,pattern = ";") %>% unlist; y[ length(y)]})
#     b$protein <- b$Protein.Group
#     set.seed(2024)
#     b$score <- sample(x = 1:nrow(b),size = nrow(b),replace = FALSE)
#     b <- b %>% arrange(q_value,score) %>% mutate(score=row_number())
    
#     cat("The number of proteins:",nrow(b),"\n")
    
#     if(is.null(out_dir)){
#       out_dir <- dirname(report_file)
#     }else{
#       # create the directory if it does not exist
#       if(!dir.exists(out_dir)){
#         dir.create(out_dir)
#       }
#     }
#     out_file <- paste(out_dir,"/",prefix,"-fdp_protein_input.tsv",sep="")
#     write_tsv(b,out_file)
#     fdp_file <- paste(out_dir,"/",prefix,"-diann_fdp_protein.csv",sep="")
    
#     if(!is.null(r)){
#       cmd <- paste("java -jar /path/to/fdrbench.jar -i ", out_file, " -level protein -o ",fdp_file, " -score 'score:0'"," -r ",r," -pick ",pick_one_protein_method,sep="")
#     }else{
#       cmd <- paste("java -jar /path/to/fdrbench.jar -i ", out_file, " -level protein -o ",fdp_file, " -score 'score:0'"," -fold ",k_fold," -pick ",pick_one_protein_method,sep="")
#     }
#     cat("Running ",cmd,"\n")
#     out <- system(cmd,intern = TRUE)
#     cat(paste(out,collapse = "\n"),"\n")
#     return(fdp_file)
#   }else if(level=="peptide" || level=="precursor"){ 
#     if(n_run>=2){
#       ## multiple runs
#       cat("Multiple runs in the report file:",n_run,"\n")
#       b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Lib.Q.Value`,`PEP`,`Protein.Group`) %>% 
#         rename(q_value=`Lib.Q.Value`,run=`Run`,peptide=`Stripped.Sequence`,mod_peptide=`Modified.Sequence`,charge=`Precursor.Charge`,protein=`Protein.Group`)
      
#       ## keep the top precursor
#       b <- b %>% group_by(peptide,mod_peptide,charge) %>% arrange(q_value,PEP) %>% filter(row_number()==1) %>% ungroup()
#     }else{
#       cat("Single run in the report file\n")
#       # b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Global.Q.Value`,`PEP`,`Protein.Group`) %>% distinct()
#       # b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Global.Q.Value`,`PEP`,`Protein.Group`) %>% 
#       #   rename(q_value=`Global.Q.Value`,run=`Run`,peptide=`Stripped.Sequence`,mod_peptide=`Modified.Sequence`,charge=`Precursor.Charge`,protein=`Protein.Group`)
#       b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Q.Value`,`PEP`,`Protein.Group`) %>% distinct()
#       b <- a %>% select(`Run`,`Stripped.Sequence`,`Modified.Sequence`,`Precursor.Charge`,`Q.Value`,`PEP`,`Protein.Group`) %>% 
#         rename(q_value=`Q.Value`,run=`Run`,peptide=`Stripped.Sequence`,mod_peptide=`Modified.Sequence`,charge=`Precursor.Charge`,protein=`Protein.Group`)
#     }
#     set.seed(2024)
#     b <- b %>% arrange(q_value,PEP) %>% mutate(score=row_number())
    
#     cat("The number of peptides:",nrow(b),"\n")
    
#     if(is.null(out_dir)){
#       out_dir <- dirname(report_file)
#     }else{
#       # create the directory if it does not exist
#       if(!dir.exists(out_dir)){
#         dir.create(out_dir)
#       }
#     }
#     out_file <- paste(out_dir,"/",prefix,"-fdp_precursor_input.tsv",sep="")
#     write_tsv(b,out_file)
#     fdp_file <- paste(out_dir,"/",prefix,"-diann_fdp_precursor.csv",sep="")
    
#     if(!is.null(pep_file)){
#       if(!is.null(r)){
#         cmd <- paste("java -jar /path/to/fdrbench.jar -i ", out_file, " -pep ", pep_file, " -level precursor -o ",fdp_file, " -r ",r," -score 'score:0'",sep="")
#       }else{
#         cmd <- paste("java -jar /path/to/fdrbench.jar -i ", out_file," -fold ",k_fold, " -pep ", pep_file, " -level precursor -o ",fdp_file, " -score 'score:0'",sep="")
#       }
#       cat("Running ",cmd,"\n")
#       out <- system(cmd,intern = TRUE)
#       cat(paste(out,collapse = "\n"),"\n")
#       return(fdp_file)
#     }else{
#       cat("No paired peptide file\n")
#     }
    
#   }
  
# }

## only peptide level is supported
run_fdp_analysis=function(report_file="",level="peptide",pep_file=NULL,prefix="test",k_fold=1,pick_one_protein_method="first",out_dir=NULL,r=NULL) {
  a <- read_tsv(report_file)
  if(level=="peptide"){         
    set.seed(2024)
    b <- a %>% arrange(q_value,score) %>% mutate(score=row_number())
    cat("The number of peptides:",nrow(b),"\n")
    cat("The number of peptides passed 1% FDR:",nrow(b %>% filter(q_value<=0.01)),"\n")
    if(is.null(out_dir)){
      out_dir <- dirname(report_file)
    }else{
      # create the directory if it does not exist
      if(!dir.exists(out_dir)){
        dir.create(out_dir)
      }
    }
    out_file <- paste(out_dir,"/",prefix,"-fdp_peptide_input.tsv",sep="")
    write_tsv(b,out_file)
    fdp_file <- paste(out_dir,"/",prefix,"-std_fdp_peptide.csv",sep="")
    
    if(!is.null(pep_file)){
      if(!is.null(r)){
        cmd <- paste("java -jar", shQuote(get_fdrbench_jar()), "-i", shQuote(out_file), "-pep", shQuote(pep_file), "-level peptide -o", shQuote(fdp_file), "-r", r, "-score 'score:0'")
      }else{
        cmd <- paste("java -jar", shQuote(get_fdrbench_jar()), "-i", shQuote(out_file), "-fold", k_fold, "-pep", shQuote(pep_file), "-level peptide -o", shQuote(fdp_file), "-score 'score:0'")
      }
      cat("Running ",cmd,"\n")
      out <- system(cmd,intern = TRUE)
      cat(paste(out,collapse = "\n"),"\n")
      return(fdp_file)
    }else{
      cat("No paired peptide file\n")
    }
    
  }
  
}


plot_fdp_fdr=function(fdp_file="",fdr_max=NULL,fig_title=NULL,scale_xy=TRUE,add_numbers=FALSE) {
  x <- read_csv(fdp_file)
  if("FDP_1B" %in% names(x)){
    dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
      rename(`Combined entrapment`=FDP,`Paired entrapment`=FDP_1B) %>%
      gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
  }else{
    dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_min) %>% distinct() %>% 
      rename(`Combined entrapment`=FDP) %>%
      gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
  }
  
  
  max_fdp <- max(c(dat$FDP,dat$q_value))
  if(!is.null(fdr_max)){
    max_fdp <- min(c(fdr_max,max_fdp))
  }
  gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
    geom_abline(slope = 1,intercept = 0,color="gray")+
    #rasterise(geom_line(), dpi = 300) + 
    geom_line()+
    xlab("FDR")+
    ylab("FDP")+
    theme_bw()+
    #geom_segment(x=0.01,xend=0.01,y=0,yend=0.01,linewidth=0.3,color="blue",linetype=2)+
    theme_pubr(base_size = 12,border = TRUE)
  
  if(scale_xy){
    gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
      xlim(0,max_fdp)+
      ylim(0,max_fdp)+
      scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
      #scale_y_pct()+
      #scale_x_pct()+
      scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
  }
  #theme(legend.position = "top",plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"))+    
  ## legend on the botton right
  ## no background color for the legend
  gg1 <- gg1 + theme(legend.position = c(0.65, 0.16),plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank())
  
  if(!is.null(fig_title)){
    gg1 <- gg1 + ggtitle(fig_title)
  }
  
  if(add_numbers){
    # add numbers on the top left size of the figure using annotation, text align to left
    y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
    #n_t <- 
    if(abs(max_fdp-0.01)<=0.02){
      # text right align
      gg1 <- gg1 + annotate("text", x = max_fdp*0.1, y = 0.9*max_fdp, label = paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep=""), color = "black", size = 3,hjust = 0)
    }else{
      gg1 <- gg1 + annotate("text", x = 0.01, y = 0.9*max_fdp, label = paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep=""), color = "black", size = 3,hjust = 0)
    }
  }
  
  # library(ggpubr)
  # options(repr.plot.width = 12, repr.plot.height = 6)
  # gg <- ggarrange(gg1,gg2,ncol = 2,  common.legend = TRUE)
  # options(jupyter.plot_mimetypes = "image/png")
  return(gg1)
  
  #pdf("fdp_fdr_lib_qvalue_protein.pdf",width = 8,height = 4)
  #print(gg)
  #dev.off()
}




plot_fdp_fdr_v2=function(fdp_file="",fdr_max=NULL,fig_title=NULL,scale_xy=TRUE,add_numbers=FALSE,r=1,fixed_fdr_max=FALSE,max_x=NA,max_y=NA,
                         color_mapping=NULL,
                         legend_position=NULL,
                         fdr_decimal_place=2,
                         return_data=FALSE,
                         text_position=NULL,
                         add_max_qvalue=FALSE) {
  x <- read_csv(fdp_file)
  if("FDP_1B" %in% names(x)){
    if(r>=2){
      dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
        rename(`Combined method`=FDP,`Matched method`=FDP_1B,`Lower bound`=FDP_min) %>%
        gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
    }else{
      dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_1B,FDP_min) %>% distinct() %>% 
        rename(`Combined method`=FDP,`Paired method`=FDP_1B,`Lower bound`=FDP_min) %>%
        gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
    }
  }else{
    dat <- x %>% mutate(FDP_min=n_p/(n_p+n_t)) %>% select(q_value,FDP,FDP_min) %>% distinct() %>% 
      rename(`Combined method`=FDP,`Lower bound`=FDP_min) %>%
      gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
  }
  
  
  max_fdp <- max(c(dat$FDP,dat$q_value))
  if(!is.null(fdr_max)){
    if(fixed_fdr_max){
      max_fdp <- fdr_max
    }else{
      max_fdp <- min(c(fdr_max,max_fdp))
    }
  }
  
  gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
    geom_abline(slope = 1,intercept = 0,color="gray")+
    #rasterise(geom_line(), dpi = 300) + 
    geom_line()+
    xlab("FDR threshold")+
    ylab("Estimated FDP")+
    theme_bw()+
    #geom_segment(x=0.01,xend=0.01,y=0,yend=0.01,linewidth=0.3,color="blue",linetype=2)+
    theme_pubr(base_size = 12,border = TRUE)
  
  if(!is.null(color_mapping)){
    gg1 <- gg1 + scale_color_manual(values = color_mapping)
  }
  
  if(scale_xy){
    
    if(!is.na(max_x) || !is.na(max_y)){
      gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")
      if(!is.na(max_x)){
        gg1 <- gg1 + xlim(0,max_x) + scale_x_continuous(labels = scales::percent,limits =c(0,max_x))
      }else{
        gg1 <- gg1 + xlim(0,max_fdp)+ scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
      }
      if(!is.na(max_y)){
        gg1 <- gg1 + ylim(0,max_y) + scale_y_continuous(labels = scales::percent,limits =c(0,max_y))
      }else{
        gg1 <- gg1 + ylim(0,max_fdp) + scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))
      }
    }else{
      gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
        xlim(0,max_fdp)+
        ylim(0,max_fdp)+
        scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
        #scale_y_pct()+
        #scale_x_pct()+
        scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
    }
  }
  #theme(legend.position = "top",plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"))+    
  ## legend on the botton right
  ## no background color for the legend
  if(is.null(legend_position)){
    gg1 <- gg1 + theme(legend.position = c(0.65, 0.16),plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank(),legend.text=element_text(size=12),legend.title=element_text(size=12))    
  }else{
    gg1 <- gg1 + theme(legend.position = legend_position,plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank(),legend.text=element_text(size=12),legend.title=element_text(size=12))    
  }
  
  
  if(!is.null(fig_title)){
    gg1 <- gg1 + ggtitle(fig_title)
  }
  
  added_numbers <- NULL
  if(add_numbers){
    # add numbers on the top left size of the figure using annotation, text align to left
    # y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
    if(fdr_decimal_place==1){
      y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.1f%%",FDP001*100))
    }else if(fdr_decimal_place==2){
      y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.2f%%",FDP001*100))
    }else{
      y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
    }
    #n_t <- 
    if(abs(max_fdp-0.01)<=0.02){
      # text right align
      added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
      if(add_max_qvalue){
        added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
      }
      if(is.null(text_position)){
        gg1 <- gg1 + annotate("text", x = max_fdp*0.1, y = 0.9*max_fdp, label = added_numbers, color = "black", size = 3,hjust = 0)
      }else{
        gg1 <- gg1 + annotate("text", x = text_position[1], y = text_position[2], label = added_numbers, color = "black", size = 3,hjust = 0)
      }
      
    }else{
      added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
      if(add_max_qvalue){
        added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
      }
      if(is.null(text_position)){
        gg1 <- gg1 + annotate("text", x = 0.01*1.05, y = 0.9*max_fdp, label = added_numbers, color = "black", size = 3,hjust = 0)
      }else{
        gg1 <- gg1 + annotate("text", x = text_position[1], y = text_position[2], label = added_numbers, color = "black", size = 3,hjust = 0)
      }
    }
  }
  
  # library(ggpubr)
  # options(repr.plot.width = 12, repr.plot.height = 6)
  # gg <- ggarrange(gg1,gg2,ncol = 2,  common.legend = TRUE)
  # options(jupyter.plot_mimetypes = "image/png")
  if(return_data){
    return(list(gg=gg1,data=dat,added_numbers=added_numbers))
  }else{
    return(gg1)
  }
  
  #pdf("fdp_fdr_lib_qvalue_protein.pdf",width = 8,height = 4)
  #print(gg)
  #dev.off()
}

plot_fdp_fdr_multiple=function(dat,added_numbers=NULL,n_row_plots=2,n_col_plots=2,fdr_max=NULL,fig_title=NULL,scale_xy=TRUE,add_numbers=FALSE,r=1,fixed_fdr_max=FALSE,max_x=NA,max_y=NA,
                               color_mapping=NULL,
                               legend_position=NULL,
                               fdr_decimal_place=2,
                               return_data=FALSE) {
  
  max_fdp <- max(c(dat$FDP,dat$q_value))
  if(!is.null(fdr_max)){
    if(fixed_fdr_max){
      max_fdp <- fdr_max
    }else{
      max_fdp <- min(c(fdr_max,max_fdp))
    }
  }
  
  gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
    geom_abline(slope = 1,intercept = 0,color="gray")+
    #rasterise(geom_line(), dpi = 300) + 
    geom_line()+
    facet_wrap(.~tool,nrow=n_row_plots,ncol=n_col_plots)+
    xlab("FDR threshold")+
    ylab("Estimated FDP")+
    theme_bw()+
    #geom_segment(x=0.01,xend=0.01,y=0,yend=0.01,linewidth=0.3,color="blue",linetype=2)+
    theme_pubr(base_size = 12,border = TRUE)
  
  if(!is.null(color_mapping)){
    gg1 <- gg1 + scale_color_manual(values = color_mapping)
  }
  
  if(scale_xy){
    
    if(!is.na(max_x) || !is.na(max_y)){
      gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")
      if(!is.na(max_x)){
        gg1 <- gg1 + xlim(0,max_x) + scale_x_continuous(labels = scales::percent,limits =c(0,max_x))
      }else{
        gg1 <- gg1 + xlim(0,max_fdp)+ scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
      }
      if(!is.na(max_y)){
        gg1 <- gg1 + ylim(0,max_y) + scale_y_continuous(labels = scales::percent,limits =c(0,max_y))
      }else{
        gg1 <- gg1 + ylim(0,max_fdp) + scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))
      }
    }else{
      gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
        xlim(0,max_fdp)+
        ylim(0,max_fdp)+
        scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
        #scale_y_pct()+
        #scale_x_pct()+
        scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
    }
  }
  #theme(legend.position = "top",plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"))+    
  ## legend on the botton right
  ## no background color for the legend
  if(is.null(legend_position)){
    gg1 <- gg1 + theme(legend.position = c(0.65, 0.16),plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank())    
  }else{
    gg1 <- gg1 + theme(legend.position = legend_position,plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),legend.background = element_blank())    
  }
  
  
  if(!is.null(fig_title)){
    gg1 <- gg1 + ggtitle(fig_title)
  }
  
  added_numbers <- NULL
  if(add_numbers){
    
  }
  
  return(gg1)
}

plot_fdp_fdr_modified <- function(fdp_file="", # the FDP (in csv format) estimation file generated by FDRBench
                      fdr_max=NULL,
                      fig_title=NULL,
                      scale_xy=TRUE,
                      add_numbers=FALSE,
                      numbers_position=NULL,
                      numbers_font_size=11,
                      r=1,
                      fixed_fdr_max=FALSE,
                      max_x=NA,
                      max_y=NA,
                      color_mapping=NULL,
                      legend_position=c(0.7, 0.16),
                      legend_font_size=11,
                      fdr_decimal_place=2,
                      return_data=FALSE,
                      add_max_qvalue=FALSE) {
  if(is.null(color_mapping)){
    color_mapping <- c(
      "Paired method" = "#009E73", 
      "Sample method" = "#E69F00", 
      "Lower bound" = "#0072B2", 
      "Combined method" = "#D55E00",
      "Matched method"  = "#CC79A7"
      )
  }
  x <- read_csv(fdp_file)
  if("paired_fdp" %in% names(x)){
    if(r>=2){
      dat <- x %>% select(q_value,combined_fdp,paired_fdp,lower_bound_fdp) %>% distinct() %>% 
        rename(`Combined method`=combined_fdp,`Matched method`=paired_fdp,`Lower bound`=lower_bound_fdp) %>%
        gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
      dat$Method <- factor(dat$Method, levels = c("Combined method","Matched method","Lower bound"))
    }else{
      dat <- x %>% select(q_value,combined_fdp,paired_fdp,lower_bound_fdp) %>% distinct() %>% 
        rename(`Combined method`=combined_fdp,`Paired method`=paired_fdp,`Lower bound`=lower_bound_fdp) %>%
        gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
      dat$Method <- factor(dat$Method, levels = c("Combined method","Paired method","Lower bound"))
    }
  }else{
    dat <- x %>% select(q_value,combined_fdp,lower_bound_fdp) %>% distinct() %>% 
      rename(`Combined method`=combined_fdp,`Lower bound`=lower_bound_fdp) %>%
      gather(key = "Method",value = "FDP",-`q_value`) %>% select(q_value,FDP,Method)
    dat$Method <- factor(dat$Method, levels = c("Combined method","Lower bound"))
  }
  
  
  max_fdp <- max(c(dat$FDP,dat$q_value))
  if(!is.null(fdr_max)){
    if(fixed_fdr_max){
      max_fdp <- fdr_max
    }else{
      max_fdp <- min(c(fdr_max,max_fdp))
    }
  }
  
  # gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
  #   geom_abline(slope = 1,intercept = 0,color="gray")+
  #   geom_line()+
  #   xlab("FDR threshold")+
  #   ylab("Estimated FDP")+
  #   theme_bw()+
  #   theme_pubr(base_size = 12,border = TRUE)

  gg1 <- ggplot(dat,aes(x=q_value,y=FDP,color=Method)) + 
    geom_abline(slope = 1, intercept = 0, color = "gray") +
    geom_line(linewidth = 0.9) +
    xlab("FDR threshold") +
    ylab("Estimated FDP") +
    scale_color_manual(values = color_mapping) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(family = "Arial"),
      axis.title = element_text(family = "Arial"),
      legend.text = element_text(size = legend_font_size, family = "Arial"),
      legend.title = element_text(size = legend_font_size, family = "Arial"),
      plot.title = element_text(family = "Arial"),
      legend.position = "inside",
      legend.position.inside = legend_position,
      legend.background = element_blank(),
      plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      panel.border = element_blank(),
      axis.line.x.top = element_blank(),
      axis.line.y.right = element_blank()
    )
  
  if(!is.null(color_mapping)){
    gg1 <- gg1 + scale_color_manual(values = color_mapping)
  }
  
  if(scale_xy){
    
    if(!is.na(max_x) || !is.na(max_y)){
      gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")
      if(!is.na(max_x)){
        gg1 <- gg1 + xlim(0,max_x) + scale_x_continuous(labels = scales::percent,limits =c(0,max_x))
      }else{
        gg1 <- gg1 + xlim(0,max_fdp)+ scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
      }
      if(!is.na(max_y)){
        gg1 <- gg1 + ylim(0,max_y) + scale_y_continuous(labels = scales::percent,limits =c(0,max_y))
      }else{
        gg1 <- gg1 + ylim(0,max_fdp) + scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))
      }
    }else{
      gg1 <- gg1 + geom_vline(xintercept = 0.01,linetype=2,color="blue")+
        xlim(0,max_fdp)+
        ylim(0,max_fdp)+
        scale_y_continuous(labels = scales::percent,limits =c(0,max_fdp))+
        scale_x_continuous(labels = scales::percent,limits =c(0,max_fdp))
    }
  }
  
  # gg1 <- gg1 + theme(
  #   text = element_text(family = "Arial"),
  #   legend.position="inside",
  #   legend.position.inside = legend_position,
  #   legend.background = element_blank(),
  #   legend.text=element_text(size=legend_font_size),
  #   legend.title=element_text(size=legend_font_size),
  #   plot.margin = unit(2*c(0.1, 0.1, 0.1, 0.1),"inches"),
  #   axis.text.y = element_text(angle = 90, hjust = 0.5)
  #   )    
  
  
  if(!is.null(fig_title)){
    gg1 <- gg1 + ggtitle(fig_title)
  }
  
  added_numbers <- NULL
  # if(add_numbers){
  #   if(fdr_decimal_place==1){
  #     y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.1f%%",FDP001*100))
  #   }else if(fdr_decimal_place==2){
  #     y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.2f%%",FDP001*100))
  #   }else{
  #     y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
  #   }
  #   if(abs(max_fdp-0.01)<=0.02){
  #     added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
  #     if(add_max_qvalue){
  #       added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
  #     }
  #     if(is.null(numbers_position)){
  #       gg1 <- gg1 + annotate("text", x = max_fdp*0.1, y = 0.9*max_fdp, label = added_numbers, color = "black", size = numbers_font_size/.pt,hjust = 0)
  #     }else{
  #       gg1 <- gg1 + annotate("text", x = numbers_position[1], y = numbers_position[2], label = added_numbers, color = "black", size = numbers_font_size/.pt,hjust = 0)
  #     }
      
  #   }else{
  #     added_numbers <- paste("Total discoveries:",nrow(x %>% filter(q_value<=0.01)),"\n",paste(y$Method,y$ratio,sep=":",collapse = "\n"),sep="")
  #     if(add_max_qvalue){
  #       added_numbers <- paste(added_numbers,"\n","Max q-value:",sprintf("%.2e",max(x$q_value)),sep="")
  #     }
  #     if(is.null(numbers_position)){
  #       gg1 <- gg1 + annotate("text", x = 0.01*1.05, y = 0.9*max_fdp, label = added_numbers, color = "black", size = numbers_font_size/.pt,hjust = 0)
  #     }else{
  #       gg1 <- gg1 + annotate("text", x = numbers_position[1], y = numbers_position[2], label = added_numbers, color = "black", size = numbers_font_size/.pt,hjust = 0)
  #     }
  #   }
  # }
  
  if(add_numbers){
    if(fdr_decimal_place==1){
      y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.1f%%",FDP001*100))
    } else if(fdr_decimal_place==2){
      y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.2f%%",FDP001*100))
    } else {
      y <- dat %>% filter(q_value<=0.01) %>% group_by(Method) %>% arrange(desc(q_value)) %>% filter(row_number()==1) %>% summarise(FDP001=max(FDP)) %>% mutate(ratio=sprintf("%.4f%%",FDP001*100))
    }
    
    added_numbers <- paste("Total discoveries:", nrow(x %>% filter(q_value<=0.01)),
                           "\n", paste(y$Method, y$ratio, sep=":", collapse = "\n"), sep="")
    if(add_max_qvalue){
      added_numbers <- paste(added_numbers, "\n", "Max q-value:", sprintf("%.2e",max(x$q_value)), sep="")
    }
    
    if(is.null(numbers_position)){
      xpos <- if(abs(max_fdp-0.01)<=0.02) max_fdp*0.1 else 0.01*1.05
      gg1 <- gg1 + annotate("text",
                            x = xpos,
                            y = 0.9*max_fdp,
                            label = added_numbers,
                            color = "black",
                            size = numbers_font_size/.pt,
                            hjust = 0,
                            family = "Arial")
    } else {
      gg1 <- gg1 + annotate("text",
                            x = numbers_position[1],
                            y = numbers_position[2],
                            label = added_numbers,
                            color = "black",
                            size = numbers_font_size/.pt,
                            hjust = 0,
                            family = "Arial")
    }
  }
  
  
  
  if(return_data){
    return(list(gg=gg1,data=dat,added_numbers=added_numbers))
  }else{
    return(gg1)
  }
  
}

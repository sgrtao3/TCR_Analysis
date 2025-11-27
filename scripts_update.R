library(tidyverse)
library(readr)
library(stringr)
library(janitor)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(randomForest)
library(scales)
library(broom)
library(patchwork)
library(umap)
library(Rtsne)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)

tcr_dir <- "G:/Ranran/syn61987835/download/TCR"
metadata_path <- "G:/Ranran/syn61987835/metadata/metadata.tsv"
out_dir <- "G:/Ranran/syn61987835/results"; dir.create(out_dir, showWarnings = FALSE)

# read metadata ----
meta <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  clean_names() %>%
  # group：Healthy/Background -> Healthy，keep other diseases
  mutate(
    group = if_else(str_detect(disease, regex("Healthy", ignore_case = TRUE)),
                    "Healthy", disease),
    age = suppressWarnings(as.numeric(age))
  )

# read all the TCR features and calculate diversity metrics ----
# helper：metrics for a single sample
read_one_and_summarise <- function(fpath) {
  # Extract participant label from filename (e.g., part_table_bfi-0000234 -> BFI-0000234)
  bn <- basename(fpath)
  m <- stringr::str_match(bn, "(?i)bfi-(\\d+)")
  id_num <- m[, 2]
  participant_label <- if (!is.na(id_num)) {
    paste0("BFI-", stringr::str_pad(id_num, width = 7, pad = "0"))
  } else {
    "BFI-NA"
  }
  # Read .tsv.gz file and clean column names
  df <- read_tsv(fpath, show_col_types = FALSE) %>% clean_names()
  # Check for essential column
  if (!("productive" %in% names(df))) {
    warning("Missing 'productive' column: ", fpath)
    return(NULL)
  }
  #cdr3_candidates <- c("junction_aa", "cdr3_aa", "cdr3")
  cdr3_candidates <- "cdr3_aa"
  cdr3_col <- cdr3_candidates[cdr3_candidates %in% names(df)][1] %||% NA_character_
  # Choose which column to define the clone identity
  # Priority: CDR3 amino acid > clone_id > sequence_id
  clone_id_col <- NA_character_
  # Normalize productive values (TRUE/FALSE, T/F, 1/0)
  prod_flag <- tolower(as.character(df$productive)) %in% c("t", "true", "1", "yes")
  df2 <- df %>% filter(prod_flag)
  # Try to use CDR3 if available and non-empty
  if (!is.na(cdr3_col)) {
    nonempty_cdr3 <- df2 %>% filter(!is.na(.data[[cdr3_col]]), .data[[cdr3_col]] != "")
    if (nrow(nonempty_cdr3) > 0) clone_id_col <- cdr3_col
  }
  # If no valid CDR3, use clone_id
  if (is.na(clone_id_col) && "clone_id" %in% names(df2)) {
    nonempty_clone <- df2 %>% filter(!is.na(clone_id), clone_id != "")
    if (nrow(nonempty_clone) > 0) clone_id_col <- "clone_id"
  }
  # If still missing, fall back to sequence_id
  if (is.na(clone_id_col) && "sequence_id" %in% names(df2)) {
    nonempty_seq <- df2 %>% filter(!is.na(sequence_id), sequence_id != "")
    if (nrow(nonempty_seq) > 0) clone_id_col <- "sequence_id"
  }
  # If no usable identifier found, return empty record
  if (is.na(clone_id_col)) {
    warning("No usable clone identifier (CDR3 / clone_id / sequence_id) found: ", fpath)
    return(tibble(
      participant_label = participant_label,
      n_sequences = 0, n_clones = 0,
      shannon = NA_real_, simpson = NA_real_,
      pielou = NA_real_, cdr3_len_mean = NA_real_,
      id_column_used = NA_character_
    ))
  }
  message("→ [", participant_label, "] Using clone identifier column: ", clone_id_col)
  # Define abundance weights (counts)
  w <- if ("duplicate_count" %in% names(df2)) {
    as.numeric(df2$duplicate_count)
  } else if ("count" %in% names(df2)) {
    as.numeric(df2$count)
  } else {
    rep(1, nrow(df2))
  }
  df2$w <- replace_na(w, 1)
  # Compute weighted mean CDR3 length if CDR3 column is available
  if (!is.na(cdr3_col)) {
    df2 <- df2 %>%
      filter(!is.na(.data[[clone_id_col]]), .data[[clone_id_col]] != "", w > 0) %>%
      mutate(cdr3_len = nchar(.data[[cdr3_col]]))
    cdr3_len_mean <- if (sum(!is.na(df2$cdr3_len)) > 0) {
      tmp <- df2 %>% group_by(.data[[clone_id_col]]) %>%
        summarise(len = first(cdr3_len), w = sum(w), .groups = "drop")
      weighted.mean(tmp$len, tmp$w)
    } else NA_real_
  } else {
    df2 <- df2 %>%
      filter(!is.na(.data[[clone_id_col]]), .data[[clone_id_col]] != "", w > 0)
    cdr3_len_mean <- NA_real_
  }
  # If file is empty after filtering, return NA metrics
  if (nrow(df2) == 0) {
    return(tibble(
      participant_label = participant_label,
      n_sequences = 0, n_clones = 0,
      shannon = NA_real_, simpson = NA_real_,
      pielou = NA_real_, cdr3_len_mean = cdr3_len_mean,
      id_column_used = clone_id_col
    ))
  }
  # Aggregate abundances by clone
  abund <- df2 %>%
    group_by(.data[[clone_id_col]]) %>%
    summarise(abundance = sum(w, na.rm = TRUE), .groups = "drop")
  # Diversity metrics
  S <- nrow(abund)
  shannon <- if (S > 0) vegan::diversity(abund$abundance, index = "shannon") else NA_real_
  simpson <- if (S > 0) vegan::diversity(abund$abundance, index = "simpson") else NA_real_
  Hmax <- if (S > 0) log(S) else NA_real_
  pielou <- if (!is.na(Hmax) && Hmax > 0) shannon / Hmax else NA_real_
  # Return sample-level summary
  tibble(
    participant_label = participant_label,
    n_sequences = nrow(df2),
    n_clones = S,
    shannon = shannon,
    simpson = simpson,
    pielou = pielou,
    cdr3_len_mean = cdr3_len_mean,
    id_column_used = clone_id_col
  )
}
# batch calculate all the samples
tcr_files <- list.files(tcr_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
message("TCR file counts：", length(tcr_files))
# comment to avoid a second running
#feat_list <- lapply(tcr_files, function(f) {              
#  tryCatch(read_one_and_summarise(f),
#           error = function(e) { warning(e$message); return(NULL) })
#})
#features <- bind_rows(feat_list)

# save
#write_csv(features, file.path(out_dir, "tcr_sample_features.csv"))
features <- read_csv(file.path(out_dir, "tcr_sample_features.csv"))


# combine metadata ----
meta_unique <- meta %>%
  distinct(participant_label, .keep_all = TRUE)
dat <- features %>%
  left_join(
    meta_unique %>% select(participant_label, group, age, sex),
    by = "participant_label"
  ) %>%
  mutate(
    group_binary = if_else(group == "Healthy", "Healthy", "Diseased"),
    group_binary = factor(group_binary, levels = c("Healthy", "Diseased"))
  )
#write_csv(dat, file.path(out_dir, "tcr_features_with_metadata.csv"))
dat <- read_csv(file.path(out_dir, "tcr_features_with_metadata.csv"))


# Create an output dir
outdir <- "EDA_plots"
dir.create(outdir, showWarnings = FALSE)

# Clonality via cumulative freq of top-N clones ----
to_logical_true <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) return(tolower(x) %in% c("true","t","1","yes","y"))
  rep(FALSE, length(x))
}
# single-file function
topN_cumfreq_one_debug <- function(tsv_gz, N = 10,
                                   clone_col = "clone_id",
                                   productive_col = "productive",
                                   fallback_clone_cols = c("cdr3_aa","junction_aa")) {
  message("=== File: ", tsv_gz, " ===")
  df <- suppressMessages(readr::read_tsv(tsv_gz, col_types = cols(.default = "c")))
  message("Columns: ", paste(names(df), collapse = ", "))
  n_raw <- nrow(df); message("n_raw: ", n_raw)
  # Handle productive filter (robust)
  if (!(productive_col %in% names(df))) {
    message("[info] No '", productive_col, "' column; skip productive filtering.")
  } else {
    keep <- to_logical_true(df[[productive_col]])
    message("table(productive raw): ", paste(capture.output(print(table(df[[productive_col]], useNA="ifany"))), collapse=" | "))
    message("sum(keep): ", sum(keep, na.rm = TRUE))
    df <- df[keep %in% TRUE, , drop = FALSE]
  }
  n_keep <- nrow(df); message("n_keep after productive filter: ", n_keep)
  # Choose clone identifier column
  use_col <- clone_col
  if (!(use_col %in% names(df))) {
    message("[warn] '", clone_col, "' missing. Trying fallbacks: ", paste(fallback_clone_cols, collapse = ", "))
    hit <- fallback_clone_cols[fallback_clone_cols %in% names(df)]
    if (length(hit)) {
      # Pick first non-empty fallback
      for (c in hit) {
        if (!all(is.na(df[[c]])) && !all(trimws(df[[c]]) == "")) { use_col <- c; break }
      }
      if (!(use_col %in% names(df))) {
        message("[fail] No usable clone identifier column found. Returning NA.")
        return(tibble(participant_label = "BFI-NA", topN = N, cumfreq = NA_real_,
                      n_raw = n_raw, n_keep = n_keep, id_col = NA_character_, reason = "no_clone_identifier"))
      }
    } else {
      message("[fail] No fallback columns present. Returning NA.")
      return(tibble(participant_label = "BFI-NA", topN = N, cumfreq = NA_real_,
                    n_raw = n_raw, n_keep = n_keep, id_col = NA_character_, reason = "no_clone_identifier"))
    }
  }
  message("Using clone identifier column: ", use_col)
  if (n_keep == 0) {
    message("[fail] No rows after filtering; returning NA.")
    return(tibble(participant_label = "BFI-NA", topN = N, cumfreq = NA_real_,
                  n_raw = n_raw, n_keep = n_keep, id_col = use_col, reason = "zero_after_filter"))
  }
  cs <- df %>% count(.data[[use_col]], name = "n") %>% arrange(desc(n))
  total <- sum(cs$n)
  message("Unique clones: ", nrow(cs), " | total sequences after filter: ", total)
  if (total == 0) {
    message("[fail] total==0; returning NA.")
    return(tibble(participant_label = "BFI-NA", topN = N, cumfreq = NA_real_,
                  n_raw = n_raw, n_keep = n_keep, id_col = use_col, reason = "total_zero"))
  }
  topN_tbl <- head(cs, N)
  message("Top ", N, " clone sizes: ", paste(topN_tbl$n, collapse = ", "))
  cumfreq <- sum(topN_tbl$n) / total
  message("cumfreq = sum(topN)/total = ", sum(topN_tbl$n), " / ", total, " = ", signif(cumfreq, 4))
  # Participant label from filename like part_table_bfi-0000234.tsv.gz
  bn <- basename(tsv_gz)
  m  <- stringr::str_match(bn, "(?i)bfi-(\\d+)")
  id <- ifelse(!is.na(m[,2]),
               paste0("BFI-", stringr::str_pad(m[,2], width = 7, side = "left", pad = "0")),
               "BFI-NA")
  message("Participant inferred: ", id)
  tibble(participant_label = id, topN = N, cumfreq = cumfreq,
         n_raw = n_raw, n_keep = n_keep, id_col = use_col, reason = "ok")
}
# test one file
#f <- "G:/Ranran/syn61987835/download/TCR/part_table_bfi-0000234.tsv.gz"
#res1 <- topN_cumfreq_one_debug(f, N = 10)
#res1
# comment to avoid running a second time
#top10_df <- purrr::map_dfr(tcr_files, ~ topN_cumfreq_one_debug(.x, N = 10)) 
#table(top10_df$reason, useNA = "ifany")
#summary(top10_df$cumfreq)
# Save it to avoid recompute
#readr::write_tsv(top10_df, file.path(outdir, "top10_cumfreq.tsv"))
# Suppose loaded it back:
top10_df <- readr::read_tsv(file.path(outdir, "top10_cumfreq.tsv"))
# Join to metadata
top10_dat <- top10_df %>%
  left_join(dat %>% select(participant_label, group_binary, age, shannon), by = "participant_label")


# merge meta data, features, and top10_dat ----
meta_use <- meta_unique %>%
  select(participant_label,disease, disease_subtype,age,sex,ancestry,specimen_time_point,study_name,available_gene_loci)
# combine top10_dat and meta
final_df <- top10_dat %>%
  left_join(features, by = "participant_label") %>%  
  left_join(meta_use, by = "participant_label") %>%   
  mutate(
    sex = factor(sex),
    ancestry = factor(ancestry),
    disease = factor(disease),
    disease_subtype = factor(disease_subtype),
    group_binary = factor(group_binary, levels = c("Healthy", "Diseased"))
  )
# Save it to avoid recompute
#readr::write_tsv(final_df, file.path(outdir, "final_df_features.tsv"))
final_df <- readr::read_tsv(file.path(outdir, "final_df_features.tsv"))

# data overview ----
df1 <- final_df %>%
  count(group_binary)
p1 <- ggplot(df1, aes(x = "", y = n, fill = group_binary)) +
  geom_col(width = 1) +
  coord_polar("y") +
  labs(title = "Healthy vs Diseased sample distribution") +
  theme_void() +
  theme(legend.title = element_blank())
p1
#ggsave(file.path(outdir, "pie_healthy_vs_diseased.png"), p1, width = 5, height = 5)
df2 <- final_df %>%
  filter(disease %in% c("HIV","Covid19","Influenza","Lupus","T1D")) %>%
  count(disease)
p2 <- ggplot(df2, aes(x = "", y = n, fill = disease)) +
  geom_col(width = 1) +
  coord_polar("y") +
  labs(title = "Distribution of Disease Categories") +
  theme_void() +
  theme(legend.title = element_blank())
p2
#ggsave(file.path(outdir,"pie_disease.png"), p2, width = 6, height = 6)
p_age <- ggplot(final_df, aes(age.x)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = .7) +
  #geom_density(color = "red", linewidth = 1) +
  theme_bw() +
  labs(title = "Age distribution", x = "Age", y = "Count")
p_age
#ggsave(file.path(outdir,"age_distribution.png"), p_age, width = 7, height = 5)
df_sex <- final_df %>%
  filter(sex %in% c("M","F")) %>%
  count(sex)
p_sex <- ggplot(df_sex, aes(x = "", y = n, fill = sex)) +
  geom_col(width = 1) +
  coord_polar("y") +
  labs(title = "Sex distribution") +
  theme_void()
p_sex
#ggsave(file.path(outdir,"sex_distribution.png"), p_sex, width = 5, height = 5)
df_tp <- final_df %>%
  filter(!is.na(specimen_time_point)) %>%
  mutate(tp_num = readr::parse_number(specimen_time_point)) %>%
  count(specimen_time_point, tp_num) %>%
  arrange(tp_num) %>%
  mutate(specimen_time_point = factor(specimen_time_point, levels = unique(specimen_time_point)))
p_tp <- ggplot(df_tp, aes(x = specimen_time_point, y = n)) +
  geom_col(fill = "gray60") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Specimen time points distribution",
    x = "Time point",
    y = "Number of samples"
  )
p_tp
#ggsave(file.path(outdir,"timepoint_distribution.png"), p_tp, width = 8, height = 4)
df_ancestry <- final_df %>%
  count(ancestry)
p_ancestry <- ggplot(df_ancestry, aes(x = "", y = n, fill = ancestry)) +
  geom_col(width = 1) +
  coord_polar("y") +
  labs(title = "Distribution of Ancestry Categories") +
  theme_void() +
  theme(legend.title = element_blank())
p_ancestry
#ggsave(file.path(outdir,"pie_ancestry.png"), p_ancestry, width = 6, height = 6)

# Age effect analysis ----
top10_dat$group_binary <- factor(
  top10_dat$group_binary,
  levels = c("Healthy", "Diseased")
)
cols <- c("Healthy" = "#1F77B4",   # blue
          "Diseased" = "#D62728")  # red
p1 <- ggplot(top10_dat, aes(age, cumfreq, color = group_binary)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  scale_color_manual(values = cols) +
  labs(title = "Age vs TCR Clonality",
       y = "Top10 cumulative frequency",
       x = "Age") +
  theme_bw(base_size = 16)
p1
#ggsave(file.path(outdir, "Age_vs_TCR_clonality.png"), p1, width = 7.5, height = 5)
p2 <- ggplot(top10_dat, aes(age, shannon, color = group_binary)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  scale_color_manual(values = cols) +
  labs(title = "Age vs TCR Diversity",
       y = "Shannon diversity",
       x = "Age") +
  theme_bw(base_size = 16)
p2
#ggsave(file.path(outdir, "Age_vs_TCR_diversity.png"), p2, width = 7.5, height = 5)
lm(cumfreq ~ age + group_binary, data = top10_dat)
lm(shannon ~ age + group_binary, data = top10_dat)



# TRBV usage ----
# tcr_files: character vector of all TCR .tsv.gz files
# e.g. tcr_files <- list.files("G:/Ranran/syn61987835/download/TCR", pattern="\\.tsv\\.gz$", full.names=TRUE)
# Extract TRBV gene from v_call (first gene, no allele)
extract_trbv <- function(v_call_vec) {
  # v_call may contain multiple genes separated by , or ;
  part1 <- str_split(v_call_vec, "[,;]", simplify = TRUE)[, 1]
  # remove trailing description after space
  part1 <- str_split(part1, "\\s+", simplify = TRUE)[, 1]
  # remove allele suffix after *
  gene  <- str_split(part1, "\\*", simplify = TRUE)[, 1]
  gene
}
# Read one TCR file and compute TRBV usage for this participant
compute_trbv_usage_one <- function(tsv_gz,
                                   productive_col = "productive",
                                   locus_col = "locus",
                                   v_col = "v_call") {
  df <- suppressMessages(readr::read_tsv(tsv_gz, col_types = cols(.default = "c")))
  # Filter productive sequences if column exists
  if (productive_col %in% names(df)) {
    keep <- to_logical_true(df[[productive_col]])
    df <- df[keep %in% TRUE, , drop = FALSE]
  }
  if (locus_col %in% names(df)) {
    df <- df %>% filter(.data[[locus_col]] %in% c("TRB","TCRB", "TRB_beta", "TRBchain"))
  }
  # If no v_call or no rows, return NA
  if (!(v_col %in% names(df)) || nrow(df) == 0) {
    # infer participant_label from filename anyway
    bn <- basename(tsv_gz)
    m  <- str_match(bn, "(?i)bfi-(\\d+)")
    pid <- ifelse(!is.na(m[,2]),
                  paste0("BFI-", str_pad(m[,2], width = 7, side = "left", pad = "0")),
                  "BFI-NA")
    return(tibble(participant_label = pid,TRBV = NA_character_,n = NA_integer_,freq = NA_real_))
  }
  # Extract TRBV gene name
  df <- df %>%
    mutate(TRBV = extract_trbv(.data[[v_col]])) %>%
    filter(!is.na(TRBV), TRBV != "")
  
  if (nrow(df) == 0) {
    bn <- basename(tsv_gz)
    m  <- str_match(bn, "(?i)bfi-(\\d+)")
    pid <- ifelse(!is.na(m[,2]),
                  paste0("BFI-", str_pad(m[,2], width = 7, side = "left", pad = "0")),
                  "BFI-NA")
    return(tibble(participant_label = pid, TRBV = NA_character_,n = NA_integer_, freq = NA_real_))
  }
  
  counts <- df %>%
    count(TRBV, name = "n") %>%
    arrange(desc(n))
  
  total_n <- sum(counts$n)
  counts <- counts %>%
    mutate(freq = n / total_n)
  
  # infer participant_label from filename
  bn <- basename(tsv_gz)
  m  <- str_match(bn, "(?i)bfi-(\\d+)")
  pid <- ifelse(!is.na(m[,2]),
                paste0("BFI-", str_pad(m[,2], width = 7, side = "left", pad = "0")),
                "BFI-NA")
  
  counts %>%
    mutate(participant_label = pid) %>%
    select(participant_label, TRBV, n, freq)
}
# test one file
#f <- "G:/Ranran/syn61987835/download/TCR/part_table_bfi-0000234.tsv.gz"
#res2 <- compute_trbv_usage_one(f)
#res2
# Run for all TCR files
#trbv_usage <- purrr::map_dfr(tcr_files, compute_trbv_usage_one)  # comment to avoid a second running
# remove rows that are all NA (no TRBV info)
#trbv_usage <- trbv_usage %>% filter(!is.na(TRBV))
#readr::write_tsv(trbv_usage, file.path(outdir, "trbv_usage.tsv"))
trbv_usage <- read_tsv(file.path(outdir, "trbv_usage.tsv"))


# confounder exploration ----
# check missing value 
final_df <- final_df %>%
  mutate(
    specimen_time_point_num = readr::parse_number(specimen_time_point)
  )
df_clono_na_summary <- final_df %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(cols = everything(),names_to = "variable", values_to = "na_count") %>%
  mutate(
    total = nrow(final_df),
    na_percent = na_count / total * 100
  ) %>%
  arrange(desc(na_count))
df_clono_na_summary
#write.csv(df_clono_na_summary, file.path(outdir, "df_clono_na_summarys.csv"))
na_df <- final_df %>% 
  summarise(
    specimen_time_point = sum(is.na(specimen_time_point)),
    ancestry = sum(is.na(ancestry)),
    sex = sum(is.na(sex)),
    age = sum(is.na(age.x))
  ) %>% 
  pivot_longer(everything(), names_to="variable", values_to="na_count") %>%
  mutate(total = nrow(final_df),
         na_percent = na_count / total * 100)
p <- ggplot(na_df, aes(x = reorder(variable, na_percent), y = na_percent)) +
  geom_col(fill = "skyblue3") +
  geom_text(aes(label = sprintf("%.1f%%", na_percent)), 
            vjust = -0.5, size=5) +
  theme_bw(base_size = 16) +
  labs(title = "Missing Value Percentage", 
       x = NULL, y = "% Missing")+
  theme(
    plot.title = element_text(size = 15, face = "bold"),   
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 5)                  
  )
p
#ggsave("Missing Value Percentage.png", p, width = 9, height = 6)

# do not need to consider specimen_time_point_num and ancestry as potential confounders due to many NA
run_univariate <- function(df, yvar){
  conf_vars <- c("age.x", "sex")
  results <- lapply(conf_vars, function(v){
    formula <- as.formula(paste0(yvar, " ~ ", v))
    m <- lm(formula, data = df)
    tidy_m <- broom::tidy(m) %>% filter(term != "(Intercept)")
    tidy_m %>%
      mutate(variable = v) %>%
      select(variable, term, estimate, p.value)
  }) %>% bind_rows()
  return(results)
}
run_multivariate <- function(df, yvar){
  formula <- as.formula(
    paste0(yvar, " ~ age.x + sex ")
  )
  m <- lm(formula, data = df)
  tidy(m)
}
df_clono <- final_df %>%
  transmute(
    cumfreq,
    shannon,
    age.x,
    sex = factor(sex),
    #ancestry = factor(ancestry),
    #specimen_time_point_num = readr::parse_number(specimen_time_point),
    #study_name = factor(study_name)
  ) %>%
  drop_na() 
# Identify valid predictors (level >= 2)
valid_pred <- sapply(df_clono, function(x){
  if (is.character(x) || is.factor(x)) {
    length(unique(x)) >= 2
  } else TRUE
})
df_clono <- df_clono[, valid_pred]
# clono Univariate
clono_uni <- run_univariate(df_clono, "cumfreq")
#readr::write_tsv(clono_uni, "clonality_univariate.tsv")
# clono Multivariate
clono_multi <- run_multivariate(df_clono, "cumfreq")
#readr::write_tsv(clono_multi, "clonality_multivariate.tsv")

# diversity Univariate
div_uni <- run_univariate(df_clono, "shannon")
#readr::write_tsv(div_uni, "diversity_univariate.tsv")
# diversity Multivariate
div_multi <- run_multivariate(df_clono, "shannon")
#readr::write_tsv(div_multi, "diversity_multivariate.tsv")

# Prepare forest plot data
df_forest_clono <- clono_multi %>%
  filter(term != "(Intercept)") %>% 
  mutate(
    term = factor(term, levels = rev(term)),
    signif = p.value < 0.05,
    CI_low = estimate - 1.96 * std.error,
    CI_high = estimate + 1.96 * std.error
  )
# clonality
p_forest_clono <- ggplot(df_forest_clono,
                         aes(x = estimate, y = term)) +
  geom_point(aes(color = signif), size = 3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = signif),
                 height = 0.18, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#d62728")) +
  labs(
    title = "Forest Plot – Multivariable Effects on Clonality",
    x = "Estimate (± 95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )
p_forest_clono
#ggsave("forest_clonality_multivariate.png",p_forest_clono,width = 6, height = 4, dpi = 300)
# diversity
df_forest_div <- div_multi %>%
  filter(term != "(Intercept)") %>% 
  mutate(
    term = factor(term, levels = rev(term)),
    signif = p.value < 0.05,
    CI_low = estimate - 1.96 * std.error,
    CI_high = estimate + 1.96 * std.error
  )
p_forest_div <- ggplot(df_forest_div,
                       aes(x = estimate, y = term)) +
  geom_point(aes(color = signif), size = 3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = signif),
                 height = 0.18, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#d62728")) +
  labs(
    title = "Forest Plot – Multivariable Effects on Diversity",
    x = "Estimate (± 95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )
p_forest_div
#ggsave("forest_diversity_multivariate.png",p_forest_div,width = 6, height = 4, dpi = 300)

# age-adgusted analysis ----
get_residual <- function(df, response, covar = "age.x") {
  f <- as.formula(paste0(response, " ~ ", covar))
  model <- lm(f, data = df)
  resid <- resid(model)
  return(resid)
}
# Age-adjusted clonality violin plot
df_adj <- final_df %>%
  filter(!is.na(cumfreq), !is.na(age.x)) %>%
  mutate(cumfreq_adj = get_residual(., "cumfreq", "age.x"))
p1 <- ggplot(df_adj, aes(disease, cumfreq_adj, fill = disease)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  stat_summary(fun = median, geom = "point", size = 2.5, color = "black") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    axis.title.y = element_text(size = 18),   
    plot.title = element_text(size = 20),   
    plot.subtitle = element_text(size = 16),  
    legend.position = "none"
  )+
  labs(title = "Age-adjusted TCR clonality",
       subtitle = "Residuals of (cumfreq ~ age)",
       y = "Clonality (age-adjusted residual)")
p1
#ggsave(file.path(outdir, "Age_adjusted_TCR_clonality_across_diseases.png"), p1, width = 7, height = 5)

# Age-adjusted diversity violin plot
df_adj <- df_adj %>%
  mutate(shannon_adj = get_residual(., "shannon", "age.x"))
p2 <- ggplot(df_adj, aes(disease, shannon_adj, fill = disease)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  stat_summary(fun = median, geom = "point", size = 2.5, color = "black") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    axis.title.y = element_text(size = 18),   
    plot.title = element_text(size = 20),   
    plot.subtitle = element_text(size = 16),  
    legend.position = "none"
  )+
  labs(title = "Age-adjusted TCR diversity",
       subtitle = "Residuals of (Shannon ~ age)",
       y = "Shannon (age-adjusted residual)")
p2
#ggsave(file.path(outdir, "Age_adjusted_TCR_diversity_across_diseases.png"), p2, width = 7, height = 5)


# TRBV heatmap
df_trbv <- trbv_usage %>%
  left_join(dat %>% select(participant_label, group, age),
            by = "participant_label") %>%
  filter(!is.na(group), !is.na(age), !is.na(freq))
# Age-adjust each TRBV gene using linear regression
adjust_one_trbv <- function(df_sub){
  fit <- lm(freq ~ age, data = df_sub)
  df_sub$freq_adj <- residuals(fit)
  return(df_sub)
}
df_adj <- df_trbv %>%
  group_by(TRBV) %>%
  group_modify(~ adjust_one_trbv(.x)) %>%
  ungroup()
# Compute group-mean age-adjusted TRBV usage
trbv_group_adj <- df_adj %>%
  group_by(group, TRBV) %>%
  summarize(freq_adj = mean(freq_adj, na.rm = TRUE), .groups = "drop")
# Pick top TRBVs for readability (based on absolute adjusted usage)
top_trbv <- trbv_group_adj %>%
  group_by(TRBV) %>%
  summarize(mean_abs = mean(abs(freq_adj), na.rm = TRUE)) %>%
  arrange(desc(mean_abs)) %>%
  slice_head(n = 30)   # keep top 30 informative TRBVs
trbv_group_filt <- trbv_group_adj %>%
  filter(TRBV %in% top_trbv$TRBV)
# heatmap
p_trbv_adj <- ggplot(trbv_group_filt, aes(x = TRBV, y = group, fill = freq_adj)) +
  geom_tile(color = "grey85") +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0
  ) +
  labs(
    title = "Age-adjusted TRBV usage across disease groups",
    x = "TRBV gene",
    y = NULL,
    fill = "Adj freq"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid = element_blank()
  )
p_trbv_adj
#ggsave(file.path(outdir, "TRBV_usage_age_adjusted_heatmap.png"),p_trbv_adj, width = 9, height = 4.8)

# Random forest to predict disease group
# numeric feature columns only
trbv_wide <- trbv_usage %>%
  filter(!is.na(TRBV)) %>%
  select(participant_label, TRBV, freq) %>%
  # if same (participant, TRBV) appears multiple times, average them
  group_by(participant_label, TRBV) %>%
  summarise(freq = mean(freq, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = TRBV,
    values_from = freq,
    values_fill = 0
  )
feat_core <- features %>%
  select(participant_label, shannon, simpson, pielou, cdr3_len_mean)
meta_use <- meta %>%
  distinct(participant_label, .keep_all = TRUE) %>%
  mutate(
    group_binary = if_else(disease == "Healthy/Background", "Healthy", "Diseased"),
    group_binary = factor(group_binary, levels = c("Healthy","Diseased")),
    sex = factor(sex)
  ) %>%
  select(participant_label, group_binary, disease, disease_subtype, age, sex, ancestry)
ana_df <- feat_core %>%
  inner_join(trbv_wide, by = "participant_label") %>%
  inner_join(meta_use,  by = "participant_label") %>%
  # remove rows with missing group
  filter(!is.na(group_binary))
rf_features <- ana_df %>%
  select(-participant_label,
         -group_binary,
         -disease_subtype,
         -sex,
         -ancestry,
         -age,
         -disease)   
# response variable (disease groups)
y_rf <- as.factor(ana_df$disease)
# convert to matrix
X_rf <- rf_features
# replace NA with 0
X_rf[is.na(X_rf)] <- 0
# Train Random Forest classifier
set.seed(2025)
rf_fit <- randomForest(
  x = X_rf,
  y = y_rf,
  ntree = 1000,
  importance = TRUE
)
print(rf_fit)   # confusion matrix, OOB error
sink("rf_output.txt")
print(rf_fit)
sink()
# Extract feature importance
imp <- importance(rf_fit, type = 2)
imp_df <- data.frame(
  feature = rownames(imp),
  MeanDecreaseGini = imp[, 1],
  row.names = NULL
) %>%
  arrange(desc(MeanDecreaseGini))
# top 20 features
top_imp <- imp_df %>% slice_head(n = 20)
# Plot feature importance
p_imp <- ggplot(top_imp,
                aes(x = reorder(feature, MeanDecreaseGini),
                    y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Random Forest Feature Importance",
    subtitle = "Top 20 predictors for distinguishing disease groups",
    x = NULL, 
    y = "Mean Decrease Gini"
  ) +
  theme_bw(base_size = 13)
p_imp
#ggsave(file.path(outdir, "RF_feature_importance_disease_groups.png"), p_imp, width = 7.5, height = 5)

# TRBV violin plot
# get all TRBV columns
trbv_cols <- grep("^TRBV", colnames(ana_df), value = TRUE)
# create empty list
age_adj_trbv <- ana_df[, c("participant_label", "disease")]
trbv_cols <- colnames(ana_df)[grepl("^TRBV", colnames(ana_df))]
for (g in trbv_cols) {
  # build model frame that keeps track of NA rows
  df_tmp <- model.frame(
    formula = ana_df[[g]] ~ ana_df$age
  )
  # lm will drop rows with NA inside df_tmp
  fit <- lm(df_tmp)
  # residuals for rows used in the model
  res_fit <- residuals(fit)
  # initialize full residual vector (all NA first)
  res_full <- rep(NA, nrow(ana_df))
  # which rows in ana_df correspond to complete cases?
  complete_idx <- complete.cases(df_tmp)
  # fill residuals back into their original positions
  res_full[complete_idx] <- res_fit
  # final step: replace NA with 0 (or keep NA if you prefer)
  res_full[is.na(res_full)] <- 0
  # store in output table
  age_adj_trbv[[g]] <- res_full
}
dim(age_adj_trbv)
plot_trbv_violin <- function(df, gene) {
  p <- ggplot(df, aes(x = disease, y = .data[[gene]], fill = disease)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", size = 2, color = "black") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none") +
    labs(
      title = paste0("Age-adjusted ", gene, " usage across diseases"),
      y = paste0(gene, " (age-adjusted residual)")
    )
  print(p)
}
TRBV14 <- plot_trbv_violin(age_adj_trbv, "TRBV14")
TRBV29_1 <- plot_trbv_violin(age_adj_trbv, "TRBV29-1")
TRBV27 <- plot_trbv_violin(age_adj_trbv, "TRBV27")
#ggsave(file.path(outdir, "TRBV14.png"),TRBV14, width = 7.5, height = 5)
#ggsave(file.path(outdir, "TRBV29_1.png"),TRBV29_1, width = 7.5, height = 5)
#ggsave(file.path(outdir, "TRBV27.png"),TRBV27, width = 7.5, height = 5)



# find top clones in each disease group ----
extract_clone_sizes <- function(tsv_gz,
                                clone_col = "clone_id",
                                productive_col = "productive") {
  df <- suppressMessages(readr::read_tsv(tsv_gz, col_types = cols(.default = "c")))
  if (productive_col %in% names(df)) {
    keep <- to_logical_true(df[[productive_col]])
    df <- df[keep %in% TRUE, , drop = FALSE]
  }
  if (nrow(df) == 0) return(NULL)
  # --- NEW: biologically meaningful clone identity ---
  df <- df %>%
    mutate(
      vjcdr3 = paste(v_call, j_call, cdr3_aa, sep = "_")
    )
  clone_freq <- df %>%
    count(vjcdr3, name = "n") %>%
    arrange(desc(n)) %>%
    rename(clone_id = vjcdr3)
  # extract participant label from file name
  bn <- basename(tsv_gz)
  m  <- str_match(bn, "(?i)bfi-(\\d+)")
  pid <- ifelse(!is.na(m[,2]),
                paste0("BFI-", str_pad(m[,2], width = 7, side = "left", pad = "0")),
                "BFI-NA")
  clone_freq %>% mutate(participant_label = pid)
}
# comment to avoid second running
#clone_sizes_with_name <- purrr::map_dfr(tcr_files, extract_clone_sizes) 
#clone_sizes_with_name2 <- clone_sizes_with_name %>% 
#  left_join(final_df %>% select(participant_label, group_binary, disease), 
#            by = "participant_label") %>% 
#  filter(!is.na(group_binary))
#readr::write_tsv(clone_sizes_with_name2, file.path(outdir, "clone_sizes_with_name.tsv"))
clone_sizes_with_name2<- read_tsv(file.path(outdir, "clone_sizes_with_name.tsv"))
top_clone_per_sample <- clone_sizes_with_name2 %>%
  group_by(participant_label) %>%
  slice_max(order_by = n, n = 1) %>% 
  ungroup()
top_clone_counts <- top_clone_per_sample %>%
  count(disease, clone_id, sort = TRUE)
top10_by_disease <- top_clone_counts %>%
  group_by(disease) %>%
  slice_max(n, n = 10)
public_top_clones <- top_clone_counts %>% 
  filter(n >= 2)
public_top_clones
#readr::write_tsv(public_top_clones, file.path(outdir, "public_top_clones.tsv"))
public_summary <- top_clone_counts %>%
  mutate(type = ifelse(n >= 2, "Public (≥2)", "Private (1)")) %>%
  count(disease, type)
p <- ggplot(public_summary, aes(x = disease, y = n, fill = type)) +
  geom_col(position = "stack") +
  theme_bw(base_size = 16) +
  scale_fill_manual(values = c("Private (1)" = "grey80", "Public (≥2)" = "#E64B35")) +
  labs(
    title = "Public vs Private Top Clonotypes",
    x = "Disease group",
    y = "Number of top clonotypes",
    fill = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
p
#ggsave(file.path(outdir, "Public_vs_Private_Top_Clonotypes.png"), p, width = 7.5, height = 5)


# cdr3 length ----
df_len <- final_df %>%
  filter(!is.na(cdr3_len_mean), !is.na(age.x))
lm_len <- lm(cdr3_len_mean ~ age.x, data = df_len)
df_len$cdr3_len_adj <- resid(lm_len)
p <- ggplot(df_len, aes(x = disease, y = cdr3_len_mean, fill = disease)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_bw(base_size = 16)
p
#ggsave(file.path(outdir, "cdr3_length_age_adjusted.png"), p, width = 12, height = 5)


# v-j pair ----
vj_df <- read_tsv(file.path(outdir, "vj_df.tsv"))
#vj_df <- clone_sizes_with_name2 %>%
#  separate(clone_id, into = c("V", "J", "CDR3"), sep = "_", remove = FALSE)
vj_usage <- vj_df %>%
  group_by(disease, V, J) %>%
  summarise(count = n(), .groups = "drop")
vj_mat <- vj_usage %>%
  pivot_wider(names_from = J, values_from = count, values_fill = 0)
ggplot(vj_usage, aes(x = J, y = V, fill = count)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = median(vj_usage$count),
    name = "Usage"
  ) +
  facet_wrap(~ disease, ncol = 3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    strip.background = element_rect(fill = "grey90")
  ) +
  labs(
    title = "V–J Pairing Usage Across Diseases",
    x = "J gene",
    y = "V gene"
  )
#readr::write_tsv(vj_df, file.path(outdir, "vj_df.tsv"))
# vj pair log-scale heatmap
vj_usage$usage_log <- log10(vj_usage$count + 1)
# top 20 V–J pairs per disease
vj_top20 <- vj_usage %>%
  group_by(disease) %>%
  slice_max(order_by = count, n = 20)
vj_top20_fig <- ggplot(vj_top20, aes(J, V, fill = count)) +
  geom_tile() +
  facet_wrap(~ disease) +
  scale_fill_viridis_c() +
  theme_minimal()
vj_top20_fig
#ggsave(file.path(outdir, "vj_top20.png"), vj_top20_fig, width = 15, height = 5)

# network 
vj_df <- vj_df %>%
  mutate(
    V = trimws(V),
    J = trimws(J),
    V = na_if(V, ""),
    J = na_if(J, "")
  ) %>%
  filter(!is.na(V), !is.na(J)) %>%
  select(disease, V, J, n)
vj_counts <- vj_df %>%
  group_by(disease, V, J) %>%
  summarise(usage = sum(n), .groups = "drop")
vj_top20 <- vj_counts %>%
  group_by(disease) %>%
  slice_max(order_by = usage, n = 20)
plot_vj_network <- function(df, disease_name) {
  df2 <- df %>%
    filter(disease == disease_name) %>%
    ungroup() %>%          
    mutate(
      V = trimws(V),
      J = trimws(J),
      V = na_if(V, ""),
      J = na_if(J, "")
    ) %>%
    filter(!is.na(V), !is.na(J))
  # edge list： V, J, usage, no disease 
  edge_df <- df2 %>% select(V, J, usage)
  # node list
  node_df <- tibble(
    name = unique(c(edge_df$V, edge_df$J)),
    type = case_when(
      grepl("^TRBV", name) ~ "V",
      grepl("^TRBJ", name) ~ "J",
      TRUE ~ "other"
    )
  ) %>%
    filter(!is.na(name))
  edge_df <- edge_df %>%
    filter(V %in% node_df$name, J %in% node_df$name)
  # build graph
  g <- graph_from_data_frame(
    d = edge_df,
    vertices = node_df,
    directed = FALSE
  )
  # plot
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = usage, color = usage), alpha = 0.7) +
    geom_node_point(aes(color = type), size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_width(range = c(0.5, 3)) +
    scale_colour_viridis_c(aesthetics = "edge_colour") +
    scale_color_manual(values = c("V" = "#1f78b4", "J" = "#33a02c", "other" = "gray")) +
    theme_void() +
    ggtitle(paste("V–J Network:", disease_name))
}
diseases <- unique(vj_top20$disease)
plots <- lapply(diseases, function(d) plot_vj_network(vj_top20, d))
plots[[1]]   
wrap_plots(plots, ncol = 3)
p <- wrap_plots(plots, ncol = 3)
#ggsave(file.path(outdir, "v_j_network.png"), p, width = 15, height = 8)
summarize_vj_modules <- function(df, top_n = 5) {
  df %>%
    group_by(disease, V, J) %>%
    summarise(total_usage = sum(usage), .groups = "drop") %>%
    arrange(disease, desc(total_usage)) %>%
    group_by(disease) %>%
    slice_max(order_by = total_usage, n = top_n) %>%
    mutate(
      module = paste0(V, " → ", J),
      rank = row_number()
    ) %>%
    select(disease, rank, module, total_usage)
}
vj_summary <- summarize_vj_modules(vj_top20, top_n = 5)
vj_summary
#readr::write_delim(vj_summary, file.path(outdir, "vj_summary.txt"))


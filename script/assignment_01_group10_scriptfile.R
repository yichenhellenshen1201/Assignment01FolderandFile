# ==============================================================================
# Title: BINF*6210 Assignment #2 - Collaborative Script Revision
# Primary Author: Yichen Shen (GitHub: yichenhellenshen1201)
#   - Role: Wrote the original Assignment #1 R script, including data import,
#     filtering, summarization, binding, and initial exploratory visualizations.
# Secondary Contributor: Purushoth Puvanendran (GitHub: Purushoth-Archive)
#   - Role: Revised comment formatting to organize outline display, integrated
#     T. nelsoni data, enhanced comparative visualizations with color coding
#     and labels, automated figure exports, and added Spearman rank correlation
#     testing.
# ==============================================================================

##***************************
## Assignment 01 - Group 10 
## 2026-11-02
##***************************

## Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())


rm(list = ls())


## Import TE annotation datasets -----

## _Import T. patagoniensis dataset
df_patagoniensis <- read_csv("../data/GCA_001447655.1.csv")
df_patagoniensis 

## _Import T. nelsoni dataset
df_nelsoni <- read_csv(file = "../data/trichinella_nelsoni_te.csv") 
df_nelsoni

## _Import T. spiralis dataset
df_spiralis = read_csv("../data/trichinella_spiralis_te.csv")
df_spiralis 


## Clean dataset columns -----
df_patagoniensis =
  df_patagoniensis |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
    no_bp_in_complement = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", ""))
  )

df_nelsoni =
  df_nelsoni |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
    no_bp_in_complement = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", ""))
  )

df_spiralis =
  df_spiralis |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
    no_bp_in_complement = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", ""))
  )


## Define non-TE classes to exclude -----
ls_non_te = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")
ls_non_te


## Sumarize data -----

## _Summarize T. spiralis TE annotations -----
df_spiralis_summary =
  df_spiralis |>
  filter(matching_class %notin% ls_non_te) |>
  group_by(matching_class) |>
  summarise(
    annotation_count = n(),
    total_length = sum(qry_width, na.rm = TRUE)
  ) |>
  mutate(species = "T. spiralis") |> 
  arrange(desc(annotation_count))


df_spiralis_summary


## _Summarize T. patagoniensis TE annotations -----
df_patagoniensis_summary =
  df_patagoniensis |>
  filter(matching_class %notin% ls_non_te) |>
  group_by(matching_class) |>
  summarise(
    annotation_count = n(),
    total_length = sum(qry_width, na.rm = TRUE)
  ) |>
  mutate(species = "T. patagoniensis") |> 
  arrange(desc(annotation_count))
df_patagoniensis_summary

## _Summarize T. nelsoni TE annotations -----
df_nelsoni_summary =
  df_nelsoni |>
  filter(matching_class %notin% ls_non_te) |>
  group_by(matching_class) |>
  summarise(
    annotation_count = n(),
    total_length = sum(qry_width, na.rm = TRUE)
  ) |>
  mutate(species = "T. nelsoni") |> 
  arrange(desc(annotation_count))
df_nelsoni_summary


## COMBINING -----

## _Combine species summaries -----
df_comparison =
  bind_rows(df_spiralis_summary,
            df_patagoniensis_summary, 
            df_nelsoni_summary)
df_comparison


## _Combine raw TE annotation datasets -----
df_TE_combined =
  bind_rows(
    df_spiralis |> mutate(species = "T. spiralis"),
    df_patagoniensis |> mutate(species = "T. patagoniensis"),
    df_nelsoni |> mutate(species = "T. nelsoni")
  )
df_TE_combined


## RANKING TE CLASSES -----

## _Rank TE classes by annotation count -----
df_comparison |>
  group_by(species) |>
  arrange(desc(annotation_count), .by_group = TRUE)


## _Rank TE classes by total annotated length -----
df_comparison |>
  group_by(species) |>
  arrange(desc(total_length), .by_group = TRUE)


# _Top 3 TE classes by total length for each species
df_comparison |> 
  group_by(species) |> 
  arrange(desc(total_length), .by_group = TRUE) |> 
  slice_head(n = 3)


## PLOTTING -----

# Create output directory if it doesn't already exist
dir.create("output", showWarnings = FALSE)

## _Plot distribution of TE annotation lengths -----
p_annotation_length = df_TE_combined |>
  filter(matching_class %notin% ls_non_te) |>
  ggplot(aes(x = qry_width, fill = species)) +
  geom_histogram(colour = "white") +
  stat_bin(
    geom = "text",
    aes(label = after_stat(ifelse(count > 0, count, ""))),
    vjust = -0.3,
    size = 2
  ) +
  scale_fill_brewer(palette = "Set2") +
  facet_grid(~ species) +
  labs(
    title = "Distribution of TE Annotation Lengths",
    x = "Annotated TE length (bp)",
    y = "Frequency",
    fill = "Species"
  )
# ADDITION: colors the 3 different species instead of leaving all graphs grey. Plus, I have added stat_bin to leave a number on the top of each histo bar. This should help with identifying our smaller values and contrast similar ones.

ggsave("output/te_length_distribution.png", plot = p_annotation_length, width = 12, height = 6, dpi = 300)

## _Plot annotation counts against total annotated lengths -----
p_TEabundance = df_comparison |> 
  ggplot(aes(x = annotation_count, y = total_length, color = species)) +
  geom_point() +
  geom_smooth() +
  facet_grid(~ species) +
  labs(title = "Comparison of TE Abundance & Length",
       x = "annotated count",
       y = "annotated sequence length(bp)")
# ADDITION: colors the 3 different species instead of leaving all graphs grey. 

ggsave("output/te_abundance_vs_length.png", plot = p_TEabundance, width = 12, height = 6, dpi = 300)

## _Plot average TE length by class -----
p_TE_avglength = df_comparison |>
  ggplot(aes(x = matching_class,
             y = total_length / annotation_count,
             fill = matching_class)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_d() +
  facet_grid(~ species) +
  labs(
    title = "Average TE Length by Class",
    x = "TE class",
    y = "Average annotated length (bp)",
    fill = "TE class"
  )
# ADDITION: Because now we have 3 graphs sharing a y-axis, to support contrasting the overwhelming amount of TE classess, I added color coding for each TE class to support visual comparison.

ggsave("output/avg_te_length_by_class.png", plot = p_TE_avglength, width = 12, height = 8, dpi = 300)


## STATISTICAL TESTING -----

## _Spearman correlation: Does TE abundance significantly correlate with total annotated length?
cor_abundance_length <- cor.test(df_comparison$annotation_count,
                                 df_comparison$total_length,
                                 method = "spearman"
                                 )
cor_abundance_length

# I chose Spearman's rank correlation to test how strongly TE class abundance rises or falls alongside total annotated sequence length. Plot 2 shows that our data is heavily skewed, because the data is not normally distributed and contains extreme outliers, a default Pearson test would be heavily distorted. Kendall is much stricter on small, pairwise rank differences and better suited for tiny sample sizes, whereas our goal is simply to test the overall upward trend across all TE classes using ranks.

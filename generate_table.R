############# SINGLE-CELL

library(dplyr)
library(gt)

import_summary <- tibble(
  "Method" = c("cellphonedb", "singlecellsignalR", "cellchat", "connectome", "natmi_specificity", "geometric_mean", "scseqcomm", "log2fc"),
  "How it computes interaction" = c(
    "Arithmetic mean of ligand/receptor expression. Multisubunit complexes are taken into account by taking the subunit with lowest expressed",
    "Regularized product score based on squared expression level of ligand and receptor",
    "Law of mass action taking into account ligand, receptors, subunits and cofactors",
    "Weight_norm is the product of normalized ligand and receptor expression.<br>Weight_scale (specificity) is the average of z-score of ligand expression in sending celltype and z-score of receptor expression in receiver celltype.",
    "Mean-expression edge weight as product of ligand and receptor expression <br>
    Specificity-based edge weight as the product of ligand and receptor specificities, defined as their mean expression in the interacting cell type divided by their total sum across all dataset cell types",
    "Geometric mean of ligand/receptor expression. Multisubunit complexes are taken into account by taking the subunit with lowest expressed",
    "Calculate ligand and receptor scores based on average cluster expression and compare it to a null distribution based on gene permutation within the cluster.<br> 
    The final score for a L-R pair is defined as the minimum of the ligand score in the sender cluster and receptor score in the receiver cluster",
    "Calculates the average one-versus-all log2 fold changes for sender and receiver celltypes."
  ),
  "Statistics" = c(
    "Label permutation",
    "Based on an estimated threshold (0.5) derived from authors benchmark",
    "Label permutation",
    "-",
    "Ranking of edge weights",
    "Label permutation from cellphonedb",
    "Manual score threshold > 0.8 as being significant according the original manuscript",
    "Manual thresholding > 1 was used"
  )
)

# Render block fixed for Markdown rendering compatibility
import_summary_table <- import_summary %>%
  gt() %>%  # FIXED: Removed groupname_col = "Method" so Markdown stays active
  fmt_markdown(columns = `How it computes interaction`) %>% # REQUIRED: Tells gt to process the <br> tag
  tab_header(
    title = md("**Table 1. Review of single cell computational methods for cell-cell communication**")
  ) %>%
  # Journal Style Sheet Rules (Nature Style Elements)
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(1.5),
    table_body.border.bottom.color = "black",
    table_body.border.bottom.width = px(1.5)
  ) %>%
  cols_align(
    align = "left", # Better readability for long character descriptions
    columns = `How it computes interaction`
  ) %>%
  cols_align(
    align = "center",
    columns = c(Method, Statistics)
  )

# Display the table in your Quarto preview
import_summary_table


import_summary_table %>%
  gtsave(
    filename = "Table1_SingleCell_Cell_Communication_Summary.pdf"
  )

library(gt)
library(dplyr)

# 1. Dataset setup
ccc_table_data <- data.frame(
  Method = c("CellChat", "CellPhoneDB", "NATMI", "SingleCellSignalR", "iTALK", "scMLnet", "ICELLNET", "NicheNet" , "Connectome", "CytoTalk"),
  `Xie (2023)` = c(TRUE, "TRUE_TOP3", "TRUE_TOP3", TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  `Wang (2022)`    = c(TRUE, TRUE, "TRUE_TOP3", TRUE, "TRUE_TOP3", TRUE, "TRUE_TOP3", TRUE, FALSE, FALSE),
  `Liu (2022)`   = c("TRUE_TOP3", "TRUE_TOP3", TRUE, TRUE, TRUE, TRUE, "TRUE_TOP3", "TRUE_TOP3", TRUE, TRUE),
  `ESICCC (2023)`  = c("TRUE_TOP3", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  check.names = FALSE # Preserves spaces in column names
)

# 2. Build the optimized gt table
gt_table <- ccc_table_data %>%
  gt(rowname_col = "Method") %>% 
  
  # Title
  tab_header(
    title = md("**Table3. Cell-cell communication tools evaluated across benchmarks and their top performers**")
  ) %>%
  
  # Transform into low-contrast HTML icons with perfectly equalized ticks
  text_transform(
    locations = cells_body(columns = c(`Xie (2023)`, `Wang (2022)`, `Liu (2022)`, `ESICCC (2023)`)), 
    fn = function(x) {
      lapply(x, function(cell_val) {
        html(
          case_when(
            # --- TOP 3 PERFORMANCE WINNER: Equalized large green tick inside a circle ---
            cell_val == "TRUE_TOP3" ~ "
              <span style='
                display: inline-block;
                border: 2px solid #499053;     /* Muted sage/olive ring */
                border-radius: 50%;            
                width: 44px;                   /* Larger container to prevent shrinking the text */
                height: 44px;
                line-height: 40px;             /* Matches container height to center the massive tick */
                text-align: center;            
                color: #499053; 
                font-weight: bold; 
                font-size: 34px;               /* IDENTICAL size to standard tick */
                box-sizing: border-box;
              '>&#10004;</span>",
            
            # --- REGULAR METHOD INCLUSION: Extra Large, low-contrast Green Tick ---
            cell_val == "TRUE" ~ "
              <span style='
                display: inline-block;
                width: 44px;                   /* Identical box footprint as circled version */
                height: 44px;
                line-height: 44px;
                color: #499053; 
                font-weight: bold; 
                font-size: 34px;               /* IDENTICAL size to top-3 tick */
                box-sizing: border-box;
              '>&#10004;</span>",
            
            # --- METHOD NOT EVALUATED: Large, low-contrast Slate/Grey Cross ---
            cell_val == "FALSE" ~ "
              <span style='
                display: inline-block;
                width: 44px; 
                height: 44px;
                line-height: 44px;
                color: #8c969e; 
                font-weight: bold; 
                font-size: 28px;
                box-sizing: border-box;
              '>&#10006;</span>",
            
            # Absolute fallback
            TRUE ~ cell_val
          )
        )
      })
    }
  ) %>%
  
  # Center align the column contents
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  # Professional academic styling
  opt_table_outline() %>%
  tab_options(
    column_labels.background.color = "#F8F9FA",
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    table.font.size = "14px"
  )

# Export directly to PDF
gt_table %>%
  gtsave(
    filename = "Table3_AcrossBenchmarksComparisson.pdf"
  )
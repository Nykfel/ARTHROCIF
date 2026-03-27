library(dplyr)
library(tidyr)

long <- read.delim(
  "C:/Users/2026mr001/Desktop/Stage/Identification cif/Rstudio/long.txt",
  check.names = FALSE
)

df <- long

# Identification automatique des colonnes dynamiques (au-delà de "taille")
idx_taille <- which(names(df) == "taille")
cols_after_taille <- names(df)[(idx_taille + 1):ncol(df)]

# Colonnes START/END uniquement
cols_start_end <- cols_after_taille[grepl("_(START|END)$", cols_after_taille)]

# Colonnes fixes
cols_fixed <- names(df)[1:idx_taille]

# Fonction de nettoyage
clean_domaine <- function(x) {
  x %>%
    gsub("\\.", "_", .) %>%
    gsub("_\\d+$", "", .)
}

# Pivot principal
df_result <- df %>%
  pivot_longer(
    cols = all_of(cols_start_end),
    names_to = c("Domain", ".value"),
    names_pattern = "(.+?)_(START|END)"
  ) %>%
  filter(!is.na(START) | !is.na(END)) %>%
  mutate(Domain = clean_domaine(Domain)) %>%
  select(all_of(cols_fixed), Domain, START, END) %>%
  arrange(across(all_of(cols_fixed)), Domain)

df <- df_result

unique(df$Domain)

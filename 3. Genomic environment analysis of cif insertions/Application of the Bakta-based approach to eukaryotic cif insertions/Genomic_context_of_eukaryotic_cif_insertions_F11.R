library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)

d <- read.delim("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/cif_like/data_CDS.txt")
e <- read.delim("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/cif_like/e.txt")

d <- d %>%
filter(Scaffold_length > 1000000)

par(mfrow = c(1, 1))

# Ajout de l'accession depuis e
d <- d %>%
  left_join(
    e %>% select(Abrv, Organism, Accession, Cif_Type),
    by = c("Souche" = "Abrv")
  )

# d <- d %>%
#   filter(Cif_Type == "I")

# Remplacement des abréviations par les noms complets
d <- d %>%
  mutate(
    Souche = Organism
  ) %>%
  select(-Organism)

# Garder uniquement les scaffolds > 100 000 pb
d <- d %>%
  filter(Scaffold_length > 99999)

# Séparation CDS annotés / protéines hypothétiques
d <- d %>%
  mutate(
    hypothetical = ifelse(
      Product %in% c(
        "hypothetical protein",
        "Conserved hypothetical exported protein",
        "(pseudo) hypothetical protein"
      ),
      1, 0
    )
  )

# Retirer les éléments qui ne sont pas des CDS
d <- d %>%
  filter(!Type %in% c(
    "ncRNA-region", "tRNA", "sorf", "ncRNA", "tmRNA",
    "rRNA", "crispr", "crispr-repeat", "crispr-spacer",
    "assembly_gap", "oriC", "gap"
  ))

all(d$Type == "cds")
unique(d$Type)

# Paramètres
bin_size <- 10000

# Tableau des scaffolds : séparation par Souche + Accession
contigs <- d %>%
  group_by(Souche, Accession) %>%
  summarise(
    Scaffold_length = max(Scaffold_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Souche, Accession, Scaffold_length)

# Longueur maximale globale
max_len <- ifelse(
  max(contigs$Scaffold_length, na.rm = TRUE) <= 1000000,
  1000000,
  max(contigs$Scaffold_length, na.rm = TRUE)
)

# Positions moyennes des couples cif par Souche + Accession
# Pour chaque numéro, on prend toutes les positions CifA/CifB disponibles,
# puis on calcule le milieu entre la valeur minimale et la valeur maximale.
cif_points <- d %>%
  group_by(Cif_Type, Souche, Accession) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    Cif_Type,
    Souche,
    Accession,
    matches("Cif[AB]_(Start|End)[0-9]+_relative")
  ) %>%
  pivot_longer(
    cols = -c(Cif_Type, Souche, Accession),
    names_to = c("cif", "coord", "id"),
    names_pattern = "(Cif[AB])_(Start|End)([0-9]+)_relative",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  group_by(Cif_Type, Souche, Accession, id) %>%
  summarise(
    min_pos = min(value, na.rm = TRUE),
    max_pos = max(value, na.rm = TRUE),
    pos = (min_pos + max_pos) / 2,
    has_CifA = any(cif == "CifA"),
    has_CifB = any(cif == "CifB"),
    .groups = "drop"
  ) %>%
  mutate(
    type_triangle = case_when(
      has_CifA & has_CifB ~ "CifA + CifB",
      !has_CifA & has_CifB ~ "CifB seul",
      has_CifA & !has_CifB ~ "CifA seul",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(type_triangle))

# Création des fenêtres par Souche + Accession
bins_by_souche <- contigs %>%
  rowwise() %>%
  mutate(
    bin_start = list(seq(1, Scaffold_length, by = bin_size))
  ) %>%
  unnest(bin_start) %>%
  ungroup() %>%
  mutate(
    bin_end = pmin(bin_start + bin_size - 1, Scaffold_length),
    bin_mid = (bin_start + bin_end) / 2,
    bin_length = bin_end - bin_start + 1
  )

grep("Cif", names(d), value = TRUE)

# Attribution des CDS à une fenêtre selon leur position de départ
d2 <- d %>%
  mutate(
    bin_start = floor((Start - 1) / bin_size) * bin_size + 1,
    CDS_length = Stop - Start + 1
  )

# Nombre de CDS annotés par fenêtre
gene_density_annotated <- d2 %>%
  filter(hypothetical == 0) %>%
  group_by(Souche, Accession, bin_start) %>%
  summarise(
    n_CDS = n(),
    .groups = "drop"
  ) %>%
  right_join(
    bins_by_souche,
    by = c("Souche", "Accession", "bin_start")
  ) %>%
  mutate(
    n_CDS = replace_na(n_CDS, 0)
  ) %>%
  arrange(Souche, Accession, bin_start)

# Découpage virtuel des CDS selon les fenêtres
d_coverage <- d %>%
  select(Souche, Accession, Start, Stop) %>%
  inner_join(
    bins_by_souche,
    by = c("Souche", "Accession"),
    relationship = "many-to-many"
  ) %>%
  filter(
    Start <= bin_end,
    Stop >= bin_start
  ) %>%
  mutate(
    Start_cut = pmax(Start, bin_start),
    Stop_cut = pmin(Stop, bin_end),
    CDS_length_cut = Stop_cut - Start_cut + 1
  )

# Coding coverage par fenêtre
coding_coverage <- d_coverage %>%
  group_by(Souche, Accession, bin_start, bin_end, bin_mid, bin_length) %>%
  summarise(
    total_cds_length = sum(CDS_length_cut, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coding_coverage = total_cds_length / bin_length
  ) %>%
  right_join(
    bins_by_souche,
    by = c("Souche", "Accession", "bin_start", "bin_end", "bin_mid", "bin_length")
  ) %>%
  mutate(
    total_cds_length = replace_na(total_cds_length, 0),
    coding_coverage = replace_na(coding_coverage, 0)
  ) %>%
  group_by(Souche, Accession, bin_start, bin_end, bin_mid, bin_length) %>%
  summarise(
    total_cds_length = sum(total_cds_length, na.rm = TRUE),
    coding_coverage = sum(coding_coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Souche, Accession, bin_start)

# Vérification optionnelle : doit retourner 0 ligne
coding_coverage %>%
  count(Souche, Accession, bin_mid) %>%
  filter(n > 1)

# Valeurs max globales
# ymax_a <- max(gene_density_annotated$n_CDS, na.rm = TRUE) * 1.05
ymax_a <- 20
ymax_b <- 1

# Fonction qui crée les 3 graphiques pour une combinaison Souche + Accession
make_souche_panel <- function(souche_name, accession_name) {
  
  panel_title <- paste0(souche_name, " | ", accession_name)
  
  d_souche <- d %>%
    filter(Souche == souche_name, Accession == accession_name)
  
  annotated_souche <- gene_density_annotated %>%
    filter(Souche == souche_name, Accession == accession_name)
  
  coverage_souche <- coding_coverage %>%
    filter(Souche == souche_name, Accession == accession_name)
  
  cif_souche <- cif_points %>%
    filter(Souche == souche_name, Accession == accession_name)
  
  p_annotated <- ggplot(annotated_souche, aes(x = bin_mid, y = n_CDS)) +
    geom_area(fill = "red", alpha = 0.5) +
    geom_line(linewidth = 0.8, color = "red") +
    scale_x_continuous(
      limits = c(1, max_len),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, ymax_a),
      expand = c(0, 0)
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = "a"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  p_coverage <- ggplot(coverage_souche, aes(x = bin_mid, y = coding_coverage)) +
    geom_area(fill = "grey40", alpha = 0.5) +
    geom_line(linewidth = 0.8, color = "grey20") +
    scale_x_continuous(
      limits = c(1, max_len),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, ymax_b),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "b"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  p_scaffold <- ggplot() +
    geom_segment(
      data = contigs %>% filter(Souche == souche_name, Accession == accession_name),
      aes(
        x = 1,
        xend = Scaffold_length,
        y = 1,
        yend = 1
      ),
      linewidth = 1,
      color = "black"
    ) +
    geom_segment(
      data = d_souche,
      aes(
        x = Start,
        xend = Stop,
        y = 1,
        yend = 1,
        color = factor(hypothetical)
      ),
      linewidth = 7
    ) +
    geom_point(
      data = cif_souche,
      aes(
        x = pos,
        y = 1.25,
        fill = type_triangle
      ),
      shape = 25,
      size = 4,
      color = "black",
      alpha = 1,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = c("0" = "red", "1" = "grey25")
    ) +
    scale_fill_manual(
      values = c(
        "CifA + CifB" = "#FF9B85",
        "CifB seul" = "#FFF385",
        "CifA seul" = "green"
      )
    ) +
    scale_x_continuous(
      limits = c(1, max_len),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.7, 1.3),
      breaks = 1,
      labels = panel_title
    ) +
    labs(
      x = NULL,
      y = "c"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  p_annotated / p_coverage / p_scaffold +
    plot_layout(heights = c(1, 1, 0.5))
}

# Création d’un bloc de 3 graphiques par combinaison Souche + Accession
souches_accessions <- contigs %>%
  select(Souche, Accession) %>%
  distinct()

plots_by_souche <- lapply(
  seq_len(nrow(souches_accessions)),
  function(i) {
    make_souche_panel(
      souches_accessions$Souche[i],
      souches_accessions$Accession[i]
    )
  }
)

# Assemblage final
figure_finale <- wrap_plots(
  plots_by_souche,
  ncol = 1
)

figure_finale

# Valeur maximale réelle du graphique A
cat(
  "Valeur maximale de CDS annotés dans le graphique A :",
  max(gene_density_annotated$n_CDS, na.rm = TRUE),
  "\n"
)

#14 par x  * 2.5 pouces


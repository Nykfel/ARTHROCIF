library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)

d <- read.delim("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/Genomes_references/data_CDS.txt")
info <- read.delim("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/Genomes_references/info.txt")
info$Sequence <- sub("/.*", "", info$Sequence)
setdiff(unique(d$Souche), unique(info$Sequence))

#d <- d %>%
#  filter(Clade != "Bacteria")

#d <- d %>%
#  filter(Clade != "Eucaryota")

d$Souche <- info$Name[match(d$Souche, info$Sequence)]
par(mfrow = c(1, 1))


# Si ton tableau d est déjà chargé, tu peux laisser cette ligne commentée
# d <- read.delim("CHEMIN/VERS/TON/NOUVEAU_FICHIER.txt")

par(mfrow = c(1, 1))

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

# Garder uniquement les vrais CDS
d <- d %>%
  filter(tolower(Type) == "cds")

# Vérifications
all(tolower(d$Type) == "cds")
unique(d$Type)

# Paramètres
bin_size <- 10000

# Tableau des scaffolds : un scaffold par Souche
contigs <- d %>%
  group_by(Clade, Souche) %>%
  summarise(
    Scaffold_length = max(Scaffold_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Scaffold_length))

# Longueur maximale globale
max_len <- ifelse(
  max(contigs$Scaffold_length, na.rm = TRUE) <= 1000000,
  1000000,
  max(contigs$Scaffold_length, na.rm = TRUE)
)

# Création des fenêtres par Souche
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

# Attribution des CDS à une fenêtre selon leur position de départ
d2 <- d %>%
  mutate(
    bin_start = floor((Start - 1) / bin_size) * bin_size + 1,
    CDS_length = Stop - Start + 1
  )

# Nombre de CDS annotés par fenêtre
gene_density_annotated <- d2 %>%
  filter(hypothetical == 0) %>%
  group_by(Clade, Souche, bin_start) %>%
  summarise(
    n_CDS = n(),
    .groups = "drop"
  ) %>%
  right_join(
    bins_by_souche,
    by = c("Clade", "Souche", "bin_start")
  ) %>%
  mutate(
    n_CDS = replace_na(n_CDS, 0)
  ) %>%
  arrange(Clade, Souche, bin_start)

# Découpage virtuel des CDS selon les fenêtres
d_coverage <- d %>%
  select(Clade, Souche, Start, Stop) %>%
  inner_join(
    bins_by_souche,
    by = c("Clade", "Souche"),
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
  group_by(Clade, Souche, bin_start, bin_end, bin_mid, bin_length) %>%
  summarise(
    total_cds_length = sum(CDS_length_cut, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coding_coverage = total_cds_length / bin_length
  ) %>%
  right_join(
    bins_by_souche,
    by = c("Clade", "Souche", "bin_start", "bin_end", "bin_mid", "bin_length")
  ) %>%
  mutate(
    total_cds_length = replace_na(total_cds_length, 0),
    coding_coverage = replace_na(coding_coverage, 0)
  ) %>%
  arrange(Clade, Souche, bin_start)

# Vérification optionnelle : doit retourner 0 ligne
coding_coverage %>%
  count(Clade, Souche, bin_mid) %>%
  filter(n > 1)

# Valeurs max globales
# ymax_a <- max(gene_density_annotated$n_CDS, na.rm = TRUE) * 1.05
ymax_a <- 30
ymax_b <- 1

# Fonction qui crée les 3 graphiques pour une Souche
make_souche_panel <- function(clade_name, souche_name) {
  
  panel_title <- paste0(souche_name, " | ", clade_name)
  
  d_souche <- d %>%
    filter(Clade == clade_name, Souche == souche_name)
  
  annotated_souche <- gene_density_annotated %>%
    filter(Clade == clade_name, Souche == souche_name)
  
  coverage_souche <- coding_coverage %>%
    filter(Clade == clade_name, Souche == souche_name)
  
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
      data = contigs %>% filter(Clade == clade_name, Souche == souche_name),
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
    scale_color_manual(
      values = c("0" = "red", "1" = "grey25")
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

# Liste des souches
souches <- contigs %>%
  select(Clade, Souche) %>%
  distinct()

# Création des graphiques
plots_by_souche <- lapply(
  seq_len(nrow(souches)),
  function(i) {
    make_souche_panel(
      souches$Clade[i],
      souches$Souche[i]
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
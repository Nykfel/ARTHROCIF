library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)

# ============================================================
# 1. IMPORTATION DES DONNÉES
# ============================================================

d <- read.delim(
  "C:/Users/2026mr001/Desktop/guipguiphp/data_CDS.txt",
  stringsAsFactors = FALSE
)

info <- read.delim(
  "C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/Genomes_references/info.txt",
  stringsAsFactors = FALSE
)

# Suppression de la partie située après "/" dans Sequence
info$Sequence <- sub("/.*", "", info$Sequence)

# Vérification des identifiants qui ne correspondent pas
setdiff(unique(d$Souche), unique(info$Sequence))

# Remplacement des identifiants par les noms des souches
d$Souche <- info$Name[match(d$Souche, info$Sequence)]

# Vérification des correspondances manquantes
d %>%
  filter(is.na(Souche))

# ============================================================
# 2. PRÉPARATION DES DONNÉES
# ============================================================

# Conservation uniquement des CDS
d <- d %>%
  filter(tolower(Type) == "cds")

# Vérifications
all(tolower(d$Type) == "cds")
unique(d$Type)

# Mise en ordre des coordonnées :
# utile si Start > Stop pour certains CDS situés sur le brin inverse
d <- d %>%
  mutate(
    CDS_start = pmin(Start, Stop),
    CDS_end = pmax(Start, Stop)
  )

# Classification des protéines hypothétiques
d <- d %>%
  mutate(
    hypothetical = ifelse(
      Product %in% c(
        "hypothetical protein",
        "Conserved hypothetical exported protein",
        "(pseudo) hypothetical protein"
      ),
      1,
      0
    )
  )

# Vérification du nombre de CDS dans chaque catégorie
d %>%
  count(hypothetical)

# ============================================================
# 3. PARAMÈTRES
# ============================================================

# Taille des fenêtres : 10 kb
bin_size <- 10000

# Limites verticales des lignes A et B
ymax_a <- 1
ymax_b <- 1

# ============================================================
# 4. TABLEAU DES SCAFFOLDS
# ============================================================

# Le code suppose ici qu'une souche correspond à un seul scaffold
contigs <- d %>%
  group_by(Clade, Souche) %>%
  summarise(
    Scaffold_length = max(Scaffold_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Scaffold_length))

# Longueur maximale commune utilisée sur les graphiques
max_len <- ifelse(
  max(contigs$Scaffold_length, na.rm = TRUE) <= 1000000,
  1000000,
  max(contigs$Scaffold_length, na.rm = TRUE)
)

# ============================================================
# 5. CRÉATION DES FENÊTRES DE 10 KB
# ============================================================

bins_by_souche <- contigs %>%
  rowwise() %>%
  mutate(
    bin_start = list(
      seq(
        from = 1,
        to = Scaffold_length,
        by = bin_size
      )
    )
  ) %>%
  unnest(bin_start) %>%
  ungroup() %>%
  mutate(
    bin_end = pmin(
      bin_start + bin_size - 1,
      Scaffold_length
    ),
    bin_mid = (bin_start + bin_end) / 2,
    bin_length = bin_end - bin_start + 1
  )

# Vérification de la longueur des fenêtres
bins_by_souche %>%
  count(bin_length)

# ============================================================
# 6. FONCTION CALCULANT LA LONGUEUR DE L'UNION DES INTERVALLES
# ============================================================

# Cette fonction fusionne les intervalles qui se chevauchent.
# Une position appartenant à plusieurs CDS n'est donc comptée
# qu'une seule fois.

union_length <- function(start, end) {
  
  # Suppression des valeurs manquantes
  keep <- !is.na(start) & !is.na(end)
  
  start <- start[keep]
  end <- end[keep]
  
  # Aucun intervalle
  if (length(start) == 0) {
    return(0)
  }
  
  # Mise dans le bon ordre des coordonnées
  interval_start <- pmin(start, end)
  interval_end <- pmax(start, end)
  
  # Tri selon le début puis la fin des intervalles
  ord <- order(interval_start, interval_end)
  
  interval_start <- interval_start[ord]
  interval_end <- interval_end[ord]
  
  # Initialisation avec le premier intervalle
  current_start <- interval_start[1]
  current_end <- interval_end[1]
  
  total_length <- 0
  
  # Fusion des intervalles chevauchants
  if (length(interval_start) > 1) {
    
    for (i in 2:length(interval_start)) {
      
      # Intervalles chevauchants ou directement adjacents
      if (interval_start[i] <= current_end + 1) {
        
        current_end <- max(
          current_end,
          interval_end[i]
        )
        
      } else {
        
        # Ajout de la longueur de l'intervalle précédent
        total_length <- total_length +
          (current_end - current_start + 1)
        
        # Début d'un nouvel intervalle
        current_start <- interval_start[i]
        current_end <- interval_end[i]
      }
    }
  }
  
  # Ajout du dernier intervalle
  total_length +
    (current_end - current_start + 1)
}

# ============================================================
# 7. FONCTION DE CALCUL DE LA COUVERTURE UNIQUE PAR FENÊTRE
# ============================================================

calculate_unique_coverage <- function(
    data,
    bins,
    annotated_only = FALSE
) {
  
  data_selected <- data
  
  # Pour la ligne A :
  # conservation uniquement des CDS non hypothétiques
  if (annotated_only) {
    
    data_selected <- data_selected %>%
      filter(hypothetical == 0)
  }
  
  # Association de chaque CDS à toutes les fenêtres qu'il traverse
  intersections <- data_selected %>%
    select(
      Clade,
      Souche,
      CDS_start,
      CDS_end
    ) %>%
    inner_join(
      bins,
      by = c("Clade", "Souche"),
      relationship = "many-to-many"
    ) %>%
    filter(
      CDS_start <= bin_end,
      CDS_end >= bin_start
    ) %>%
    mutate(
      # Partie du CDS réellement comprise dans la fenêtre
      Start_cut = pmax(CDS_start, bin_start),
      Stop_cut = pmin(CDS_end, bin_end)
    )
  
  # Fusion des régions chevauchantes dans chaque fenêtre
  coverage <- intersections %>%
    group_by(
      Clade,
      Souche,
      bin_start,
      bin_end,
      bin_mid,
      bin_length
    ) %>%
    summarise(
      covered_bp = union_length(
        Start_cut,
        Stop_cut
      ),
      .groups = "drop"
    ) %>%
    right_join(
      bins,
      by = c(
        "Clade",
        "Souche",
        "bin_start",
        "bin_end",
        "bin_mid",
        "bin_length"
      )
    ) %>%
    mutate(
      covered_bp = replace_na(covered_bp, 0),
      
      # Proportion de positions couvertes par au moins un CDS
      coverage = covered_bp / bin_length
    ) %>%
    arrange(
      Clade,
      Souche,
      bin_start
    )
  
  return(coverage)
}

# ============================================================
# 8. LIGNE A : COUVERTURE PAR LES CDS NON HYPOTHÉTIQUES
# ============================================================

annotated_coverage <- calculate_unique_coverage(
  data = d,
  bins = bins_by_souche,
  annotated_only = TRUE
)

# Signification :
# proportion des positions de chaque fenêtre couvertes
# par au moins un CDS non hypothétique

# ============================================================
# 9. LIGNE B : COUVERTURE PAR TOUS LES CDS
# ============================================================

total_coding_coverage <- calculate_unique_coverage(
  data = d,
  bins = bins_by_souche,
  annotated_only = FALSE
)

# Signification :
# proportion des positions de chaque fenêtre couvertes
# par au moins un CDS, protéines hypothétiques incluses

# ============================================================
# 10. VÉRIFICATIONS DES COUVERTURES
# ============================================================

# Aucune couverture ne doit être inférieure à 0
# ou supérieure à 1

annotated_coverage %>%
  filter(
    coverage < 0 |
      coverage > 1
  )

total_coding_coverage %>%
  filter(
    coverage < 0 |
      coverage > 1
  )

# La couverture des CDS non hypothétiques ne doit jamais
# dépasser la couverture totale

coverage_check <- annotated_coverage %>%
  select(
    Clade,
    Souche,
    bin_start,
    annotated_coverage = coverage
  ) %>%
  left_join(
    total_coding_coverage %>%
      select(
        Clade,
        Souche,
        bin_start,
        total_coverage = coverage
      ),
    by = c(
      "Clade",
      "Souche",
      "bin_start"
    )
  ) %>%
  filter(
    annotated_coverage >
      total_coverage + 1e-12
  )

coverage_check

# Cette commande doit normalement renvoyer 0 ligne

# Vérification des éventuels doublons de fenêtres

annotated_coverage %>%
  count(
    Clade,
    Souche,
    bin_start
  ) %>%
  filter(n > 1)

total_coding_coverage %>%
  count(
    Clade,
    Souche,
    bin_start
  ) %>%
  filter(n > 1)

# ============================================================
# 11. FONCTION DE CRÉATION DU PANNEAU D'UNE SOUCHE
# ============================================================

make_souche_panel <- function(
    clade_name,
    souche_name
) {
  
  panel_title <- paste0(
    souche_name,
    " | ",
    clade_name
  )
  
  # CDS de la souche pour la ligne C
  d_souche <- d %>%
    filter(
      Clade == clade_name,
      Souche == souche_name
    )
  
  # Couverture des CDS non hypothétiques pour la ligne A
  annotated_souche <- annotated_coverage %>%
    filter(
      Clade == clade_name,
      Souche == souche_name
    )
  
  # Couverture totale des CDS pour la ligne B
  coverage_souche <- total_coding_coverage %>%
    filter(
      Clade == clade_name,
      Souche == souche_name
    )
  
  # Informations sur le scaffold
  contig_souche <- contigs %>%
    filter(
      Clade == clade_name,
      Souche == souche_name
    )
  
  # ----------------------------------------------------------
  # LIGNE A
  # Proportion de positions couvertes par au moins un CDS
  # non hypothétique dans chaque fenêtre de 10 kb
  # ----------------------------------------------------------
  
  p_annotated <- ggplot(
    annotated_souche,
    aes(
      x = bin_mid,
      y = coverage
    )
  ) +
    geom_area(
      fill = "red",
      alpha = 0.5
    ) +
    geom_line(
      linewidth = 0.8,
      color = "red"
    ) +
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
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(
        angle = 0,
        vjust = 0.5
      )
    )
  
  # ----------------------------------------------------------
  # LIGNE B
  # Proportion de positions couvertes par au moins un CDS
  # dans chaque fenêtre de 10 kb
  # ----------------------------------------------------------
  
  p_coverage <- ggplot(
    coverage_souche,
    aes(
      x = bin_mid,
      y = coverage
    )
  ) +
    geom_area(
      fill = "grey40",
      alpha = 0.5
    ) +
    geom_line(
      linewidth = 0.8,
      color = "grey20"
    ) +
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
      axis.title.y = element_text(
        angle = 0,
        vjust = 0.5
      )
    )
  
  # ----------------------------------------------------------
  # LIGNE C
  # Positions exactes des CDS sur le scaffold
  # Rouge : CDS non hypothétiques
  # Gris : CDS hypothétiques
  # ----------------------------------------------------------
  
  p_scaffold <- ggplot() +
    
    # Scaffold complet
    geom_segment(
      data = contig_souche,
      aes(
        x = 1,
        xend = Scaffold_length,
        y = 1,
        yend = 1
      ),
      linewidth = 1,
      color = "black"
    ) +
    
    # CDS
    geom_segment(
      data = d_souche,
      aes(
        x = CDS_start,
        xend = CDS_end,
        y = 1,
        yend = 1,
        color = factor(hypothetical)
      ),
      linewidth = 7
    ) +
    
    scale_color_manual(
      values = c(
        "0" = "red",
        "1" = "grey25"
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
      axis.title.y = element_text(
        angle = 0,
        vjust = 0.5
      )
    )
  
  # Assemblage des trois lignes
  p_annotated /
    p_coverage /
    p_scaffold +
    plot_layout(
      heights = c(1, 1, 0.5)
    )
}

# ============================================================
# 12. LISTE DES SOUCHES
# ============================================================

souches <- contigs %>%
  select(
    Clade,
    Souche
  ) %>%
  distinct()

# ============================================================
# 13. CRÉATION DES PANNEAUX
# ============================================================

plots_by_souche <- lapply(
  seq_len(nrow(souches)),
  function(i) {
    
    make_souche_panel(
      clade_name = souches$Clade[i],
      souche_name = souches$Souche[i]
    )
  }
)

# ============================================================
# 14. ASSEMBLAGE FINAL
# ============================================================

figure_finale <- wrap_plots(
  plots_by_souche,
  ncol = 1
)

figure_finale

# ============================================================
# 15. VALEURS MAXIMALES OBSERVÉES
# ============================================================

cat(
  "Couverture maximale des CDS non hypothétiques :",
  max(
    annotated_coverage$coverage,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Couverture maximale de l'ensemble des CDS :",
  max(
    total_coding_coverage$coverage,
    na.rm = TRUE
  ),
  "\n"
)

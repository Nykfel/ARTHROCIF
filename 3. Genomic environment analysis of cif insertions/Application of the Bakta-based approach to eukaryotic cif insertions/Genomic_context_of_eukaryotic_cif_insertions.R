library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)

# ============================================================
# 1. IMPORTATION DES DONNÉES
# ============================================================

d <- read.delim(
  "C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/cif_like/data.txt",
  stringsAsFactors = FALSE
)

e <- read.delim(
  "C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/cif_like/e.txt",
  stringsAsFactors = FALSE
)

par(mfrow = c(1, 1))


# ============================================================
# 2. PREMIER FILTRE SUR LA LONGUEUR DES SCAFFOLDS
# ============================================================

# Ce filtre est conservé tel qu'il figurait dans ton script.
# Il implique que seuls les scaffolds strictement supérieurs
# à 1 Mb sont conservés.
d <- d %>%
  filter(Scaffold_length > 1000000)


# ============================================================
# 3. AJOUT DES INFORMATIONS DU TABLEAU e
# ============================================================

# Préparation du tableau utilisé pour la jointure
e_join <- e %>%
  select(
    Abrv,
    Organism,
    Accession,
    Cif_Type
  ) %>%
  distinct()

# Vérification : une même abréviation est-elle associée
# à plusieurs lignes différentes dans e ?
e_duplicates <- e_join %>%
  count(Abrv, name = "n") %>%
  filter(n > 1)

if (nrow(e_duplicates) > 0) {
  warning(
    "Certaines valeurs de Abrv correspondent à plusieurs lignes dans e. ",
    "La jointure peut donc dupliquer certaines lignes de d."
  )

  print(e_duplicates)
}

# Ajout du nom de l'organisme, de l'accession et du type de cif
d <- d %>%
  left_join(
    e_join,
    by = c("Souche" = "Abrv")
  )

# Filtre facultatif sur un type particulier de cif
# d <- d %>%
#   filter(Cif_Type == "I")


# ============================================================
# 4. REMPLACEMENT DES ABRÉVIATIONS PAR LES NOMS COMPLETS
# ============================================================

d <- d %>%
  mutate(
    Souche = Organism
  ) %>%
  select(-Organism)

# Vérification des correspondances manquantes
missing_metadata <- d %>%
  filter(
    is.na(Souche) |
      is.na(Accession)
  ) %>%
  distinct(Souche, Accession)

if (nrow(missing_metadata) > 0) {
  warning(
    "Certaines lignes n'ont pas reçu de nom d'organisme ",
    "ou d'accession après la jointure."
  )

  print(missing_metadata)
}


# ============================================================
# 5. SECOND FILTRE SUR LA LONGUEUR DES SCAFFOLDS
# ============================================================

# Ce filtre est également conservé.
# Il est toutefois redondant avec le filtre > 1 000 000
# appliqué précédemment.
d <- d %>%
  filter(Scaffold_length > 99999)


# ============================================================
# 6. CONSERVATION DES CDS
# ============================================================

# Retrait explicite des types non CDS présents dans le script initial
d <- d %>%
  filter(
    !Type %in% c(
      "ncRNA-region",
      "tRNA",
      "sorf",
      "ncRNA",
      "tmRNA",
      "rRNA",
      "crispr",
      "crispr-repeat",
      "crispr-spacer",
      "assembly_gap",
      "oriC",
      "gap"
    )
  )

# Conservation stricte des lignes de type CDS
# Cette deuxième étape garantit qu'un éventuel autre type
# non prévu dans la liste précédente ne sera pas conservé.
d <- d %>%
  filter(tolower(Type) == "cds")

# Vérifications
all(tolower(d$Type) == "cds")
unique(d$Type)

if (nrow(d) == 0) {
  stop(
    "Aucun CDS ne reste après l'application des filtres."
  )
}


# ============================================================
# 7. PRÉPARATION DES COORDONNÉES DES CDS
# ============================================================

# CDS_start et CDS_end restent dans le bon ordre,
# y compris pour les CDS situés sur le brin inverse
d <- d %>%
  mutate(
    CDS_start = pmin(Start, Stop),
    CDS_end = pmax(Start, Stop)
  )

# Vérification des coordonnées manquantes
invalid_coordinates <- d %>%
  filter(
    is.na(CDS_start) |
      is.na(CDS_end) |
      is.na(Scaffold_length)
  )

if (nrow(invalid_coordinates) > 0) {
  warning(
    nrow(invalid_coordinates),
    " CDS possèdent des coordonnées ou une longueur de scaffold manquantes. ",
    "Ces lignes sont retirées."
  )

  d <- d %>%
    filter(
      !is.na(CDS_start),
      !is.na(CDS_end),
      !is.na(Scaffold_length)
    )
}


# ============================================================
# 8. SÉPARATION DES CDS NON HYPOTHÉTIQUES ET HYPOTHÉTIQUES
# ============================================================

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
print(
  d %>%
    count(hypothetical)
)


# ============================================================
# 9. PARAMÈTRES GÉNÉRAUX
# ============================================================

# Taille des fenêtres : 10 kb
bin_size <- 10000

# Les lignes A et B représentent désormais des proportions
# nécessairement comprises entre 0 et 1.
ymax_a <- 1
ymax_b <- 1


# ============================================================
# 10. TABLEAU DES SCAFFOLDS
# ============================================================

# Un scaffold est défini par la combinaison :
# Souche + Accession
contigs <- d %>%
  group_by(
    Souche,
    Accession
  ) %>%
  summarise(
    Scaffold_length = max(
      Scaffold_length,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Souche,
    Accession,
    Scaffold_length
  )

if (nrow(contigs) == 0) {
  stop(
    "Aucun scaffold ne reste après la préparation des données."
  )
}


# ============================================================
# 11. DÉFINITION DE L'ÉCHELLE HORIZONTALE COMMUNE
# ============================================================

max_scaffold_length <- max(
  contigs$Scaffold_length,
  na.rm = TRUE
)

# Si tous les scaffolds mesurent au maximum 1 Mb,
# l'axe horizontal s'étend jusqu'à 1 Mb.
#
# Si au moins un scaffold dépasse 1 Mb,
# l'axe s'étend jusqu'à la longueur du scaffold le plus long.
max_len <- if (
  max_scaffold_length <= 1000000
) {
  1000000
} else {
  max_scaffold_length
}

cat(
  "Longueur maximale réelle des scaffolds :",
  max_scaffold_length,
  "pb\n"
)

cat(
  "Limite maximale utilisée pour l'axe des x :",
  max_len,
  "pb\n"
)


# ============================================================
# 12. POSITIONS DES CIF SUR LES SCAFFOLDS
# ============================================================

# Pour chaque identifiant numérique :
# - récupération des coordonnées disponibles pour CifA et CifB ;
# - recherche de la position minimale et maximale ;
# - calcul du milieu de l'ensemble CifA/CifB ;
# - détermination de la présence de CifA, CifB ou des deux.

cif_points <- d %>%
  group_by(
    Cif_Type,
    Souche,
    Accession
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    Cif_Type,
    Souche,
    Accession,
    matches(
      "Cif[AB]_(Start|End)[0-9]+_relative"
    )
  ) %>%
  pivot_longer(
    cols = -c(
      Cif_Type,
      Souche,
      Accession
    ),
    names_to = c(
      "cif",
      "coord",
      "id"
    ),
    names_pattern =
      "(Cif[AB])_(Start|End)([0-9]+)_relative",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  group_by(
    Cif_Type,
    Souche,
    Accession,
    id
  ) %>%
  summarise(
    min_pos = min(
      value,
      na.rm = TRUE
    ),
    max_pos = max(
      value,
      na.rm = TRUE
    ),
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

# Affichage des colonnes relatives aux cif pour vérification
print(
  grep(
    "Cif",
    names(d),
    value = TRUE
  )
)


# ============================================================
# 13. CRÉATION DES FENÊTRES DE 10 KB
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
    bin_mid = (
      bin_start + bin_end
    ) / 2,
    bin_length =
      bin_end - bin_start + 1
  )

# Vérification de la longueur des fenêtres
# La majorité doit mesurer 10 000 pb.
# Seule la dernière fenêtre de chaque scaffold peut être plus courte.
print(
  bins_by_souche %>%
    count(bin_length) %>%
    arrange(bin_length)
)


# ============================================================
# 14. FONCTION CALCULANT L'UNION DES INTERVALLES
# ============================================================

# Cette fonction fusionne les intervalles qui se chevauchent
# ou qui sont directement adjacents.
#
# Une position appartenant à plusieurs CDS n'est ainsi comptée
# qu'une seule fois.

union_length <- function(start, end) {

  # Retrait des coordonnées manquantes
  keep <- !is.na(start) & !is.na(end)

  start <- start[keep]
  end <- end[keep]

  # Aucun intervalle disponible
  if (length(start) == 0) {
    return(0)
  }

  # Sécurisation de l'ordre des coordonnées
  interval_start <- pmin(start, end)
  interval_end <- pmax(start, end)

  # Tri des intervalles selon leur position de départ,
  # puis selon leur position de fin
  interval_order <- order(
    interval_start,
    interval_end
  )

  interval_start <-
    interval_start[interval_order]

  interval_end <-
    interval_end[interval_order]

  # Initialisation avec le premier intervalle
  current_start <- interval_start[1]
  current_end <- interval_end[1]

  total_length <- 0

  # Fusion progressive des intervalles
  if (length(interval_start) > 1) {

    for (
      i in seq.int(
        from = 2,
        to = length(interval_start)
      )
    ) {

      # Chevauchement ou continuité directe
      if (
        interval_start[i] <=
          current_end + 1
      ) {

        current_end <- max(
          current_end,
          interval_end[i]
        )

      } else {

        # Ajout de l'intervalle fusionné précédent
        total_length <-
          total_length +
          (
            current_end -
              current_start +
              1
          )

        # Initialisation du nouvel intervalle
        current_start <-
          interval_start[i]

        current_end <-
          interval_end[i]
      }
    }
  }

  # Ajout du dernier intervalle fusionné
  total_length +
    (
      current_end -
        current_start +
        1
    )
}


# ============================================================
# 15. FONCTION DE CALCUL DE LA COUVERTURE UNIQUE
# ============================================================

# annotated_only = TRUE :
# couverture par les CDS non hypothétiques uniquement.
#
# annotated_only = FALSE :
# couverture par tous les CDS.

calculate_unique_coverage <- function(
    data,
    bins,
    annotated_only = FALSE
) {

  data_selected <- data

  # Pour la ligne A, seuls les CDS non hypothétiques
  # sont conservés.
  if (annotated_only) {

    data_selected <- data_selected %>%
      filter(hypothetical == 0)
  }

  # Association de chaque CDS à toutes les fenêtres
  # qu'il chevauche.
  intersections <- data_selected %>%
    select(
      Souche,
      Accession,
      CDS_start,
      CDS_end
    ) %>%
    inner_join(
      bins,
      by = c(
        "Souche",
        "Accession"
      ),
      relationship = "many-to-many"
    ) %>%
    filter(
      CDS_start <= bin_end,
      CDS_end >= bin_start
    ) %>%
    mutate(
      # Partie du CDS effectivement comprise
      # dans la fenêtre considérée
      Start_cut = pmax(
        CDS_start,
        bin_start
      ),
      Stop_cut = pmin(
        CDS_end,
        bin_end
      )
    )

  # Fusion des segments qui se chevauchent
  # à l'intérieur de chaque fenêtre.
  coverage <- intersections %>%
    group_by(
      Souche,
      Accession,
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
        "Souche",
        "Accession",
        "bin_start",
        "bin_end",
        "bin_mid",
        "bin_length"
      )
    ) %>%
    mutate(
      # Les fenêtres sans CDS reçoivent une couverture nulle
      covered_bp = replace_na(
        covered_bp,
        0
      ),

      # Proportion de positions couvertes par au moins un CDS
      coverage =
        covered_bp /
        bin_length
    ) %>%
    arrange(
      Souche,
      Accession,
      bin_start
    )

  return(coverage)
}


# ============================================================
# 16. LIGNE A : COUVERTURE PAR LES CDS NON HYPOTHÉTIQUES
# ============================================================

annotated_coverage <- calculate_unique_coverage(
  data = d,
  bins = bins_by_souche,
  annotated_only = TRUE
)

# La ligne A représente désormais :
# nombre de positions couvertes par au moins un CDS
# non hypothétique / nombre total de positions de la fenêtre.


# ============================================================
# 17. LIGNE B : COUVERTURE PAR TOUS LES CDS
# ============================================================

total_coding_coverage <- calculate_unique_coverage(
  data = d,
  bins = bins_by_souche,
  annotated_only = FALSE
)

# La ligne B représente désormais :
# nombre de positions couvertes par au moins un CDS
# / nombre total de positions de la fenêtre.
#
# Les CDS hypothétiques et non hypothétiques sont inclus.


# ============================================================
# 18. VÉRIFICATIONS DES COUVERTURES
# ============================================================

# Les valeurs doivent rester comprises entre 0 et 1.
invalid_annotated_coverage <-
  annotated_coverage %>%
  filter(
    coverage < 0 |
      coverage > 1
  )

invalid_total_coverage <-
  total_coding_coverage %>%
  filter(
    coverage < 0 |
      coverage > 1
  )

cat(
  "Nombre de valeurs invalides pour la ligne A :",
  nrow(invalid_annotated_coverage),
  "\n"
)

cat(
  "Nombre de valeurs invalides pour la ligne B :",
  nrow(invalid_total_coverage),
  "\n"
)

# Vérification des éventuelles fenêtres dupliquées
duplicated_annotated_bins <-
  annotated_coverage %>%
  count(
    Souche,
    Accession,
    bin_start
  ) %>%
  filter(n > 1)

duplicated_total_bins <-
  total_coding_coverage %>%
  count(
    Souche,
    Accession,
    bin_start
  ) %>%
  filter(n > 1)

cat(
  "Nombre de fenêtres dupliquées pour la ligne A :",
  nrow(duplicated_annotated_bins),
  "\n"
)

cat(
  "Nombre de fenêtres dupliquées pour la ligne B :",
  nrow(duplicated_total_bins),
  "\n"
)

# La couverture par les CDS non hypothétiques
# ne doit jamais dépasser la couverture totale.
coverage_comparison <- annotated_coverage %>%
  select(
    Souche,
    Accession,
    bin_start,
    annotated_coverage = coverage
  ) %>%
  left_join(
    total_coding_coverage %>%
      select(
        Souche,
        Accession,
        bin_start,
        total_coverage = coverage
      ),
    by = c(
      "Souche",
      "Accession",
      "bin_start"
    )
  ) %>%
  filter(
    annotated_coverage >
      total_coverage + 1e-12
  )

cat(
  paste0(
    "Nombre de fenêtres dans lesquelles la couverture ",
    "annotée dépasse la couverture totale : "
  ),
  nrow(coverage_comparison),
  "\n"
)


# ============================================================
# 19. FONCTION DE CRÉATION D'UN PANNEAU
# ============================================================

make_souche_panel <- function(
    souche_name,
    accession_name
) {

  panel_title <- paste0(
    souche_name,
    " | ",
    accession_name
  )

  # CDS de la combinaison Souche + Accession
  d_souche <- d %>%
    filter(
      Souche == souche_name,
      Accession == accession_name
    )

  # Couverture par les CDS non hypothétiques
  annotated_souche <-
    annotated_coverage %>%
    filter(
      Souche == souche_name,
      Accession == accession_name
    )

  # Couverture par tous les CDS
  coverage_souche <-
    total_coding_coverage %>%
    filter(
      Souche == souche_name,
      Accession == accession_name
    )

  # Positions des cif
  cif_souche <- cif_points %>%
    filter(
      Souche == souche_name,
      Accession == accession_name
    )

  # Informations sur le scaffold
  contig_souche <- contigs %>%
    filter(
      Souche == souche_name,
      Accession == accession_name
    )


  # ----------------------------------------------------------
  # LIGNE A
  #
  # Proportion de positions de chaque fenêtre de 10 kb
  # couvertes par au moins un CDS non hypothétique.
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
      limits = c(
        1,
        max_len
      ),
      expand = c(
        0,
        0
      )
    ) +
    scale_y_continuous(
      limits = c(
        0,
        ymax_a
      ),
      expand = c(
        0,
        0
      )
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
  #
  # Proportion de positions de chaque fenêtre de 10 kb
  # couvertes par au moins un CDS, quelle que soit
  # son annotation.
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
      limits = c(
        1,
        max_len
      ),
      expand = c(
        0,
        0
      )
    ) +
    scale_y_continuous(
      limits = c(
        0,
        ymax_b
      ),
      expand = c(
        0,
        0
      )
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
  #
  # Distribution exacte des CDS et position des cif
  # le long du scaffold.
  # ----------------------------------------------------------

  p_scaffold <- ggplot() +

    # Trait noir correspondant à la longueur totale du scaffold
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

    # CDS non hypothétiques en rouge
    # CDS hypothétiques en gris
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

    # Triangles indiquant les positions des cif
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
      values = c(
        "0" = "red",
        "1" = "grey25"
      )
    ) +

    scale_fill_manual(
      values = c(
        "CifA + CifB" = "white",
        "CifB seul" = "grey20",
        "CifA seul" = "green"
      )
    ) +

    scale_x_continuous(
      limits = c(
        1,
        max_len
      ),
      expand = c(
        0,
        0
      )
    ) +

    scale_y_continuous(
      limits = c(
        0.7,
        1.3
      ),
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


  # Assemblage vertical des trois lignes
  p_annotated /
    p_coverage /
    p_scaffold +
    plot_layout(
      heights = c(
        1,
        1,
        0.5
      )
    )
}


# ============================================================
# 20. LISTE DES COMBINAISONS SOUCHE + ACCESSION
# ============================================================

souches_accessions <- contigs %>%
  select(
    Souche,
    Accession
  ) %>%
  distinct()


# ============================================================
# 21. CRÉATION DES PANNEAUX
# ============================================================

plots_by_souche <- lapply(
  seq_len(
    nrow(souches_accessions)
  ),
  function(i) {

    make_souche_panel(
      souche_name =
        souches_accessions$Souche[i],

      accession_name =
        souches_accessions$Accession[i]
    )
  }
)


# ============================================================
# 22. ASSEMBLAGE DE LA FIGURE FINALE
# ============================================================

figure_finale <- wrap_plots(
  plots_by_souche,
  ncol = 1
)

figure_finale


# ============================================================
# 23. VALEURS MAXIMALES RÉELLES
# ============================================================

cat(
  paste0(
    "Couverture maximale par les CDS non hypothétiques ",
    "dans le graphique A : "
  ),
  max(
    annotated_coverage$coverage,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  paste0(
    "Couverture maximale par l'ensemble des CDS ",
    "dans le graphique B : "
  ),
  max(
    total_coding_coverage$coverage,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 24. DIMENSIONS CONSEILLÉES POUR L'EXPORT
# ============================================================

# Largeur : 14 pouces
# Hauteur : 2,5 pouces par combinaison Souche + Accession

figure_width <- 14

figure_height <-
  nrow(souches_accessions) * 2.5

cat(
  "Largeur conseillée :",
  figure_width,
  "pouces\n"
)

cat(
  "Hauteur conseillée :",
  figure_height,
  "pouces\n"
)

# Export facultatif
#
# ggsave(
#   filename = "figure_finale_cif_like.png",
#   plot = figure_finale,
#   width = figure_width,
#   height = figure_height,
#   units = "in",
#   dpi = 300,
#   limitsize = FALSE
# )

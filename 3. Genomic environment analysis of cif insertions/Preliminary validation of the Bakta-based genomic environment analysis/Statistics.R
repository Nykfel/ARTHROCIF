#Unfinished

library(dplyr)
library(ggplot2)
library(patchwork)
library(betareg)
library(lmtest)
library(coin)

#========================
# 1. Import des données
#========================

info <- read.delim("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/Genomes_references/info.txt")
info$Sequence <- sub("/.*", "", info$Sequence)
info <- info[info$Clade != "Insert", ]
info$Clade <- factor(info$Clade)

# Conversion éventuelle en proportions [0,1]
if (max(info$GC, na.rm = TRUE) > 1) {
  info$GC <- info$GC / 100
}

if (max(info$coding_density, na.rm = TRUE) > 1) {
  info$coding_density <- info$coding_density / 100
}

# Proportions à l'échelle de la souche
info <- info %>%
  mutate(
    prop_hypothetical = hypotheticals / CDSs,
    prop_annotated = (CDSs - hypotheticals) / CDSs
  )

# Vérifications rapides
range(info$GC, na.rm = TRUE)
range(info$coding_density, na.rm = TRUE)
range(info$prop_annotated, na.rm = TRUE)
table(info$Clade)

#========================
# 2. Fonctions utilitaires
#========================

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  if (p < 0.001) return(paste0("p = ", format(p, scientific = TRUE, digits = 2)))
  paste0("p = ", sprintf("%.3f", p))
}

make_letters_two_groups <- function(data, response, pval, alpha = 0.05) {
  
  tmp <- data %>%
    group_by(Clade) %>%
    summarise(
      centre = median(.data[[response]], na.rm = TRUE),
      .groups = "drop"
    )
  
  if (pval < alpha) {
    tmp <- tmp %>%
      arrange(desc(centre)) %>%
      mutate(letter = c("a", "b")[seq_len(n())]) %>%
      arrange(Clade)
  } else {
    tmp <- tmp %>%
      mutate(letter = "a")
  }
  
  tmp$y <- 0.92
  tmp
}

make_panel <- function(data, response, ylab, pval, title_panel, colors) {
  
  lab_df <- make_letters_two_groups(
    data = data,
    response = response,
    pval = pval
  )
  
  ggplot(data, aes(x = Clade, y = .data[[response]], fill = Clade)) +
    geom_boxplot(
      width = 0.55,
      alpha = 0.75,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.2,
      alpha = 0.9,
      color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3.2,
      fill = "white",
      color = "black",
      stroke = 0.7
    ) +
    geom_text(
      data = lab_df,
      aes(x = Clade, y = y, label = letter),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 5
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    labs(
      title = title_panel,
      x = NULL,
      y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

#========================
# 3. Tests statistiques
#========================

# A. GC ~ Clade : Wilcoxon exact
test_gc <- wilcox_test(GC ~ Clade, data = info, distribution = exact())
p_gc <- pvalue(test_gc)

# B. Coding density ~ Clade : Wilcoxon exact
# La régression bêta a été écartée : le paramètre de précision φ diffère
# fortement entre modèle nul (φ = 2.36) et modèle complet (φ = 13.13),
# ce qui invalide le test du rapport de vraisemblance.
test_cd <- wilcox_test(coding_density ~ Clade, data = info, distribution = exact())
p_cd <- pvalue(test_cd)

# C. Proportion de CDS annotés ~ Clade : Wilcoxon exact
test_annot <- wilcox_test(prop_annotated ~ Clade, data = info, distribution = exact())
p_annot <- pvalue(test_annot)

# Affichage des tests
test_gc
test_cd
test_annot

#========================
# 4. Palette de couleurs
#========================

clade_levels <- levels(info$Clade)

palette_clade <- setNames(
  c("#4C78A8", "#F58518")[seq_along(clade_levels)],
  clade_levels
)

#========================
# 5. Création des panneaux
#========================

p1 <- make_panel(
  data = info,
  response = "GC",
  ylab = "Taux de GC",
  pval = p_gc,
  title_panel = "GC",
  colors = palette_clade
)

p2 <- make_panel(
  data = info,
  response = "coding_density",
  ylab = "Densité codante",
  pval = p_cd,
  title_panel = "Densité codante",
  colors = palette_clade
)

p3 <- make_panel(
  data = info,
  response = "prop_annotated",
  ylab = "Proportion de CDS annotés",
  pval = p_annot,
  title_panel = "Proportion de CDS annotés",
  colors = palette_clade
)

#========================
# 6. Figure finale
#========================

figure_finale <- (p1 | p2 | p3) +
  plot_annotation(
    tag_levels = "A")
#caption = "Points = souches individuelles ; carré blanc = moyenne ; boxplots = distribution par clade.\nDes lettres différentes indiquent une différence significative entre clades au seuil alpha = 0,05."


figure_finale

#========================
# 7. Export
#========================

ggsave(
  filename = "Figure_clade_resume_tests.pdf",
  plot = figure_finale,
  width = 12,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = "Figure_clade_resume_tests.png",
  plot = figure_finale,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 600
)

#========================
# Résumé descriptif par clade
#========================

resume_stats <- info %>%
  group_by(Clade) %>%
  summarise(
    n = n(),
    
    GC_moyenne = mean(GC, na.rm = TRUE),
    GC_ecart_type = sd(GC, na.rm = TRUE),
    
    densite_codante_moyenne = mean(coding_density, na.rm = TRUE),
    densite_codante_ecart_type = sd(coding_density, na.rm = TRUE),
    
    prop_CDS_annotes_moyenne = mean(prop_annotated, na.rm = TRUE),
    prop_CDS_annotes_ecart_type = sd(prop_annotated, na.rm = TRUE),
    
    .groups = "drop"
  )

resume_stats


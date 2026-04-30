library(dplyr)
library(ggplot2)
library(patchwork)

par(mfrow = c(1,1))
d <- read.delim2("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/Tentative/df.txt")

Cif_pos <- read.delim("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/Tentative/Cifs_pos.txt")
library(dplyr)

#permet de retrovuer les vrai positions des cifs sur la portion retenue
Cif_pos <- Cif_pos %>%
  mutate(across(
    .cols = starts_with("Cif"),
    .fns = ~ . - (Marker_start - 1),
    .names = "{.col}_relative"
  ))

d <- d %>%
  left_join(Cif_pos, by = c("Souche" = "Bakta"))

#d <- d[d$Clade != "HGT", ]
#d <- d[d$Souche != "BFE", ]

#Permet de séparer les CDSs en deux catégories, ceux reconnus et les hypotheticals
d$hypothetical <- ifelse(
  d$Product %in% c(
    "hypothetical protein",
    "Conserved hypothetical exported protein",
    "(pseudo) hypothetical protein"
  ),
  1, 0
)

d <- d[!d$Type %in% c(
  "ncRNA-region", "tRNA", "sorf", "ncRNA", "tmRNA",
  "rRNA", "crispr", "crispr-repeat", "crispr-spacer",
  "assembly_gap", "oriC"
), ]


# Sert à vérifier ma méthode de calcul, résultats concordent avec le reste

density_df <- d %>%
  group_by(Souche, Scaffold_length) %>%
  summarise(
    total_cds_length = sum(Stop - Start + 1),
    coding_density = total_cds_length / unique(Scaffold_length),
    .groups = "drop"
  )

boxplot(density_df$coding_density ~ density_df$Souche)

##############
bin_size <- 10000
max_len <- max(d$Scaffold_length)

d2 <- d %>%
  mutate(
    bin = floor((Start - 1) / bin_size) * bin_size + 1,
    bin_mid = bin + bin_size / 2
  )

gene_density <- d2 %>%
  group_by(Souche, bin_mid) %>%
  summarise(n_CDS = n(), .groups = "drop")

contigs <- d %>%
  distinct(Souche, Scaffold_length)

p_density <- ggplot(gene_density,
                    aes(x = bin_mid, y = n_CDS)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Souche, ncol = 1)+
  scale_x_continuous(limits = c(1, max_len)) +
  labs(
    title = "Densité de CDS par fenêtre de 10 kb",
    x = "Position sur le scaffold",
    y = "Nombre de CDS"
  ) +
  theme_minimal()

p_density

## avec cifs d'indiqués 

cif_points <- bind_rows(
  d %>% select(Souche, pos = CifA_Start1_relative) %>% mutate(type = "CifA"),
  d %>% select(Souche, pos = CifB_Start1_relative) %>% mutate(type = "CifB"),
  d %>% select(Souche, pos = CifA_Start2_relative) %>% mutate(type = "CifA"),
  d %>% select(Souche, pos = CifB_Start2_relative) %>% mutate(type = "CifB"),
  d %>% select(Souche, pos = CifA_Start3_relative) %>% mutate(type = "CifA"),
  d %>% select(Souche, pos = CifB_Start3_relative) %>% mutate(type = "CifB")
) %>%
  filter(!is.na(pos))

ggplot(gene_density,
                    aes(x = bin_mid, y = n_CDS)) +
  
  geom_line(linewidth = 0.8) +
  
  geom_point(
    data = cif_points,
    aes(x = pos, y = Inf, fill = type),
    shape = 25,       # triangle inversé
    size = 3,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = c("CifA" = "pink", "CifB" = "red"),
    name = "Gènes cif"
  ) +
  
  facet_wrap(~ Souche, ncol = 1) +
  scale_x_continuous(limits = c(1, max_len)) +
  labs(
    title = "Densité de CDS par fenêtre de 10 kb",
    x = "Position relative sur le scaffold",
    y = "Nombre de CDS"
  ) +
  theme_minimal()


#### N CDSs annotés

gene_density_hypo <- d2 %>%
  filter(hypothetical == 0) %>%
  group_by(Souche, bin_mid) %>%
  summarise(n_CDS = n(), .groups = "drop")

density_Cds_annotated = ggplot(gene_density_hypo,
                               aes(x = bin_mid, y = n_CDS)) +
  geom_area(fill = "red", alpha = 0.5) +   # ← surface
  geom_line(linewidth = 0.8, color = "red") +
  facet_wrap(~ Souche, ncol = 1) +
  scale_x_continuous(limits = c(1, max_len)) +
  labs(
    title = "Densité de CDS bactériens par fenêtre de 50 kb",
    x = "Position sur le scaffold",
    y = "Nombre de CDS hypothétiques"
  ) +
  theme_minimal()

density_Cds_annotated

## avec cifs
density_Cds_annotated <- ggplot(gene_density_hypo,
                                aes(x = bin_mid, y = n_CDS)) +
  geom_area(fill = "red", alpha = 0.5) +
  geom_line(linewidth = 0.8, color = "red") +
  
  geom_point(
    data = cif_points,
    aes(x = pos, y = Inf, fill = type),
    shape = 25,
    size = 3,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = c("CifA" = "pink", "CifB" = "red"),
    name = "Gènes cif"
  ) +
  
  facet_wrap(~ Souche, ncol = 1) +
  scale_x_continuous(limits = c(1, max_len)) +
  labs(
    title = "Densité de CDS bactériens par fenêtre de 50 kb",
    x = "Position relative sur le scaffold",
    y = "Nombre de CDS annotés"
  ) +
  theme_minimal()

density_Cds_annotated

#####

p_scaffold <- ggplot() +
  geom_segment(data = contigs,
               aes(x = 1, xend = Scaffold_length,
                   y = Souche, yend = Souche),
               linewidth = 1, color = "black") +
  geom_segment(data = d,
               aes(x = Start, xend = Stop,
                   y = Souche, yend = Souche,
                   color = factor(hypothetical)),
               linewidth = 3) +
  scale_color_manual(
    values = c("0" = "red", "1" = "grey25"),
    labels = c("0" = "CDS annoté", "1" = "protéine hypothétique"),
    name = "Type de CDS"
  ) +
  scale_x_continuous(limits = c(1, max_len)) +
  labs(
    title = "Répartition des CDS annotés et hypothétiques sur les scaffolds",
    x = "Position sur le scaffold",
    y = "Souche"
  ) +
  theme_minimal()

## AVEC CIFS

p_scaffold <- ggplot() +
  geom_segment(data = contigs,
               aes(x = 1, xend = Scaffold_length,
                   y = Souche, yend = Souche),
               linewidth = 1, color = "black") +
  
  geom_segment(data = d,
               aes(x = Start, xend = Stop,
                   y = Souche, yend = Souche,
                   color = factor(hypothetical)),
               linewidth = 3) +
  
  geom_point(
    data = cif_points,
    aes(x = pos, y = Souche, fill = type),
    shape = 25,
    size = 3,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = c("0" = "red", "1" = "grey25"),
    labels = c("0" = "CDS annoté", "1" = "protéine hypothétique"),
    name = "Type de CDS"
  ) +
  
  scale_fill_manual(
    values = c("CifA" = "pink", "CifB" = "red"),
    name = "Gènes cif"
  ) +
  
  scale_x_continuous(limits = c(1, max_len)) +
  labs(
    title = "Répartition des CDS et positions relatives des gènes cif",
    x = "Position relative sur le scaffold",
    y = "Souche"
  ) +
  theme_minimal()

p_scaffold

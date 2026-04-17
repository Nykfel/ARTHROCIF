library(dplyr)
library(ggplot2)

dom_cols <- c(
  "AAA_ATPase_like"             = "#2CA25F",
  "PD-(D/E)XK_nuclease"         = "#D73027",
  "ULP1"                        = "#8ECDF4",
  "OTU-like_cysteine_protease"  = "orange",
  "DUF3491"                     = "gray65",
  "RTX_toxin"                   = "#003366",
  "Ankyrin_repeat"              = "#40E0D0",
  "Pore_forming_toxin_TcdA/B"   = "pink",
  "Latrotoxin"                  = "#D12A83",
  "Unknown"                     = "gray90",
  "RNA-binding-like"            = "purple",
  "Apoptosis_regulator-like"    = "#FFD814",
  "Salivary-gland_toxin"        = "#365729",
  "TPR_repeats"                 = "#513D73"
)

plot_domain_architecture <- function(df,
                                     order_sequences = NULL,
                                     box_halfheight  = 0.28,
                                     mark_halfheight = 0.35,
                                     show_legend     = TRUE) {
  
  needed <- c("sequence", "taille", "Domain", "START", "END")
  missing_needed <- setdiff(needed, names(df))
  
  if (length(missing_needed) > 0) {
    stop("Il manque les colonnes obligatoires : ",
         paste(missing_needed, collapse = ", "))
  }
  
  df2 <- df %>%
    mutate(
      sequence = as.character(sequence),
      Domain   = ifelse(tolower(as.character(Domain)) == "unknown", "Unknown", as.character(Domain)),
      START    = as.numeric(START),
      END      = as.numeric(END),
      taille   = as.numeric(taille)
    ) %>%
    mutate(
      length_domaine = END - START,
      ordre_seq = match(sequence, unique(sequence))
    ) %>%
    arrange(ordre_seq, sequence, desc(length_domaine)) %>%
    select(-length_domaine, -ordre_seq)
  
  proteins <- df2 %>%
    group_by(sequence) %>%
    summarise(
      length = {
        x <- suppressWarnings(max(taille, na.rm = TRUE))
        y <- max(END, na.rm = TRUE)
        ifelse(is.infinite(x) | is.na(x), y, x)
      },
      .groups = "drop"
    )
  
  seq_levels <- if (!is.null(order_sequences)) {
    rev(order_sequences)
  } else if ("Ordre" %in% names(df2)) {
    df2 %>%
      distinct(sequence, Ordre) %>%
      arrange(Ordre) %>%
      pull(sequence) %>%
      rev()
  } else {
    rev(unique(df2$sequence))
  }
  
  proteins <- proteins %>%
    mutate(protein = factor(sequence, levels = seq_levels))
  
  backbone <- proteins %>%
    transmute(
      protein,
      x = 1,
      xend = length,
      y = as.numeric(protein),
      yend = as.numeric(protein)
    )
  
  df_domains <- df2 %>% filter(!Domain %in% c("STOP", "START"))
  df_marks   <- df2 %>% filter(Domain %in% c("STOP", "START"))
  
  domains_plot <- df_domains %>%
    transmute(
      protein = factor(sequence, levels = seq_levels),
      domain  = Domain,
      xmin    = START,
      xmax    = END
    ) %>%
    mutate(
      y    = as.numeric(protein),
      ymin = y - box_halfheight,
      ymax = y + box_halfheight
    )
  
  marks_plot <- df_marks %>%
    transmute(
      protein = factor(sequence, levels = seq_levels),
      x       = START,
      mark    = factor(Domain, levels = c("START", "STOP"))
    ) %>%
    mutate(
      y  = as.numeric(protein),
      y0 = y - mark_halfheight,
      y1 = y + mark_halfheight
    )
  
  pal <- c(
    dom_cols,
    setNames(
      rep("grey80", length(setdiff(unique(domains_plot$domain), names(dom_cols)))),
      setdiff(unique(domains_plot$domain), names(dom_cols))
    )
  )
  
  ggplot() +
    geom_segment(
      data = backbone,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.2, color = "grey40", lineend = "round"
    ) +
    geom_rect(
      data = domains_plot,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = domain),
      color = NA, linewidth = 0.3
    ) +
    geom_segment(
      data = marks_plot,
      aes(x = x, xend = x, y = y0, yend = y1, color = mark),
      linewidth = 0.25, lineend = "butt"
    ) +
    scale_fill_manual(name = "Domain", values = pal, drop = FALSE) +
    scale_color_manual(
      name   = "Domain",
      values = c("START" = "#006400", "STOP" = "black"),
      limits = c("START", "STOP"),
      drop   = FALSE,
      guide  = guide_legend(override.aes = list(linewidth = 1.2))
    ) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 1)) +
    scale_y_continuous(
      breaks = seq_along(seq_levels),
      labels = seq_levels,
      expand = expansion(mult = c(0.12, 0.12))
    ) +
    labs(x = "Position (aa)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 8),
      legend.position    = if (show_legend) "right" else "none",
      plot.margin        = margin(5.5, 5.5, 5.5, 20)
    )
}

p <- plot_domain_architecture(df)
p

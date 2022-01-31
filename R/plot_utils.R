my_pretty_theme <- theme_minimal(base_family = "Lato", base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          # Bold, bigger title
          #plot.title = element_text(face = "bold", size = rel(1.7)),
          # Plain, slightly bigger subtitle that is grey
          plot.subtitle = element_text(face = "plain", size = rel(1.3), color = "grey70"),
          # Italic, smaller, grey caption that is left-aligned
          plot.caption = element_text(face = "italic", size = rel(0.7), 
                                      color = "grey70", hjust = 0),
          # Bold legend titles
          legend.title = element_text(face = "bold"),
          # Bold, slightly larger facet titles that are left-aligned for the sake of repetition
          strip.text = element_text(face = "bold", size = rel(1.1), hjust = 0),
          # Add a light grey background to the facet titles, with no borders
          strip.background = element_rect(fill = "grey90", color = NA),
          # Add a thin grey border around all the plots to tie in the facet titles
          panel.border = element_rect(color = "grey90", fill = NA))

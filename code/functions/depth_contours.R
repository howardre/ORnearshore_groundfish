## Depth distribution by year

depth_contours <- function(subset, year_range, title){
  subset %>% filter(year %in% (year_range)) %>%
    filter(CPUE > 0) %>%
    ggplot(aes(x = depth, y = lat)) +
    stat_density_2d(aes(fill = stat(level)),
                    geom = "polygon",
                    show.legend = F,
                    bins = 11) +
    scale_fill_viridis_c(option = "F", direction = -1) +
    scale_x_continuous(expand = c(0, 0), limits = c(-210, 10)) +
    scale_y_continuous(expand = c(0, 0), limits = c(40, 49)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey84"),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 9),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          text = element_text(family = "serif")) +
    labs(x = " ", y = " ", title = title) +
    coord_cartesian(ylim = c(42, 47), xlim = c(-200, 0))
}

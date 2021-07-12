## Depth distribution by year

depth_contours <- function(subset, year_upper, color, title){
  subset %>% filter(year <= year_upper) %>% filter(CPUE > 0) %>%
    ggplot(aes(x = depth, y = lat)) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon", show.legend = F) +
    scale_fill_distiller(palette = color, direction = -1) +
    scale_x_continuous(expand = c(0, 0), limits = c(-200, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(40, 49)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey84"),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          text = element_text(family = "serif")) +
    labs(x = "Depth (m)", y = "Latitude", title = title) +
    coord_cartesian(ylim = c(42, 47), xlim = c(-200, 0))
}

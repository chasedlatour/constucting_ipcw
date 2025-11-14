pretty plot coord_equal(
  +
    labs(y = "Cumulative Incidence",
         x = "Years since December 6, 1995",
         color = NULL) +
    scale_x_continuous(breaks = seq(0, end_time*365, by = 365),
                       labels = seq(0, end_time, by = 1), 
                       expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 0.35), 
                       expand = c(0,0)) +
    scale_color_grey(start = 0, end = 0.6,
                     labels = c("IPCW", "Naive")) +
    theme_classic(base_size = 16) + 
    theme(legend.position = "inside",
          legend.justification = c(0.2, 0.8),
          legend.background = element_rect(
            color = "black", 
            fill = "white",
            linewidth = 0.5,
            linetype = "solid"
          ))
  
)
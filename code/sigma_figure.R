
library(geomtextpath)

sigma <- 1.5

d <- tibble(x = seq(4, 16, .1),
            z = (x - 10) / sigma,
            d = dnorm(z) / dnorm(0))

p <- ggplot(d, aes(x, d)) +
      geom_area(fill = "gray70") +
      
      geom_textpath(label = "historic frequency distribution across species CA range", hjust = .5, vjust = 1.5) +
      
      annotate(geom = "line", x = c(5.5, 7), y = .75, color = "red",
               arrow = arrow(angle = 20, ends = "both", type = "closed", length = unit(.1, "in"))) +
      annotate(geom = "text", x = 6.25, y = .75, color = "red", vjust = -.5,
               label = "sigma = 1") +
      
      annotate(geom = "line", x = c(11.5, 11.5 + sigma*3), y = .75, color = "red",
               arrow = arrow(angle = 20, ends = "both", type = "closed", length = unit(.1, "in"))) +
      annotate(geom = "text", x = 11.5 + sigma*3/2, y = .75, color = "red", vjust = -.5,
               label = "sigma = 3") +
      
      scale_x_continuous(
            "e.g. winter mininum temperature (deg C)", 
            breaks = c(4, 7, 10, 13, 16),
            sec.axis = sec_axis(~ (. - 10) / 1.5, breaks = -4:4,
                                name = "\nstandard deviations with respect to species CA range")) +
      theme_minimal() +
      theme(panel.grid.major.x = element_line(color = "black"),
            panel.grid.minor.x = element_line(color = "black"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank())

ggsave("www/sigma.jpg", p, width = 8, height = 3, units = "in")

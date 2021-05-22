library(tidyverse)
theme_set(theme_bw())

### plotting growth for each plant over time (not standardized, but probably fine for this since they are split by plant, and more intuitive being just cm units)
growth_tab <- read.table("data/all-plant-average-heights.tsv", sep = "\t", header = TRUE, na.strings="")
growth_for_plotting <- growth_tab %>% pivot_longer(!Days, names_to = "Plant", values_to = "Height (cm)")

# split
ggplot(growth_for_plotting, aes(x = Days, y = `Height (cm)`, color = Plant, group = Plant)) + geom_point() + geom_smooth(aes(fill = Plant), alpha = 0.1) +
    facet_wrap(~Plant, nrow=3, ncol=1, scales = "free_y") + theme(legend.position = "none") + ggtitle("Height over time")

# together (probably fine together)
ggplot(growth_for_plotting, aes(x = Days, y = `Height (cm)`, color = Plant, group = Plant)) + geom_point() + geom_smooth(aes(fill = Plant), alpha = 0.1) +
    theme(legend.position = "bottom") + ggtitle("Height over time")


tank_tab <- read.table("data/all-tanks-measured-data.tsv", sep = "\t", header = TRUE, na.strings = "", stringsAsFactors=FALSE)

tank_tab$Tank <- factor(tank_tab$Tank)
tank_tab$Ammonia <- as.numeric(tank_tab$Ammonia)
tank_tab$Nitrite <- as.numeric(tank_tab$Nitrite)
tank_tab$Nitrate <- as.numeric(tank_tab$Nitrate)
tank_tab$pH <- as.numeric(tank_tab$pH)
tank_tab$pH <- as.numeric(tank_tab$pH)

tank_tab_for_plotting <- tank_tab %>% pivot_longer(cols = c("Ammonia", "Nitrite", "Nitrate", "pH", "Temp"), names_to = "Measure")
tank_tab_for_plotting$Measure <- factor(tank_tab_for_plotting$Measure)
str(tank_tab_for_plotting)

ggplot(tank_tab_for_plotting, aes(x = Day, y = value, color = Measure)) +
  geom_point() + geom_smooth(aes(fill = Measure), alpha = 0.25) + facet_grid(Measure ~ Tank, scales="free_y") +
  guides(color = FALSE, fill = FALSE) +
    labs(x = "Days", y = "Value")


### by water (values standardized by plant type)
water_tab <- read.table("data/std-height-averages-by-water.tsv", sep = "\t", header = TRUE, na.strings="")
water_for_plotting <- water_tab %>% pivot_longer(!Days, names_to = "Water source", values_to = "Height (z-score)")

ggplot(water_for_plotting, aes(x = Days, y = `Height (z-score)`, color = `Water source`, group = `Water source`)) + geom_point() + geom_smooth(aes(fill = `Water source`), alpha = 0.1) +
    theme(legend.position = "bottom") + ggtitle("Plant-height over time")

### by lighting (values standardized by plant type)
lighting_tab <- read.table("data/std-height-averages-by-lighting.tsv", sep = "\t", header = TRUE, na.strings="")
lighting_for_plotting <- lighting_tab %>% pivot_longer(!Days, names_to = "Lighting source", values_to = "Height (z-score)")

ggplot(lighting_for_plotting, aes(x = Days, y = `Height (z-score)`, color = `Lighting source`, group = `Lighting source`)) + geom_point() + geom_smooth(aes(fill = `Lighting source`), alpha = 0.1) +
    theme(legend.position = "bottom") + ggtitle("Plant-height over time")


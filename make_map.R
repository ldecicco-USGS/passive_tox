library(tidyverse)
library(toxEval)
library(maptools)
library(maps)
library(nhdplusTools)
library(sf)

path_to_data <- Sys.getenv("PASSIVE_PATH")

options(dplyr.summarise.inform = FALSE)
source(file = "read_chemicalSummary.R")
source(file = "R/report/combo_plot2.R")
source(file = "R/report/stack_plots.R")

site_info <- tox_list$chem_site

site_info <- site_info %>% 
  bind_rows(filter(site_info, SiteID == "04231600"))

site_info$SiteID[duplicated(site_info$SiteID)] <- "04232007"

lake <- st_read("Great_Lakes-shp/Great_Lakes.shp")

glri_states <- map_data("state") %>% 
  filter(region %in% c("wisconsin", "illinois", "michigan",
                       "ohio","indiana", "minnesota",
                       "pennsylvania", "new york"))

site_info <- site_info %>% 
  filter(!is.na(dec_lat),
         !is.na(dec_lon))

p <- sf::st_as_sf(data.frame(lon = site_info$dec_lon,
                             lat = site_info$dec_lat), 
                  coords = c("lon", "lat"), crs = 4326)


# bbox <- sf::st_bbox(c(xmin = min(site_info$dec_lon, na.rm = TRUE)-3, ymin = min(site_info$dec_lat, na.rm = TRUE)-3,
#                       xmax = max(site_info$dec_lon, na.rm = TRUE)+3, ymax = max(site_info$dec_lat, na.rm = TRUE)+3),
#                     crs = "+proj=longlat +datum=WGS84 +no_defs")

bbox <- sf::st_bbox(c(xmin = min(glri_states$long, na.rm = TRUE), ymin = min(glri_states$lat, na.rm = TRUE),
                      xmax = max(glri_states$long, na.rm = TRUE), ymax = max(glri_states$lat, na.rm = TRUE)),
                    crs = "+proj=longlat +datum=WGS84 +no_defs")

data <- plot_nhdplus(bbox = bbox, streamorder = 4, actually_plot = FALSE)

# saveRDS(data, file = "plot_nhd_out_3.rds")
data <- readRDS("plot_nhd_out_2.rds")

flowlines <- data$flowline

indexes <- get_flowline_index(flowlines,
                              p,
                              search_radius = 0.1,
                              max_matches = 1)

indexes <- left_join(sf::st_sf(id = c(1:nrow(p)),
                               geom = sf::st_geometry(p)),
                     indexes, by = "id")

flowlines_small <- flowlines %>% 
  filter(REACHCODE %in% indexes$REACHCODE)

flowlines_med <- flowlines %>% 
  filter(TerminalPa %in% unique(flowlines_small$TerminalPa))

site_info <- prep_site_list(site_info)

graphData <- graph_chem_data(chemical_summary = chemicalSummary)

gd <- graphData %>% 
  group_by(site) %>% 
  summarize(maxmax = sum(meanEAR)) %>% 
  left_join(site_info, by = c("site"="SiteID")) %>% 
  arrange(`Short Name`)

gdplot <- ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group),
               data = glri_states, fill = "grey90",
               color = "white") +
  geom_sf(data = lake, fill = "lightblue") +
  geom_point(data = gd,
             aes(x = dec_lon, y = dec_lat,
                 fill = maxmax),
             size = 3, pch = 21, color = "black") +
  scale_fill_gradient(low = "grey", high = "red") +
  # scale_fill_gradient(low = "white", high = "red") +
  geom_sf(data = flowlines_med , color = "blue") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        plot.background = element_rect(fill = "transparent",colour = NA)
        )

gdplot
ggsave(gdplot, file = "map3.png", bg = "transparent")

library(grid)

color_map <- class_colors(tox_list)

gd_stack <- graphData %>% 
  left_join(site_info, by = c("site"="SiteID"))

top4 <- gd_stack %>% 
  group_by(Class) %>% 
  summarise(sumEARs = sum(meanEAR)) %>% 
  arrange(desc(sumEARs))

gd_stack$Class <- as.character(gd_stack$Class)
gd_stack$Class[!gd_stack$Class %in% as.character(top4$Class[1:5])] <- "Other (14)"

gd_stack$Class <- factor(gd_stack$Class, 
                         levels = c(as.character(top4$Class[1:5]), "Other (14)"))


names(color_map)[names(color_map) == "Other"] <- "Other (14)"

stack1 <- ggplot() +
  geom_col(data = gd_stack, 
           aes(x=site, y=meanEAR, fill = Class),
           position = position_stack())  +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ylab("Exposure-Activity Ratio (EAR)") +
  facet_grid(. ~ site_grouping , switch = "x",
             scales="free", space="free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values = color_map, drop=TRUE) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        # legend.position = c(.95, .95),
        # legend.justification = c("right", "top"),
        # legend.box.just = "right",
        # legend.margin = margin(6, 6, 6, 6),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
stack1

ggsave(stack1, file = "stack1.png", bg = "transparent")

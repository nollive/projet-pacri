library(ggplot2)
library(gganimate)
library(gifski)
library(dplyr)
library(ggrepel)
library(tidyr)



#loc_path <- file.path("dev-localization-nodscov2.RData")
wd <- getwd()
loc_path <- file.path(wd, "..", "..", "out", "loc-nodscov2","dev-localization-nodscov2.RData")
load(loc_path)


new_begin_date <- begin_date + (5*60*60) ## OFFSET TO START THE DAY AT 12AM
new_end_date <- new_begin_date + (24*60*60)
print(new_end_date - new_begin_date)
t_begin <- 5*60*2 ## OFFSET t
t_end <- (t_begin + 24*60*2) - 1


distribute_points <- function(n) {
  theta <- seq(0, 2 * pi, length.out = n + 1)[-1]
  radius <- 0.2  # Radius of the circle where the points will be plotted
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  return(data.frame(offset_x = x, offset_y = y))
}


rooms_coords <- rooms %>%
  mutate(
    x = case_when(
      grepl("001", id) ~ id_room,
      grepl("Corridor", room) ~ 1,
      room == "Office" ~ 1,
      room == "Nursing station" ~ 3,
      room == "Medical Restroom"~ 1,
      room == "Paramedical Restroom"~ 3,
      TRUE ~ NA_real_
    ),
    y = case_when(
      grepl("001", id) ~ 1,
      grepl("Corridor", room) ~ 2,
      room == "Office" ~ 3,
      room == "Nursing station" ~ 3,
      grepl("Restroom", room) ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(localization = id_room) %>%
  distinct(id_room, .keep_all = TRUE) %>%
  select(localization, x, y)


data <- do.call(rbind, global_localization) %>% filter(between(time, t_begin, t_end))
data <- do.call(rbind, global_localization) %>% filter(time < 50)
data <- data %>%
  filter(localization != -1) %>%
  left_join(rooms_coords, by = c("localization")) %>%
  left_join(admission, by = c("id")) %>%
  mutate(cat = ifelse(is.na(cat), "patient", cat)) %>%
  filter(cat != "patient")



max_n <- data %>%
  group_by(localization, time) %>%
  summarise(n = n(), .groups = 'drop') %>%
  summarise(max_n = max(n)) %>%
  pull(max_n)



apply_offsets <- function(group_by_data) {
  n <- nrow(group_by_data)
  offsets <- distribute_points(n)
  group_by_data <- group_by_data %>%
    mutate(offset_x = offsets$offset_x, offset_y = offsets$offset_y)
  return(group_by_data)
}


data <- data %>%
  group_by(localization, time) %>%
  group_modify(~ apply_offsets(.x)) %>%
  ungroup() %>%
  mutate(x_adj = x + offset_x, y_adj = y + offset_y)



segments_data <- data %>%
  group_by(id) %>%
  arrange(time) %>%
  mutate(x_end = lead(x_adj),
         y_end = lead(y_adj)) %>%
  filter(!is.na(x_end) & !is.na(y_end)) %>%
  ungroup()


# Plot
p <- ggplot(data %>% filter(time > 10), aes(x = x_adj, y = y_adj, color = cat, group = id)) +
  geom_point(size = 5) +
  #geom_segment(data = segments_data, aes(xend = x_end, yend = y_end), alpha = 0.5) +
  geom_text_repel(aes(label = id), size = 3, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  scale_color_discrete(name = "Category") +
  scale_y_continuous(breaks = 1:4, labels = c("Patient room", "Corridor", "Office \n Nursing station", "Restroom")) +
  scale_x_continuous(breaks = 1:(sum(rooms_coords$y == 1)+1), labels = c(paste0("Chambre ", 1:(sum(rooms_coords$y == 1)+1)))) +
  labs(title = 'Individual position at time : {begin_date + (frame_time * 30)}, at time step : {frame_time}', x = 'Patient room', y = 'Type of room') +
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 12, hjust = 0.5), # 
    axis.text = element_text(size = 18), 
    axis.title = element_text(size = 18), 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22) 
  )+
  coord_cartesian(clip = "off")

#print(p)
animation <- p +
  transition_time(time) +
  ease_aes('linear')


nframes <- nrow(data%>%distinct(time))
fps <- 50
duration <- nframes / fps
animate(animation, renderer = gifski_renderer("deplacement.gif"), width = 2000, height = 1333, fps = 50, nframes = nframes, duration = duration)


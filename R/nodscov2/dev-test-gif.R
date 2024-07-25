library(ggplot2)
library(gganimate)
library(gifski)
library(dplyr)
library(ggrepel)
library(tidyr)

wd <- getwd()
loc_path <- file.path(wd, "out", "loc-nodscov2","dev-localization-nodscov2.RData")
load(loc_path)

distribute_points <- function(n) {
  theta <- seq(0, 2 * pi, length.out = n + 1)[-1]
  radius <- 0.2  # Radius of the circle where the points will be plotted
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  return(data.frame(offset_x = x, offset_y = y))
}


rooms <- data.frame(
  id_room = c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,  19, 20, 21, 22, 23),
  room = c("Patient Room 1", "Patient Room 2", "Patient Room 3", "Patient Room 4", "Patient Room 5", 
           "Patient Room 6", "Patient Room 7", "Patient Room 8", "Patient Room 9", "Patient Room 10", 
           "Patient Room 11", "Patient Room 12", "Patient Room 13", "Patient Room 14", "Patient Room 15", 
           "Patient Room 16", "Medical Restroom", "Paramedical Restroom", "Nursing Station", "Office", "Corridor")
)



rooms_coords <- rooms %>%
  mutate(
    x = case_when(
      room == "Medical Restroom" ~ 1,
      room == "Paramedical Restroom" ~ 1,
      room == "Office" ~ 2.5,
      room == "Nursing Station" ~ 2.5,
      room == "Patient Room 1" ~ 4,
      room == "Patient Room 2" ~ 5,
      room == "Patient Room 3" ~ 6,
      room == "Patient Room 4" ~ 7,
      room == "Patient Room 5" ~ 4,
      room == "Patient Room 6" ~ 5,
      room == "Patient Room 7" ~ 6,
      room == "Patient Room 8" ~ 7,
      room == "Patient Room 9" ~ 4,
      room == "Patient Room 10" ~ 5,
      room == "Patient Room 11" ~ 6,
      room == "Patient Room 12" ~ 7,
      room == "Patient Room 13" ~ 4,
      room == "Patient Room 14" ~ 5,
      room == "Patient Room 15" ~ 6,
      room == "Patient Room 16" ~ 7,
      room == "Corridor" ~ 1.75
    ),
    y = case_when(
      room == "Medical Restroom" ~ 5.5,
      room == "Office" ~ 5.5,
      room == "Paramedical Restroom" ~ 4,
      room == "Nursing Station" ~ 4,
      room == "Patient Room 1" ~ 5.5,
      room == "Patient Room 2" ~ 5.5,
      room == "Patient Room 3" ~ 5.5,
      room == "Patient Room 4" ~ 5.5,
      room == "Patient Room 5" ~ 4,
      room == "Patient Room 6" ~ 4,
      room == "Patient Room 7" ~ 4,
      room == "Patient Room 8" ~ 4,
      room == "Patient Room 9" ~ 2.5,
      room == "Patient Room 10" ~ 2.5,
      room == "Patient Room 11" ~ 2.5,
      room == "Patient Room 12" ~ 2.5,
      room == "Patient Room 13" ~ 1,
      room == "Patient Room 14" ~ 1,
      room == "Patient Room 15" ~ 1,
      room == "Patient Room 16" ~ 1,
      room == "Corridor" ~ (2.5 + 1)/2
    )
  ) %>%
  mutate(localization = id_room) %>%
  distinct(id_room, .keep_all = TRUE) %>%
  select(localization, x, y, room)




new_begin_date <- begin_date + (5*60*60) ## OFFSET TO START THE DAY AT 12AM
new_end_date <- new_begin_date + (24*60*60)
print(new_end_date - new_begin_date)
t_begin <- 5*60*2 ## OFFSET t
t_end <- (t_begin + 24*60*2) - 1


cat_medical <- c("physician", "ext physician", "reeducation staff", "other")
cat_paramedical <- c("nurse", "student nurse", "aux nurse")

data <- do.call(rbind, global_localization) %>%
  filter(between(time, t_begin, t_end)) %>%
  filter(time < 620) %>%
  filter(localization != -1) %>%
  left_join(rooms_coords, by = c("localization")) %>%
  left_join(admission, by = c("id")) %>%
  mutate(cat = ifelse(is.na(cat), "patient", cat))  %>%
  filter(cat != "patient") %>%
  mutate(type = ifelse(cat %in% cat_medical, 'Medical', 'Paramedical'))




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


##########
## Plot ##
##########
## Color palette
## ALPHA = 1
pal = c('Medical' = "#5CD6D6", 'Paramedical' = "#A9B9E8", 'Patient' = "#FFA766", 'Room' = "#666699")


p <- ggplot() +
  geom_rect(data = rooms_coords, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5), color = "black", fill = NA) +
  geom_text(data = rooms_coords, aes(x = x, y = (y + 0.6), label = room), size = 10) +
  labs(title = 'Individual position at time: {begin_date + (frame_time * 30)},  at time step: {frame_time}') +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 45, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 35),
    legend.title = element_text(size = 40),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off")

animation <- p +
  geom_point(data = data, aes(x = x_adj, y = y_adj, colour = type, group = id), size = 5) +
  geom_text_repel(data = data, aes(x = x_adj, y = y_adj, label = id), size = 5, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
  scale_colour_manual(values = pal, name = "Type") +
  transition_time(time) +
  ease_aes('linear')

nframes <- nrow(data%>%distinct(time))
fps <- 50
duration <- nframes / fps
animate(animation, renderer = gifski_renderer("trajectories-gif.gif"), width = 3000, height = 1333, fps = 50, nframes = nframes, duration = duration)

## .mp4 conversion with ffmpeg
## ffmpeg -i trajectories-gif.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" trajectories.mp4









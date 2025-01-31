################################################################################
##
##    CODE TO CREATE A GIF FROM THE RECONSTRUCTED LOCATION DATA
##
################################################################################
# Clean environment
rm(list = ls())

# Libraries
library(ggplot2)
library(gganimate)
library(gifski)
library(dplyr)
library(ggrepel)
library(tidyr)

# Helper functions
source("R/nodscov2/helper-functions.R")

# Load data 
loc_path <- file.path("data", "data-synthetic-graphs","loc", "nodscov2-reconstructed-locations.rda")
load(loc_path)

# Create room coordinates
rooms <- rooms %>%
  distinct(id_room, room) %>%
  mutate(id_room = as.integer(id_room)) %>%
  mutate(room = ifelse(id_room<=17, paste("Patient Room", room), room)) %>%
  mutate(room = case_when(
    room == "Nursing station" ~ "Nursing Station",
    room == "Medical Restroom" ~ "Medical Staff Room",
    room == "Paramedical Restroom" ~ "Paramedical Staff Room", 
    .default = room
  ))

rooms_coords <- rooms %>%
  mutate(
    x = case_when(
      room == "Medical Staff Room" ~ 1,
      room == "Paramedical Staff Room" ~ 1,
      room == "Office" ~ 2.5,
      room == "Nursing Station" ~ 2.5,
      room == "Patient Room 1" ~ 4,
      room == "Patient Room 2" ~ 5,
      room == "Patient Room 3" ~ 6,
      room == "Patient Room 4" ~ 7,
      room == "Patient Room 5" ~ 4,
      room == "Patient Room 6" ~ 5,
      room == "Patient Room 8" ~ 6,
      room == "Patient Room 9" ~ 7,
      room == "Patient Room 10" ~ 4,
      room == "Patient Room 11" ~ 5,
      room == "Patient Room 12" ~ 6,
      room == "Patient Room 13" ~ 7,
      room == "Patient Room 14" ~ 4,
      room == "Patient Room 15" ~ 5,
      room == "Patient Room 16" ~ 6,
      room == "Patient Room 17" ~ 7,
      room == "Corridor" ~ 1.75
    ),
    y = case_when(
      room == "Medical Staff Room" ~ 5.5,
      room == "Office" ~ 5.5,
      room == "Paramedical Staff Room" ~ 4,
      room == "Nursing Station" ~ 4,
      room == "Patient Room 1" ~ 5.5,
      room == "Patient Room 2" ~ 5.5,
      room == "Patient Room 3" ~ 5.5,
      room == "Patient Room 4" ~ 5.5,
      room == "Patient Room 5" ~ 4,
      room == "Patient Room 6" ~ 4,
      room == "Patient Room 8" ~ 4,
      room == "Patient Room 9" ~ 4,
      room == "Patient Room 10" ~ 2.5,
      room == "Patient Room 11" ~ 2.5,
      room == "Patient Room 12" ~ 2.5,
      room == "Patient Room 13" ~ 2.5,
      room == "Patient Room 14" ~ 1,
      room == "Patient Room 15" ~ 1,
      room == "Patient Room 16" ~ 1,
      room == "Patient Room 17" ~ 1,
      room == "Corridor" ~ (2.5 + 1)/2
    )
  ) %>%
  mutate(location = id_room) %>%
  distinct(id_room, .keep_all = TRUE) %>%
  select(location, x, y, room)

#
# new_begin_date <- begin_date + (5*60*60) ## OFFSET TO START THE DAY AT 12AM
# new_end_date <- new_begin_date + (24*60*60)
# print(new_end_date - new_begin_date)
t_begin <- 1 #5*60*2 ## OFFSET t
t_end <- (t_begin + 48*60*2) - 1 #as.numeric(difftime(end_date, begin_date, units = "secs"))/30

# Get location data
data <- do.call(rbind, global_location) %>%
  filter(between(time, t_begin, t_end)) %>%
  # filter(time < 620) %>%
  filter(location != -1) %>%
  left_join(admission, by = c("id")) %>%
  filter(cat != "Patient") %>%
  mutate(location = as.integer(location)) %>%
  left_join(rooms_coords, by = c("location")) %>%
  rename(type=cat)

max_n <- data %>%
  group_by(location, time) %>%
  summarise(n = n(), .groups = 'drop') %>%
  summarise(max_n = max(n)) %>%
  pull(max_n)

data <- data %>%
  group_by(location, time) %>%
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

# Plots-------------------------------------------------------------------------
## ALPHA = 1
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
animate(animation, renderer = gifski_renderer("fig/loc-reconstruction/trajectories-gif.gif"), width = 3000, height = 1333, fps = 50, nframes = nframes, duration = duration)

## .mp4 conversion with ffmpeg
## ffmpeg -i trajectories-gif.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" trajectories.mp4

############ OLD #########

## OLD VERSION FOR ROOMS COORDINATES
# rooms_coords <- rooms %>%
#   mutate(
#     x = case_when(
#       grepl("001", id) ~ id_room,
#       grepl("Corridor", room) ~ 1,
#       room == "Office" ~ 1,
#       room == "Nursing station" ~ 3,
#       room == "Medical Restroom"~ 1,
#       room == "Paramedical Restroom"~ 3,
#       TRUE ~ NA_real_
#     ),
#     y = case_when(
#       grepl("001", id) ~ 1,
#       grepl("Corridor", room) ~ 2,
#       room == "Office" ~ 3,
#       room == "Nursing station" ~ 3,
#       grepl("Restroom", room) ~ 4,
#       TRUE ~ NA_real_
#     )
#   ) %>%
#   mutate(location = id_room) %>%
#   distinct(id_room, .keep_all = TRUE) %>%
#   select(location, x, y)

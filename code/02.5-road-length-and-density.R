# This script calculates the total road length and road density (km/km²) for each
# study area, including a 3km buffer zone. The analysis uses the previously processed
# OpenStreetMap road data to quantify the extent of road infrastructure in each region.

library(tidyverse)
library(sf)
library(units)


# Load data
place_roads <- st_read("data/place_roads_2025-04-22-2.gpkg")
place_buff <- st_read("data/place_buffer3km.gpkg")


# Calculate road length and density ------------------------------------------

# Calculate road length using map
road_length <- place_roads |>
  group_split(name) |>
  map_dfr(~ tibble(
    name = first(.x$name),
    length = set_units(sum(st_length(.x$geom)), "km")
  ))

# Join the road length with the place buffer data
road_density <- place_buff |>
  mutate(
    area = set_units(st_area(geom), "km^2") # ,
    # area = as.numeric(area)
  ) |>
  left_join(road_length, by = "name") |>
  mutate(
    length = length, # as.numeric(length),
    road_density = length / area # This will give km of road per km²
  ) |>
  st_drop_geometry() |>
  tibble()

road_density |> View()

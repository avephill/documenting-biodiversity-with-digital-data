# This script performs Nearest Neighbor Index (NNI) analysis on verified species
# occurrences to assess spatial clustering patterns, comparing distributions between
# iNaturalist observations and herbarium specimens. The script also generates density
# plots to visualize spatial patterns across the study regions.

library(tidyverse)
library(duckdb)
library(spatialEco)
library(duckdbfs)
library(patchwork)
library(cowplot)
library(ggspatial)
library(sf)
library(spatstat)
library(stars)
library(tictoc)
library(multidplyr)


# NNI Analysis ---------------------------------------


verified_occ <- st_read("data/verified_occurrences_dedup_2025-04-21.gpkg")
# verified_occ |> head() |> View()
# verified_occ |>
#   filter(place_name == "One Tam", decimallatitude > 38.05) |>
#   slice_head(n = 100) |>
#   pull(decimallatitude)

cluster <- new_cluster(2)

tic()
verified_nni_spec <- verified_occ |>
  # filter(place_name == "One Tam") |>
  group_split(species, place_name, basisofrecord) %>%
  map_dfr(~ {
    # browser()
    if (nrow(.x) > 1) {
      nni_results <- nni(.x, win = "extent")
      .x |>
        distinct(species, place_name, basisofrecord) |>
        cbind(as_tibble(nni_results))
    }
  })
toc()

# With family?
tic()
verified_nni_family <- verified_occ |>
  # filter(place_name == "One Tam") |>
  group_split(family, place_name, basisofrecord) %>%
  map_dfr(~ {
    # browser()
    if (nrow(.x) > 1) {
      nni_results <- nni(.x, win = "extent")
      .x |>
        distinct(family, place_name, basisofrecord) |>
        cbind(as_tibble(nni_results))
    }
  })
toc()

complete_nni <- verified_nni_family |>
  filter(is.finite(NNI)) |>
  as_tibble() |>
  filter(p < .05)

nni_ttest <- complete_nni |>
  mutate(NNI_log = log(NNI)) |>
  filter(is.finite(NNI_log)) |>
  group_split(place_name) |>
  map_dfr(~ {
    nni_inat <- .x |>
      filter(basisofrecord == "HUMAN_OBSERVATION") |>
      pull(NNI_log)
    nni_herb <- .x |>
      filter(basisofrecord == "PRESERVED_SPECIMEN") |>
      pull(NNI_log)

    # browser()

    t_result <- t.test(nni_inat, nni_herb)

    .x |>
      distinct(place_name) |>
      cbind(as_tibble(t_result[["p.value"]])) |>
      rename(p_value = value)
  })

nni_figuredf <- complete_nni |>
  left_join(nni_ttest)

nni_figuredf |>
  write_csv(sprintf(
    "results/nni_results_%s.csv",
    format(Sys.Date(), "%Y-%m-%d")
  ))

old_test <- read_csv("results/nni_results.csv")

# t.test(log(NNI$NNI.h[which(NNI$p.h < 0.05 & NNI$p.o < 0.05)]), log(NNI$NNI.o[which(NNI$p.h < 0.05 & NNI$p.o < 0.05)])) # test for a difference bewteen herbarium and inat







# PDF Density Analysis ---------------------------------------
# This is scrap for now, unless we go back to it

densityMaker <- function(input_occ, place_name) {
  place_owin <- places |>
    st_transform(3310) |>
    filter(name == place_name) |>
    st_geometry() |>
    as.owin()

  as.ppp(
    input_occ |>
      st_transform(3310) |>
      st_filter(st_as_sf(place_owin) |> st_set_crs(3310)),
    W = place_owin
  ) |>
    density.ppp(dimyx = 200) |>
    st_as_stars() |>
    st_set_crs(3310) |>
    st_as_sf() |>
    mutate(geom_text = st_as_text(geometry))
}

fullDensCompiler <- function(full_occ, place_name) {
  # browser()
  inat_tamdens <- densityMaker(full_occ |> filter(basisofrecord == "HUMAN_OBSERVATION"), place_name)

  herb_tamdens <- densityMaker(full_occ |> filter(basisofrecord == "PRESERVED_SPECIMEN"), place_name)

  all_tamdens <- densityMaker(full_occ, place_name)

  joinem <- all_tamdens |>
    rename(sum_vals = v) |>
    left_join(inat_tamdens |> as_tibble() |>
      select(inat_vals = v, geom_text), by = "geom_text") |>
    left_join(herb_tamdens |> as_tibble() |>
      select(herb_vals = v, geom_text), by = "geom_text") |>
    mutate(place_name = place_name)
}

full_dens <- map(places |> pull(name), function(x) fullDensCompiler(full_occ, x)) |> bind_rows()


full_dens_long <- full_dens |>
  mutate(diff_vals = herb_vals - inat_vals) |>
  pivot_longer(
    cols = c(sum_vals, inat_vals, herb_vals, diff_vals),
    names_to = "variable",
    values_to = "value"
  ) |>
  mutate(
    sqrt_value = sqrt(value),
    value_km2 = value * 1e6
  )



trails <- st_read("data/place_roads.gpkg") |>
  filter(!st_geometry_type(geom) %in% c("POLYGON", "MULTIPOLYGON")) |>
  st_filter(places)
# trails_grouped <- trails |>
#   group_by(name) |>
#   summarize(geom = st_union(geom))

full_dens_plotready <- full_dens_long |>
  mutate(variable_label = case_when(
    variable == "diff_vals" ~ "Herbarium - iNat",
    variable == "herb_vals" ~ "Herbarium",
    variable == "inat_vals" ~ "iNat",
    variable == "sum_vals" ~ "Herbarium + iNat"
  ))

facet_plots <- full_dens_plotready |>
  filter(variable != "diff_vals") |>
  group_split(place_name, variable_label) |>
  map(~ {
    tam_most_plot <- ggplot() +
      geom_sf(
        data = .x,
        aes(fill = value_km2), color = NA
      ) +
      geom_sf(
        data = trails |>
          filter(name == unique(.x$place_name)),
        color = "white",
        linewidth = .05 # , size = .01
      ) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        legend.position = "right",
        legend.key.size = unit(.2, "cm"), # Increase size of legend keys
        legend.text = element_text(size = 5), # Increase size of legend text
        legend.title = element_text(size = 7),
        plot.margin = margin(10, 10, 10, 10) # Add margin around each panel
      ) +
      coord_sf() +
      scale_fill_viridis_c(option = "magma") +
      # facet_grid(cols = vars(variable_label)) +
      labs(
        fill = expression("Obs/km"^2),
        title = .x$variable_label
      )
  })


# x <- plot_grid(
#   plotlist = facet_plots,
#   # rel_widths = c(2, 1, 1),
#   nrow = 3, ncol = 3

#   # labels = "auto"
# )
# ggsave("results/testing.png", x, width = 12, height = 12, units = "in")

# combined_plot <- wrap_plots(
#   facet_plots,
#   nrow = 3
# )
# combined_plot


# ggsave("results/density_plots1.png", combined_plot, width = 7, height = 4)

diff_plots <-
  full_dens_plotready |>
  filter(variable == "diff_vals") |>
  group_split(place_name) |>
  map(~ {
    tam_diff_plot <- ggplot() +
      geom_sf(
        data = .x,
        aes(fill = value), color = NA
      ) +
      geom_sf(
        data = trails |>
          filter(name == unique(.x$place_name)),
        color = "white",
        linewidth = .05 # , size = .01
      ) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        legend.position = "right",
        legend.key.size = unit(.2, "cm"), # Increase size of legend keys
        legend.text = element_text(size = 5), # Increase size of legend text
        legend.title = element_text(size = 7),
        plot.margin = margin(10, 10, 10, 10) # Add margin around each panel
      ) +
      # scale_fill_gradient2(
      #   low = "blue", # Color for low values
      #   mid = "white", # Color for midpoint
      #   high = "red", # Color for high values
      #   midpoint = 0, # Center the scale at 0
      #   guide = guide_colorbar(title = "Difference (Herb - iNat)") # Add legend title
      # ) +
      scale_fill_gradient2(
        low = "#D73027", # Red
        mid = "#4575B4", # Neutral blue
        high = "#1A9850", # Green
        midpoint = 0,
        guide = guide_colorbar(title = "Difference (Herb - iNat)") # Add legend title
      ) +
      labs(
        fill = expression("Obs/km"^2),
        title = .x$variable_label
      ) +
      coord_sf()
  })


# Define target positions for diff_plots
target_positions <- c(4, 8, 12)
total_length <- 12
# Create big_plotlist using map()
big_plotlist <- map(seq_len(total_length), function(i) {
  if (i %in% target_positions) {
    # Insert from diff_plots
    diff_plots[[which(target_positions == i)]]
  } else {
    # Insert from facet_plots
    facet_plots[[i - sum(i > target_positions)]]
  }
})
x <- plot_grid(
  plotlist = big_plotlist,
  # rel_widths = c(2, 1, 1),
  nrow = 3, ncol = 4

  # labels = "auto"
)



ggsave("results/density_plots.png", x, width = 12, height = 10)

library(tidyverse)
library(sf)

rovbase_data <- readxl::read_excel("data_raw/MSKO13062022110710932.xlsx", guess_max = Inf) %>% 
  filter(`Art (Prøve)` == "Björn", 
         lubridate::year(Funnetdato) %in% c(2010, 2016, 2021), lubridate::month(Funnetdato) %in% 8:10, 
         Prøvetype == "Spillning",
         Fylke == "Norrbottens län (S)")

rovbase_data %>% mutate(År = lubridate::year(Funnetdato)) %>% 
  count(År, Prøvestatus) %>% 
  mutate(Prøvestatus = ifelse(is.na(Prøvestatus), "Status saknas", Prøvestatus)) %>% 
  pivot_wider(names_from = Prøvestatus, values_from = n, values_fill = 0) %>% 
  mutate(`Antal prover` = `Inget resultat` + `Resultat erhållet` + `Status saknas`,
         `Resultat erhållet` = paste0(`Resultat erhållet`, " (", round(100 *  `Resultat erhållet` / `Antal prover`), "%)")) %>% 
  select(År, `Antal prover`, `Resultat erhållet`, `Inget resultat`, `Status saknas`) %>% 
  write_csv("summary_table.csv")

cities <- tibble(addr = c("Kiruna", "Gällivare", "Jokkmokk", "Boden, Sweden", "Luleå", "Piteå"),
                 city = c("Kiruna", "Gällivare", "Jokkmokk", "Boden", "Luleå", "Piteå")) %>% 
  tidygeocoder::geocode(addr) %>% 
  st_as_sf(crs = 4326, coords = c("long", "lat")) %>% 
  st_transform(crs = 3006)
save(cities, file = "cities.Rdata")

data <- rovbase_data %>% 
  mutate(id = str_sub(Individ, 1, 8),
         year = lubridate::year(Funnetdato),
         week = lubridate::isoweek(Funnetdato)) %>%
  select(id, year, week, sex = Kjønn, date = Funnetdato, lat = `Nord (UTM33/SWEREF99 TM)`, lon = `Øst (UTM33/SWEREF99 TM)`) %>% 
  filter(!is.na(id), sex %in% c("Hona", "Hane")) %>% 
  group_by(year, id, sex) %>% 
  mutate(mlon = mean(lon), mlat = mean(lat))

write_csv(data, "data_BD.csv")

map2 <- st_read("LanSweref99TM/Lan_Sweref99TM_region.shp") %>% 
  filter(LnKod == 25)
bbox <- st_bbox(map2)
bbox["xmin"] <- bbox["xmin"] - 10000
bbox["ymin"] <- bbox["ymin"] - 10000
bbox["xmax"] <- bbox["xmax"] + 10000
bbox["ymax"] <- bbox["ymax"] + 10000


grid <- st_make_grid(bbox, cellsize = c(10000, 10000)) %>% 
  st_as_sf() %>% 
  mutate(in_region = st_intersects(x, st_union(map2), sparse = FALSE) %>% as.logical(),
         distance = st_distance(x, st_union(map2))) %>% 
  filter(as.numeric(distance) < 10000) %>% 
  rename(grid = x) %>% 
  mutate(grid_id = 1:n()) %>% 
  rowwise() %>% 
  mutate(center = st_centroid(grid),
         lat = st_coordinates(center)[2], 
         lon = st_coordinates(center)[1])
grid_counts <- 
  home_region %>% st_as_sf(coords = c("mlon", "mlat"), crs = st_crs(grid)) %>% 
  st_join(grid) %>% 
  as_tibble() %>% 
  count(year, grid_id, sex) %>% 
  right_join(expand_grid(year = unique(home_region$year), 
                         grid_id = 1:max(grid$grid_id),
                         sex = c("Hona", "Hane")), by = c("year", "grid_id", "sex")) %>% 
  mutate(n = replace_na(n, 0)) %>% 
  arrange(year, grid_id) %>% 
  left_join(grid %>% as_tibble() %>% select(grid_id, lon, lat, in_region), by = "grid_id")
grid_counts_total <- bind_rows(grid_counts,
                               grid_counts %>% group_by(year, grid_id, lon, lat) %>% 
                                 summarise(n = sum(n), .groups = "drop") %>% 
                                 mutate(sex = "Båda"))

write_csv(grid_counts_total, "grid_counts.csv")

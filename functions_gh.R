plot_spatial <- function(fit, plot_var){
  library(sf)
  map <- st_read("LanSweref99TM/Lan_Sweref99TM_region.shp", quiet = TRUE) %>% 
    filter(LnKod == 25)
  map_k <- st_read("KommunSweref99TM/Kommun_Sweref99TM_region.shp", quiet = TRUE) %>% 
    mutate(Lan = str_sub(KnKod, 1, 2)) %>% 
    filter(Lan == "25")
  bbox <- st_make_grid(st_as_sf(fit$predictions, crs = st_crs(map), coords = c("lon", "lat")), n = c(1, 1))
  cover <- st_difference(bbox, map)
  fit[["predictions"]] %>% ggplot() + 
    geom_contour_filled(aes(x = lon, y = lat, z = {{plot_var}})) +
    geom_point(data = fit[["data"]], aes(x = lon, y = lat), color = "white", size = .05, ) +
    geom_sf(data = map, fill = NA, color = "white", size = .1, stroke = 0) +
    theme_void() + theme(legend.position = "bottom") + labs(fill = "") +
    geom_sf(data = map_k, size = .1, color = "white", fill = NA) +
geom_sf(data = cover, fill = "white") 
}

fit_spatial_2021 <- function(data, grid_counts){
  suppressPackageStartupMessages(library(mgcv))

  intensity_fit <- gam(n ~ s(lon, lat, k = 50), 
                       data = grid_counts %>% filter(in_region == TRUE), 
                       family = "poisson")
  agg_data <- data %>% 
    select(week, id, lat = mlat, lon = mlon) %>% 
    distinct() %>% 
    count(id, lon, lat)
  
  effort_fit <- gam(n ~ s(lon, lat, k = 50),
                    data = agg_data,
                    family = countreg::ztpoisson)
  predictions <- grid_counts %>% 
    mutate(mu = exp(predict(effort_fit, newdata = grid_counts)),
           p = 1 - exp(-mu),
           lambda = exp(predict(intensity_fit, newdata = grid_counts)) / p
    )
  list(predictions = predictions, data = agg_data, intensity_fit = intensity_fit, effort_fit = effort_fit)
}

pretty_ci <- function(l, u, loglin = FALSE){
  if (loglin == FALSE){
    paste0("(", round(l), ", ", round(u), ")")
  }
  else
  {
    paste0("(", signif((exp(l)-1)*100, 3), ", ", signif((exp(u)-1)*100, 3), ")")
  }
}

process_data <- function(data, model){
  data %>% group_by(id, sex, week) %>% 
    summarise(n = 1, .groups = "drop") %>% 
    ungroup() %>% 
    arrange(week) %>%
    pivot_wider(id_cols = c("id", "sex"), 
                names_from = week, values_from = n, 
                values_fill = c(n = 0)) %>% 
    unite("ch", -c("id", "sex"), sep = "") %>% 
    RMark::process.data(model = model, groups = "sex")
}

fit_models <- function(data){
  p.time.mixture <- list(formula = ~time + mixture, share = TRUE)
  p.time <- list(formula = ~time, share = TRUE)
  p.time.sex <- list(formula = ~time + sex, share = TRUE)
  pi.1 <- list(formula = ~1)
  pi.sex <- list(formula = ~sex)
  f0.sex <-  list(formula = ~sex)
  
  fit1 <- process_data(data, "FullHet") %>% 
    RMark::mark(model = "HetClosed", model.parameters = list(pi = pi.1, p = p.time.mixture, f0 = f0.sex), 
                output = FALSE, silent = TRUE, delete = TRUE)
  fit2 <- process_data(data, "FullHet") %>% 
    RMark::mark(model = "HetClosed", model.parameters = list(pi = pi.sex, p = p.time.mixture, f0 = f0.sex), 
                output = FALSE, silent = TRUE, delete = TRUE)
  fit3 <- process_data(data, "Closed") %>% 
    RMark::mark(model = "Closed", model.parameters = list(p = p.time.sex, f0 = f0.sex), 
                output = FALSE, silent = TRUE, delete = TRUE)
  fit4 <- process_data(data, "Closed") %>% 
    RMark::mark(model = "Closed", model.parameters = list(p = p.time, f0 = f0.sex), 
                output = FALSE, silent = TRUE, delete = TRUE)
  
  all_fit <- tibble(fit = list(fit1, fit2, fit3, fit4))
    
  table <- all_fit %>% rowwise() %>% 
    mutate(nm = fit[["results"]][["derived"]][["N Population Size"]][["estimate"]][1] %>% round(),
           nf = fit[["results"]][["derived"]][["N Population Size"]][["estimate"]][2] %>% round(),
           n2f = 2 * nf,
           nfm = nm + nf,
           nm_l = fit[["results"]][["derived"]][["N Population Size"]][["lcl"]][1],
           nm_u = fit[["results"]][["derived"]][["N Population Size"]][["ucl"]][1],
           nf_l = fit[["results"]][["derived"]][["N Population Size"]][["lcl"]][2],
           nf_u = fit[["results"]][["derived"]][["N Population Size"]][["ucl"]][2],
           n2f_l = 2 * nf_l,
           n2f_u = 2 * nf_u,
           tot_var = fit[["results"]][["derived.vcv"]][["N Population Size"]] %>% sum(),
           nfm_l = nfm / exp(1.96 * sqrt(log(1 + tot_var / nfm^2))),
           nfm_u = nfm * exp(1.96 * sqrt(log(1 + tot_var / nfm^2))),
           AICc = fit[["results"]][["AICc"]],
           model = fit[["model.name"]] %>% str_remove_all("~|c\\(\\)")
           ) %>% 
    ungroup() %>% 
    mutate(dAICc = round(AICc - min(AICc), 1))
  table
}



rootogram_ztp <- function(object){
  count(tibble(captures = object$y), captures) %>% 
    right_join(tibble(captures = 1:max(object$y)), by = "captures") %>% 
    mutate(n = ifelse(is.na(n), 0, n)) %>% 
    rowwise() %>% 
    mutate(freq = sum(countreg::dztpois(captures, lambda = exp(predict(object))))) %>% 
    ggplot(aes(x = captures, y = sqrt(freq))) + 
    geom_rect(aes(xmin = captures-.4, xmax = captures+.4, ymax = sqrt(freq), ymin = sqrt(freq) - sqrt(n)), fill = "grey") +
    geom_point() + geom_line() +
    geom_hline(yintercept = 0) +
    theme_bw() + scale_x_continuous(breaks = 1:max(object$y)) + labs(x = "")
}




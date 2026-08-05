suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(tidyverse))

# uploading landscape and climate data------------------------
tif_metad <- c("other.tsv", "worldclim.tsv" ) %>% 
  map_dfr(~readr::read_delim(paste0("rasters/", .x), "\t", show_col_types = FALSE))
rst <- tif_metad$file %>% 
  paste0("./rasters/", .) %>% 
  terra::rast()

tif_metad <- tif_metad %>% 
  mutate (name_rus = c ("Высота, м н.у.м.",
                        "Степень антропогенного воздейтсвия, гга",
                        "Среднегодовая температура, °C",
                        "Среднесуточная амплитуда температур, °C",
                        "Изотермальность, %",
                        "Температурная сезонность, %",
                        "Максимальная температура наиболее теплого месяца,°C",
                        "Минимальная температура наиболее холодного месяца, °C",
                        "Годовая амплитуда температур, °C",
                        "Средняя температура наиболее влажного квартала, °C",
                        "Средняя температура наиболее сухого квартала, °C",
                        "Средняя температура наиболее теплого квартала, °C",
                        "Средняя температура наиболее холодного квартала, °C",
                        "Годовое количество осадков, мм",
                        "Количество осадков в наиболее влажный месяц, мм",
                        "Количество осадков в наиболее сухой месяц, мм",
                        "Сезонность осадков, %",
                        "Количество осадков в наиболее влажный квартал, мм",
                        "Количество осадков в наиболее сухой квартал, мм",
                        "Количество осадков в наиболее теплый квартал, мм",
                        "Количество осадков в наиболее холодный квартал, мм"))

# data family------------------------------
tbl_fam <- readRDS("res.RDS") %>%
  filter(
    family %in% c('Linyphiidae', 'Lycosidae', 'Gnaphosidae', 'Thomisidae', 'Theridiidae', 'Araneidae', 'Salticidae',  'Tetragnathidae', 'Philodromidae'),
    occurrencestatus == 'present',
    !is.na(decimallongitude), 
    !is.na(decimallatitude)
  ) %>% 
  select(family, decimallongitude, decimallatitude) %>% 
  # distinct() %>%
  st_as_sf(crs=4326, coords = c("decimallongitude", "decimallatitude"))

# data species------------------------------
tbl_spc <- readRDS("res.RDS") %>%
  filter(
    str_detect(taxonrank, "SPEC"),
    occurrencestatus == 'present',
    !is.na(decimallongitude), 
    !is.na(decimallatitude)
  ) %>% 
  select(genus, specificepithet, decimallongitude, decimallatitude) %>% 
  unite("spec", genus, specificepithet, sep = " ") %>% 
  # distinct() %>%
  st_as_sf(crs=4326, coords = c("decimallongitude", "decimallatitude"))

top100 <- tbl_spc %>% 
  st_drop_geometry() %>% 
  count(spec, sort = T) %>% 
  filter(n >= 100) %>% 
  pull(spec) 

tbl_spc <- filter(tbl_spc, spec %in% top100) # mass species

# result tables-----------------------------------------------

# family
res_fam <- tbl_fam %>% 
  terra::vect() %>% 
  terra::extract(rst, .) %>% 
  select(-ID)

res_fam <- cbind(
  st_drop_geometry(tbl_fam),
  res_fam, 
  mutate_all(res_fam, ~.x^2) %>% `colnames<-`(paste0(colnames(.), "_2"))
) %>% 
  as_tibble() 

# species
res_sp <- tbl_spc %>% 
  terra::vect() %>% 
  terra::extract(rst, .) %>% 
  select(-ID)

res_sp <- cbind(
  st_drop_geometry(tbl_spc),
  res_sp, 
  mutate_all(res_sp, ~.x^2) %>% `colnames<-`(paste0(colnames(.), "_2"))
) %>% 
  as_tibble() 

# models families---------------------------
models2 <- expand_grid(
  fct =  res_fam %>% 
    select(-1:-2) %>% 
    select(- ends_with("_2")) %>%
    colnames(),
  fam = sort(unique(tbl_fam$family))
) %>% 
  filter(str_detect(fct, "landc", negate = TRUE))

tmp <- lapply(
  # 1:nrow(models2[1:3,]), function(x){
  1:nrow(models2), function(x){
    # input <- 
    fct   <- models2$fct[x]
    
    input <- res_fam %>% 
      filter(family == models2$fam[x]) %>% 
      select(pred = all_of(fct)) %>% 
      #filter_out(is.na(pred)) %>% 
      filter(!is.na(pred)) %>% 
      mutate(
        interval = cut(pred, breaks = 15, include.lowest = TRUE), 
        interval = str_remove_all(interval, "\\(|\\]|\\)|\\[")
      )
    
    res <- input %>% 
      mutate(mid_pred = map_dbl(input$interval, ~mean(as.numeric(str_split(.x, ",")[[1]])))) %>% 
      group_by(mid_pred) %>% 
      summarise(n_rec = n(), .groups = "drop") %>% 
      mutate(n_rec2 = n_rec/sum(n_rec)*100, mid_pred_2 = mid_pred^2, .after = mid_pred) %>% 
      lm(formula = as.formula("n_rec2 ~ mid_pred + mid_pred_2"), data = .)
    sm_res <- summary(res)
    if(nrow(sm_res$coefficients) != 3) {
      return(NULL)
    } else {
      return(list(
        models = res, 
        summaries = sm_res,
        input_data = res$model[,1:2],
        fct_min = min(input$pred), 
        fct_max = max(input$pred), 
        # fct_steps = sort(unique(map_dbl(input$interval, ~mean(as.numeric(str_split(.x, ",")[[1]]))))), 
        r2 = sm_res$r.squared,
        A       = sm_res$coefficients[3,1],
        A_pval  = sm_res$coefficients[3,4],
        B       = sm_res$coefficients[2,1],
        B_pval  = sm_res$coefficients[2,4],
        C       = sm_res$coefficients[1,1],
        shapiro_test = shapiro.test(res$residuals)$p.value
      ))
    }
  }
) %>% 
  transpose()

models2 <- tmp %>% 
  as_tibble() %>% 
  mutate_at(4:12, flatten_dbl) %>% 
  tibble(models2, .)


# change factors names
models2 <- models2 %>%
  mutate( x0 = -B/(2*A), # x-coordinate of the parabola's vertex
    fct = case_when(
    substr(fct, 1, 12) %in% substr(tif_metad$file, 1, 12) ~ tif_metad$name_rus[match(substr(fct, 1, 12), substr(tif_metad$file, 1, 12))],
    TRUE ~ fct  
  ))

# selected models:
models_selected <- models2 %>% 
  filter(A_pval <= 0.05, A < 0, shapiro_test > 0.05, fct != "Изотермальность, %") %>% #  B_pval <= 0.05
  arrange(desc(r2)) %>% 
  slice(1:12)

mod_fam <- models_selected %>% 
  select(-models, -summaries, -C, -input_data) %>% 
  mutate_at(c("A_pval", "B_pval", "shapiro_test", "A"), ~round(.x, 4)) %>% 
  mutate_at(c("r2"), ~round(.x, 2)) %>% 
  mutate_at(c("B", "fct_min", "fct_max"), ~round(.x, 1)) 

# visualization family -------------------------------------
df_tmp <- models2 %>% 
  filter(A_pval <= 0.05, A < 0, shapiro_test > 0.05, fct != "Изотермальность, %") %>% #  B_pval <= 0.05
  arrange(desc(r2)) %>% 
  slice(1:12) %>% 
  unite(model, fam, fct, sep = " ~ ")

df_observed <- df_tmp %>% 
  split(.$model) %>% 
  map_dfr(~.x$input_data %>% `[[`(1) %>% mutate(model = .x$model)) %>% 
  separate(model, "Family", sep = " ~ ", extra = "drop", remove = F)

df_predicted <- df_tmp %>% 
  # slice(1:2) %>% 
  split(.$model) %>% 
  lapply(function(xx){
    d <- tibble(
      mid_pred = seq(xx$fct_min, xx$fct_max, length.out = 35),
      mid_pred_2 = mid_pred^2
    ) 
    mutate(d, n_rec2 = predict(xx$models[[1]], newdata = d))
  }) %>% 
  map_dfr(rbind, .id = "model") %>% 
  filter(n_rec2 >= 0) %>% 
  separate(model, "Family", sep = " ~ ", extra = "drop", remove = F)

graph_fam <- ggplot(df_observed, aes(mid_pred, n_rec2, color = Family)) + 
  geom_point() + 
  geom_line(data = df_predicted) + 
  facet_wrap(~model,
             ncol = 3,
             labeller = labeller(model = function(x) {
               before <- str_extract(x, "^[^~]+")
               after  <- str_extract(x, "(?<=~).*")
               # Используем str_wrap для разбивки текста на строки
               wrapped_text <- str_wrap(after, width = 35)  # Ширина в символах
               paste(wrapped_text, paste0("Семейство: ", before), sep = "\n")
             }),
             scales = "free") +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  labs(x = NULL, y = "Доля находок, %", 
       color = "Семейства") + 
  theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "#dafdce"),
        strip.text = element_text(size = 8, lineheight = 0.9, hjust = 0.5))  # Уменьшение междустрочного интервала

print(graph_fam)

ggsave("plotfamfinal.pdf", plot = graph_fam, device = cairo_pdf, width = 9, height = 13.5)

# models species---------------------------
models3 <- expand_grid(
  fct =  res_sp %>% 
    select(-1:-3) %>% 
    select(- ends_with("_2")) %>%
    colnames(),
  species = sort(unique(tbl_spc$spec))
) %>% 
  mutate(frm =   paste0("n_rec2 ~ ", fct, " + ", fct, "_2"))


tmp_sp <- lapply(
  # 1:nrow(models3[1:300,]), function(x){
  1:nrow(models3), function(x){
    # input <- 
    fct   <- models3$fct[x]
    
    input <- res_sp %>% 
      filter(spec == models3$species[x]) %>% 
      select(pred = all_of(fct)) %>% 
      #filter_out(is.na(pred)) %>% 
      filter(!is.na(pred)) %>% 
      mutate(
        interval = cut(pred, breaks = 15, include.lowest = TRUE), 
        interval = str_remove_all(interval, "\\(|\\]|\\)|\\[")
      )
    
    res <- input %>% 
      mutate(mid_pred = map_dbl(input$interval, ~mean(as.numeric(str_split(.x, ",")[[1]])))) %>% 
      group_by(mid_pred) %>% 
      summarise(n_rec = n(), .groups = "drop") %>% 
      mutate(n_rec2 = n_rec/sum(n_rec)*100, mid_pred_2 = mid_pred^2, .after = mid_pred) %>% 
      lm(formula = as.formula("n_rec2 ~ mid_pred + mid_pred_2"), data = .)
    sm_res <- summary(res)
    if(nrow(sm_res$coefficients) != 3) {
      return(NULL)
    } else {
      return(list(
        models = res, 
        summaries = sm_res,
        input_data = res$model[,1:2],
        fct_min = min(input$pred), 
        fct_max = max(input$pred), 
        # fct_steps = sort(unique(map_dbl(input$interval, ~mean(as.numeric(str_split(.x, ",")[[1]]))))), 
        r2 = sm_res$r.squared,
        A       = sm_res$coefficients[3,1],
        A_pval  = sm_res$coefficients[3,4],
        B       = sm_res$coefficients[2,1],
        B_pval  = sm_res$coefficients[2,4],
        C       = sm_res$coefficients[1,1],
        shapiro_test = shapiro.test(res$residuals)$p.value
      ))
    }
  }
) %>% 
  transpose()

models3 <- tmp_sp %>% 
  as_tibble() %>% 
  mutate_at(4:12, flatten_dbl) %>% 
  tibble(models3, .)

# change factors names
models3 <- models3 %>%
  mutate(x0 = -B/(2*A), # x-coordinate of the parabola's vertex
    fct = case_when(
    substr(fct, 1, 12) %in% substr(tif_metad$file, 1, 12) ~ tif_metad$name_rus[match(substr(fct, 1, 12), substr(tif_metad$file, 1, 12))],
    TRUE ~ fct  
  ))

# selected models:
models_selected <- models3 %>% 
  filter(A_pval <= 0.05, A < 0, shapiro_test > 0.05, fct != "Изотермальность, %") %>% #  B_pval <= 0.05
  arrange(desc(r2)) %>%
  slice(1:12)

mod_sp <- models_selected %>% 
  select(-models, -summaries, -C, -input_data, -frm) %>% 
  mutate_at(c("A"), ~round(.x, 4)) %>% 
  mutate_at(c("A_pval", "B_pval", "shapiro_test"), ~round(.x, 3)) %>% 
  mutate_at(c("r2"), ~round(.x, 2)) %>% 
  mutate_at(c("B", "fct_min", "fct_max"), ~round(.x, 1)) 

# write.table("export_fam/tsv", sep = "\t")
rm(tmp, tmp_sp)

# visualization species ---------------------------------
df_tmp_sp <- models3 %>% 
  filter(A_pval <= 0.05, A < 0, shapiro_test > 0.05, fct != "Изотермальность, %") %>% #  B_pval <= 0.05
  arrange(desc(r2)) %>% 
  slice(1:12) %>% 
  unite(model, species, fct, sep = " ~ ")

df_observed_sp <- df_tmp_sp %>% 
  split(.$model) %>% 
  map_dfr(~.x$input_data %>% `[[`(1) %>% mutate(model = .x$model)) %>% 
  separate(model, "Species", sep = " ~ ", extra = "drop", remove = F)

df_predicted_sp <- df_tmp_sp %>% 
  # slice(1:2) %>% 
  split(.$model) %>% 
  lapply(function(xx){
    d <- tibble(
      mid_pred = seq(xx$fct_min, xx$fct_max, length.out = 35),
      mid_pred_2 = mid_pred^2
    ) 
    mutate(d, n_rec2 = predict(xx$models[[1]], newdata = d))
  }) %>% 
  map_dfr(rbind, .id = "model") %>% 
  filter(n_rec2 >= 0) %>% 
  separate(model, "Species", sep = " ~ ", extra = "drop", remove = F)


graph_sp <- ggplot(df_observed_sp, aes(mid_pred, n_rec2, color = Species)) + 
  geom_point() + 
  geom_line(data = df_predicted_sp) + 
  facet_wrap(~model,
             ncol = 3,
             labeller = labeller(model = function(x) {
               before <- str_extract(x, "^[^~]+")
               after  <- str_extract(x, "(?<=~).*")
               # Используем str_wrap для разбивки текста на строки
               wrapped_text <- str_wrap(after, width = 35)  # Ширина в символах
               paste(wrapped_text, paste0("Вид: ", before), sep = "\n")
             }),
             scales = "free") +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  labs(x = NULL, y = "Доля находок, %", 
       color = "Виды") + 
  theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "#dafdce"),
        strip.text = element_text(size = 9))

ggsave("plotspfinal.pdf", plot = graph_sp, device = cairo_pdf, height = 13.5, width = 9)

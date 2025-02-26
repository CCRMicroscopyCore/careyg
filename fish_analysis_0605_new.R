# New FISH analysis

# 02/14/25
# Author: Andy D. Tran, CCR Microscopy Core, LCBG, CCR, NCI

# Libraries and themes---------------------------------------------------------

library(tidyverse)
library(ggsci)

theme <- theme(panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               panel.background=element_blank(),
               axis.line=element_line(color="black", linewidth=1),
               axis.ticks=element_line(color="black", linewidth=1),
               text=element_text(size=18),
               plot.title=element_text(size=24),
               axis.text=element_text(color="black"),
               plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Define paths-----------------------------------------------------------------

base_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Grace Carey/output/06_05_Imaris'

output_path <- '/Volumes/LECIMAGE/Analysis/[NCI] [LGI] Sergio Ruiz/Grace Carey/results/06_05_Imaris'

# Download data----------------------------------------------------------------
### FISH spots-----------------------------------------------------------------
fish_list <- list.files(base_path, pattern = '_fish_Statistics')

df_tmp <- tibble(image = character(),
                 fish_id = integer(),
                 cell_id = integer()
                 )

for(fish_name in fish_list){
  fish_path <- file.path(base_path, fish_name)
  csv_list <- list.files(fish_path, pattern = '_Median_Ch=4_Img=1.csv')
  csv_path <- file.path(fish_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(fish_name, '_fish_Statistics', '')) %>% 
    rename(fish_id = contains('ID'),
           cell_id = contains('Median')) %>% 
    #mutate(fish_id = as.integer(),
    #       cell_id = as.integer()) %>% 
    select(image, fish_id, cell_id)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_fish_cellid <- df_tmp

df_tmp <- tibble(image = character(),
                 fish_id = integer(),
                 pos_x = numeric(),
                 pos_y = numeric(),
                 pos_z = numeric())

for(fish_name in fish_list){
  fish_path <- file.path(base_path, fish_name)
  csv_list <- list.files(fish_path, pattern = '_Position.csv')
  csv_path <- file.path(fish_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(fish_name, '_fish_Statistics', '')) %>% 
    rename(fish_id = contains('ID'),
           pos_x = contains('Position X'),
           pos_y = contains('Position Y'),
           pos_z = contains('Position Z')) %>% 
    select(image, fish_id, pos_x, pos_y, pos_z)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_fish_position <- df_tmp

df_tmp <- tibble(image = character(),
                 fish_id = integer(),
                 intensity = numeric())

for(fish_name in fish_list){
  fish_path <- file.path(base_path, fish_name)
  csv_list <- list.files(fish_path, pattern = '_Intensity_Mean_Ch=3_Img=1.csv')
  csv_path <- file.path(fish_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(fish_name, '_fish_Statistics', '')) %>% 
    rename(fish_id = contains('ID'),
           intensity = contains('Mean')) %>% 
    select(image, fish_id, intensity)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_fish <- full_join(df_fish_cellid, df_fish_position) %>%
  full_join(df_tmp) %>% 
  mutate(fish_id = as.integer(fish_id),
         cell_id = as.integer(cell_id)) %>% 
  filter(cell_id > 0) %>% 
  mutate(fish_id = fish_id + 1)

any(is.na(df_fish))
any(df_fish$fish_id == 0)

df_fish <- df_fish %>% 
  rename(fish_pos_x = pos_x,
         fish_pos_y = pos_y,
         fish_pos_z = pos_z,
         fish_intensity = intensity)

### Nucleoli-------------------------------------------------------------------

nucleolin_list <- list.files(base_path, pattern = '_nucleolin_Statistics')

df_tmp <- tibble(image = character(),
                 nucleolin_id = integer(),
                 cell_id = integer()
          )

for(nucleolin_name in nucleolin_list){
  nucleolin_path <- file.path(base_path, nucleolin_name)
  csv_list <- list.files(nucleolin_path, pattern = '_Median_Ch=4_Img=1.csv')
  csv_path <- file.path(nucleolin_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(nucleolin_name, '_nucleolin_Statistics', '')) %>% 
    rename(nucleolin_id = contains('ID'),
           cell_id = contains('Median')) %>% 
    select(image, nucleolin_id, cell_id)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_nucleolin_cellid <- df_tmp

df_tmp <- tibble(image = character(),
                 nucleolin_id = integer(),
                 pos_x = numeric(),
                 pos_y = numeric(),
                 pos_z = numeric())

for(nucleolin_name in nucleolin_list){
  nucleolin_path <- file.path(base_path, nucleolin_name)
  csv_list <- list.files(nucleolin_path, pattern = '_Position.csv')
  csv_path <- file.path(nucleolin_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(nucleolin_name, '_nucleolin_Statistics', '')) %>% 
    rename(nucleolin_id = contains('ID'),
           pos_x = contains('Position X'),
           pos_y = contains('Position Y'),
           pos_z = contains('Position Z')) %>% 
    select(image, nucleolin_id, pos_x, pos_y, pos_z)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_nucleolin_position <- df_tmp

df_tmp <- tibble(image = character(),
                 nucleolin_id = integer(),
                 intensity = numeric())

for(nucleolin_name in nucleolin_list){
  nucleolin_path <- file.path(base_path, nucleolin_name)
  csv_list <- list.files(nucleolin_path, pattern = '_Intensity_Mean_Ch=2_Img=1.csv')
  csv_path <- file.path(nucleolin_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(nucleolin_name, '_nucleolin_Statistics', '')) %>% 
    rename(nucleolin_id = contains('ID'),
           intensity = contains('Mean')) %>% 
    select(image, nucleolin_id, intensity)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_nucleolin_intensity <- df_tmp

df_tmp <- tibble(image = character(),
                 nucleolin_id = integer(),
                 volume = numeric())

for(nucleolin_name in nucleolin_list){
  nucleolin_path <- file.path(base_path, nucleolin_name)
  csv_list <- list.files(nucleolin_path, pattern = 'Volume.csv')
  csv_path <- file.path(nucleolin_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(nucleolin_name, '_nucleolin_Statistics', '')) %>% 
    rename(nucleolin_id = contains('ID'),
           volume = contains('Volume')) %>% 
    select(image, nucleolin_id, volume)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_nucleolin <- full_join(df_nucleolin_cellid, df_nucleolin_position) %>% 
  full_join(df_nucleolin_intensity) %>% 
  full_join(df_tmp) %>% 
  mutate(nucleolin_id = as.integer(nucleolin_id),
         cell_id = as.integer(cell_id)) %>% 
  filter(cell_id > 0) %>% 
  mutate(nucleolin_id = nucleolin_id + 1)

any(is.na(df_nucleolin))
any(df_nucleolin$nucleolin_id == 0)

df_nucleolin <- df_nucleolin %>% 
  rename(nuc_pos_x = pos_x,
         nuc_pos_y = pos_y,
         nuc_pos_z = pos_z,
         nuc_intensity = intensity)

### Cells----------------------------------------------------------------------

cell_list <- list.files(base_path, pattern = '_dapi_Statistics')

df_tmp <- tibble(image = character(),
                 temp_id = integer(),
                 cell_id = integer())

for(cell_name in cell_list){
  cell_path <- file.path(base_path, cell_name)
  csv_list <- list.files(cell_path, pattern = '_Median_Ch=4_Img=1.csv')
  csv_path <- file.path(cell_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(cell_name, '_dapi_Statistics', '')) %>% 
    rename(temp_id = contains('ID'),
           cell_id = contains('Median')) %>% 
    select(image, temp_id, cell_id)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_cell_cellid <- df_tmp

df_tmp <- tibble(image = character(),
                 temp_id = integer(),
                 volume = numeric())

for(cell_name in cell_list){
  cell_path <- file.path(base_path, cell_name)
  csv_list <- list.files(cell_path, pattern = '_Volume.csv')
  csv_path <- file.path(cell_path, csv_list[1])
  print(csv_list[1])
  csv_tmp <- read_csv(csv_path, skip = 3) %>% 
    mutate(image = str_replace(cell_name, '_dapi_Statistics', '')) %>% 
    rename(temp_id = contains('ID'),
           volume = 'Volume') %>% 
  select(image, temp_id, volume)
  df_tmp <- full_join(df_tmp, csv_tmp)
}

df_cell <- full_join(df_cell_cellid, df_tmp) %>% 
  filter(cell_id > 0)

any(is.na(df_cell))

df_cell <- df_cell %>% 
  rename(cell_volume = volume) %>% 
  select(image, cell_id, cell_volume)

### Combine data---------------------------------------------------------------

df_fish2 <- full_join(df_cell, df_fish) %>% 
  mutate(cell_line = case_when(
    grepl('E14', image) ~ 'E14',
    grepl('529', image) ~ '529',
    grepl('CMV184', image) ~ 'CMV184',
    grepl('CMV45', image) ~ 'CMV45',
    grepl('CMV118', image) ~ 'CMV118',
    grepl('CMV60', image) ~ 'CMV60'
  )) %>% 
  mutate(IAA = case_when(
    grepl('ctrl', image) ~ 'Ctrl',
    grepl('24h', image) ~ '24h',
    grepl('48h', image) ~ '48h'
  )) %>% 
  filter(cell_volume > 100)

df_fish3 <- df_fish2 %>% 
  group_by(image, cell_id) %>% 
  arrange(desc(fish_intensity)) %>% 
  slice(1:2) %>% 
  ungroup()

df_fish_test <- df_fish3 %>% 
  group_by(image, cell_line, IAA, cell_id) %>% 
  mutate(fish_count = ifelse(!is.na(fish_id), 1, 0)) %>% 
  summarise(cell_volume = mean(cell_volume),
            fish_count = sum(fish_count))

max(df_fish_test$fish_count)

df_fish_screen <- df_fish3 %>% 
  filter(!is.na(fish_id)) %>% 
  group_by(image, cell_id) %>% 
  summarise()

df_nucleolin2 <- full_join(df_cell, df_nucleolin) %>% 
  mutate(cell_line = case_when(
    grepl('E14', image) ~ 'E14',
    grepl('529', image) ~ '529',
    grepl('CMV184', image) ~ 'CMV184',
    grepl('CMV45', image) ~ 'CMV45',
    grepl('CMV118', image) ~ 'CMV118',
    grepl('CMV60', image) ~ 'CMV60'
  )) %>% 
  mutate(IAA = case_when(
    grepl('ctrl', image) ~ 'Ctrl',
    grepl('24h', image) ~ '24h',
    grepl('48h', image) ~ '48h'
  )) %>% 
  filter(cell_volume > 100)

df_nucleolin_screen <- df_nucleolin2 %>% 
  filter(!is.na(nucleolin_id)) %>% 
  group_by(image, cell_id) %>% 
  summarise() 

# Measure shortest distance to nucleoli----------------------------------------

distance_calc <- function(x1, y1, z1, x2, y2, z2){
  x_dist <- abs(x2 - x1)
  y_dist <- abs(y2 - y1)
  z_dist <- abs(z2 - z1)
  dist <- sqrt(x_dist^2 + y_dist^2 + z_dist^2)
  return (dist)
}

find_shortest_distance <- function(point, tmp_nuc_df){
  distances <- pmap_dbl(list(tmp_nuc_df$nuc_pos_x, tmp_nuc_df$nuc_pos_y, 
                             tmp_nuc_df$nuc_pos_z),
                    ~ distance_calc(point$fish_pos_x, point$fish_pos_y, point$fish_pos_z,
                                    ..1, ..2, ..3))
  return(min(distances))
}

df_cell_screen <- full_join(df_fish_screen, df_nucleolin2) %>% 
  filter(!is.na(nucleolin_id)) %>% 
  group_by(image, cell_id) %>% 
  summarise()

df_fish_dist <- full_join(df_cell_screen, df_fish3) %>% 
  filter(!is.na(fish_id))
df_nucleolin_dist <- full_join(df_cell_screen, df_nucleolin2) %>% 
  filter(!is.na(nucleolin_id))

df_dist <- tibble(image = character(),
                  cell_id = integer(),
                  cell_volume = numeric(),
                  fish_id = integer(),
                  fish_pos_x = numeric(),
                  fish_pos_y = numeric(),
                  fish_pos_z = numeric(),
                  fish_intensity = numeric(),
                  cell_line = character(),
                  IAA = character(),
                  shortest_distance = numeric())
image_list <- unique(df_fish_dist$image)

for(image_name in image_list){
  print(image_name)
  df_fish_image <- df_fish_dist %>% 
    filter(image == image_name)
  df_nucleolin_image <- df_nucleolin_dist %>% 
    filter(image == image_name)
  cell_list <- unique(df_fish_image$cell_id)
  for(cell in cell_list){
    print(cell)
    df_fish_tmp <- df_fish_image %>% 
      filter(cell_id == cell)
    df_nucleolin_tmp <- df_nucleolin_image %>% 
      filter(cell_id == cell)
    df_dist_tmp <- df_fish_tmp %>% 
      rowwise() %>% 
      mutate(shortest_distance = find_shortest_distance(cur_data(), 
                                                        df_nucleolin_tmp)) %>% 
      ungroup()
    df_dist <- full_join(df_dist, df_dist_tmp)
  }
}

df_dist <- df_dist %>% 
  filter(shortest_distance <= 10) %>% 
  filter(!is.infinite(shortest_distance))

df_dist$cell_line <- factor(df_dist$cell_line, levels = c('E14',
                                                          '529',
                                                          'CMV184',
                                                          'CMV45',
                                                          'CMV118',
                                                          'CMV60'))
df_dist$IAA <- factor(df_dist$IAA, levels = c('Ctrl',
                                              '24h',
                                              '48h'))

unique(df_dist$cell_line)
unique(df_dist$IAA)

ggplot(df_dist, aes(x = IAA, y = shortest_distance, fill = cell_line)) +
  geom_violin(width = 0.5, col = 'black', alpha = 0.75, 
              show.legend = F) +
  geom_boxplot(width = 0.1, col = 'black', alpha = 1,
               show.legend = F, outlier.size = -1) +
  facet_grid(cols = vars(cell_line)) +
  theme +
  labs(title = 'FISH spot min distance to heterochromatin', x = 'IAA treatment',
       y = 'Distance [um]')

plotname <- 'FISH_spot_heterochromatin_distance.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16, 
       height = 8, dpi = 150)

# Plot nucleolin intensity-----------------------------------------------------

df_nucleolin_plot <- df_nucleolin2 %>% 
  filter(!is.na(nucleolin_id))

df_nucleolin_plot$cell_line <- factor(df_nucleolin_plot$cell_line, 
                                      levels = c('E14',
                                                 '529',
                                                 'CMV184',
                                                 'CMV45',
                                                 'CMV118',
                                                 'CMV60'))

df_nucleolin_plot$IAA <- factor(df_nucleolin_plot$IAA, 
                                levels = c('Ctrl',
                                           '24h',
                                           '48h'))

unique(df_nucleolin_plot$cell_line)
unique(df_nucleolin_plot$IAA)

ggplot(df_nucleolin_plot, aes(x = IAA, y = nuc_intensity, fill = cell_line)) +
  geom_violin(width = 0.5, col = 'black', alpha = 0.75, 
              show.legend = F) +
  geom_boxplot(width = 0.1, col = 'black', alpha = 1,
               show.legend = F, outlier.size = -1) +
  facet_grid(cols = vars(cell_line)) +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = 'l') +
  theme +
  labs(title = 'Nucleoli intensity', x = 'IAA treatment',
       y = 'Intensity [AU]')

plotname <- 'Nucleoli_intensity.png'
ggsave(plotname, plot = last_plot(), path = output_path, width = 16, 
       height = 8, dpi = 150)

# Export data------------------------------------------------------------------

df_out <- df_dist %>% 
  select(image, cell_line, IAA, cell_id, fish_id, fish_intensity, shortest_distance)

csvname <- 'fish_distance_nucleoli.csv'
csvpath <- file.path(output_path, csvname)
write_csv(df_out, csvpath)

df_out <- df_nucleolin_plot %>% 
  select(image, cell_line, IAA, cell_id, nucleolin_id, nuc_intensity, volume)

csvname <- 'nucleoli_intensity.csv'
csvpath <- file.path(output_path, csvname)
write_csv(df_out, csvpath)

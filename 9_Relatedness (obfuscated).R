library(flextable)
library(dplyr)
library(here)
library(PNWColors)
library(ggplot2)
library(patchwork)
library(sf)
library(MMWRweek)
library(stringr)
source(here::here("From_to_functions.R"))
source(here::here("Creating_networks.R"))

Tobfuscate <- readRDS(here::here("Objects/Talent obfuscate.rds"))
Aobfuscate <- readRDS(here::here("Objects/Ashland obfuscate.rds"))
append <- zzz #PII.. obfuscated
Aobfuscate <- rbind(Aobfuscate, append)

append <- data.frame(Tsite = "Wagner Creek", Tlabel = "T29")
Tobfuscate <- rbind(Tobfuscate, append)

ash <- readRDS(here::here("Objects/kinship/ash full graph.RDS")) %>%
  gene_dist() %>%
  dplyr::mutate(
    Dist = as.numeric(Dist),
    To_site = zzz, # PII... obfuscated
    From_site = zzz # PII... obfuscated
  )
tal <- readRDS(here::here("Objects/kinship/tal full graph.RDS")) %>%
  gene_dist() %>%
  dplyr::mutate(
    Dist = as.numeric(Dist),
    date_diff = as.numeric(date_diff),
    epi_wk = MMWRweek(To_Date)$MMWRweek
  )

tal <- tal %>%
  left_join(Tobfuscate, by = c("From_site" = "Tsite")) %>%
  dplyr::rename(From_label = "Tlabel") %>%
  left_join(Tobfuscate, by = c("To_site" = "Tsite")) %>%
  dplyr::rename(To_label = "Tlabel") %>%
  dplyr::mutate(
    From_obf = paste(From_label, str_extract(From, "_[0-9]{0,2}[A-Z]{0,1}$")),
    To_obf = paste(To_label, str_extract(To, "_[0-9]{0,2}[A-Z]{0,1}$")),
  )

ash <- ash %>%
  left_join(Aobfuscate, by = c("From_site" = "Asite")) %>%
  dplyr::rename(From_label = "Alabel") %>%
  left_join(Aobfuscate, by = c("To_site" = "Asite")) %>%
  dplyr::rename(To_label = "Alabel") %>%
  dplyr::mutate(
    From_obf = paste(From_label, str_extract(From, "_[0-9]{1,2}[A-Z]{0,1}")),
    To_obf = paste(To_label, str_extract(To, "_[0-9]{1,2}[A-Z]{0,1}")),
  )

tal_types <- as.data.frame(table(tal$From_obf, tal$Type))
names(tal_types) <- c("Individual", "Relatedness", "Count")

tal_types <- tal_types %>%
  filter(Relatedness != "MZ") %>%
  mutate(T_num = as.numeric(str_extract(Individual, "(?<=T)\\d+"))) %>%
  arrange(T_num, Individual) %>%
  mutate(Individual = factor(Individual, levels = unique(Individual)))


ash_types <- as.data.frame(table(ash$From_obf, ash$Type))
names(ash_types) <- c("Individual", "Relatedness", "Count")

ash_types <- ash_types %>%
  filter(Relatedness != "MZ") %>%
  mutate(A_num = as.numeric(str_extract(Individual, "(?<=T)\\d+"))) %>%
  arrange(A_num, Individual) %>%
  mutate(Individual = factor(Individual, levels = unique(Individual)))

ash_count_sum <- ash_types %>% 
  group_by(Individual) %>% 
  dplyr::mutate(Count = sum(Count))

tal_count_sum <- tal_types %>% 
  group_by(Individual) %>% 
  dplyr::mutate(Count = sum(Count))

y_max <- max(ash_count_sum$Count, tal_count_sum$Count)

# trying to keep these colors in line with JRs descriptive figures...
my_palette <- c("First" = "#E8B84B", 
                "Second" = "#F2C6D0",
                "Third" = "#4B2E83")

p1 <- ggplot(tal_types, aes(x = Individual, y = Count, fill = Relatedness)) +
  geom_col() +
  scale_fill_manual(values = my_palette) +
  scale_y_continuous(limits = c(0, y_max)) +
  labs(x = NULL, y = "Count", fill = "Close kin category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(ash_types, aes(x = Individual, y = Count, fill = Relatedness)) +
  geom_col() +
  scale_fill_manual(values = my_palette) +
  scale_y_continuous(limits = c(0, y_max)) +
  labs(x = NULL, y = NULL, fill = "Close kin category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank())

p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

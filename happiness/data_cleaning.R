# Code to transform raw data from WHR or WDI

WHR <- read.csv("happiness/WHR.csv")
WDI <- read.csv("happiness/WDI.csv")
ladder <- read.csv("happiness/ladder.csv")

# WDI
df_avg <- WDI %>%
  mutate(across(starts_with("X"), as.numeric)) %>% # Convert all year columns to numeric
  rowwise() %>%
  mutate(Average = mean(c_across(starts_with("X")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(-starts_with("X")) %>%
  select(-c(Country.Code, Series.Code)) %>%
  filter(!(is.na(Series.Name) | Series.Name == "")) %>%
  slice(1:1302) %>%
  pivot_wider(
    names_from = Series.Name, 
    values_from = Average     
  )

write_csv(df_avg, "happiness/WDI.csv")


# WHR
df_avg <- WHR %>%
  filter(year >= 2019 & year <= 2023) %>% 
  select(-tail(names(.), 2)) %>%
  select(-year) %>%
  group_by(Country.name) %>%
  summarize(across(2:7, \(x) mean(x, na.rm = TRUE)), .groups = 'drop') %>%
  rename(Country.Name = Country.name)

write_csv(df_avg, "happiness/WHR.csv")


# Combined_dat
Combined_dat <- inner_join(WHR, WDI)
Combined_dat <- Combined_dat %>% 
  filter(!if_any(-c(9), is.na)) %>%
  select(-9) %>%
  select(-c(Generosity, Life.expectancy.at.birth..total..years.)) %>%
  mutate(Log.Scientific.and.technical.journal.articles = log(Scientific.and.technical.journal.articles + 1)) %>%
  select(-Scientific.and.technical.journal.articles)

ladder <- ladder %>%
  rename(Country.Name = Country.name) %>%
  select(c(Country.Name, Ladder.score))
Combined_dat <- inner_join(Combined_dat, ladder)

write_csv(Combined_dat, "happiness/Combined_dat.csv")


# Scaled_dat
min_vals <- Combined_dat %>%
  select(2:10) %>%
  summarise(across(everything(), min, na.rm = TRUE)) %>%
  unlist()

max_vals <- Combined_dat %>%
  select(2:10) %>%
  summarise(across(everything(), max, na.rm = TRUE)) %>%
  unlist()

Scaled_Centered <- Combined_dat %>% 
  mutate(across(2:10, ~ (. - min_vals[cur_column()]) / (max_vals[cur_column()] - min_vals[cur_column()]))) %>%
  mutate(Log.GDP.per.capita = Log.GDP.per.capita - mean(Log.GDP.per.capita, na.rm = TRUE))
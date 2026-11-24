# --Bg value ---
df_sero_sub[, c("IgA-MFI","IgG-MFI","SiteIII-MPE8", "Site II-Motavizumab",
                "Site V-Suptavumab", "Site Φ-Nirsevimab",
                "Site IV-Clesrovimab")] %>%
  mutate(across(ends_with('MFI'),~ifelse(.x>1,.x,NA))) %>%
  mutate(across(starts_with('Site'),~ifelse(.x>0.01,.x,NA))) %>%
  map_dfr(~quantile(.x,c(1:10/10),na.rm = T))

df_sero_sub %>%
  select(`IgA-MFI`, `IgG-MFI`, `SiteIII-MPE8`, `Site II-Motavizumab`,
         `Site V-Suptavumab`, `Site Φ-Nirsevimab`, `Site IV-Clesrovimab`) %>%
  mutate(across(ends_with('MFI'),~ifelse(.x>1,.x,NA))) %>%
  mutate(across(starts_with('Site'),~ifelse(.x>0.01,.x,NA))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  facet_wrap(~Variable, scales = "free", ncol = 3) +
  scale_x_log10() +  # Log transform due to wide range
  labs(title = "Distribution of Antibody/Antigen Measurements",
       x = "Value (log10 scale)", y = "Count") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8, face = "bold"))

BG_VALUE<-df_sero_sub[, c("IgA-MFI","IgG-MFI","SiteIII-MPE8", "Site II-Motavizumab",
                          "Site V-Suptavumab", "Site Φ-Nirsevimab",
                          "Site IV-Clesrovimab")] %>%
  mutate(across(ends_with('MFI'),~ifelse(.x>1,.x,NA))) %>%
  mutate(across(starts_with('Site'),~ifelse(.x>0.01,.x,NA))) %>%
  map_dbl(~quantile(.x,c(0.1),na.rm = T))

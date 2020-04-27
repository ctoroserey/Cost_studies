
summaryOC <- list()
summaryOC$all <- dataWth %>%
  mutate(Block = case_when(
    Block %in% c(1, 2) ~ 1,
    Block %in% c(3, 4) ~ 2,
    Block %in% c(5, 6) ~ 3
  )) %>%
  group_by(SubjID, Cost, Half) %>% #group_by(SubjID, Cost, Block) %>%
  do(optimizeOCModel(., simplify = T)) %>%
  ungroup()

# plot
(summaryOC$plot <- ggplot(summaryOC$all, aes(Cost, Gamma, fill = Cost)) +
  geom_hline(yintercept = 0.7, alpha = 0.9, color = "gray40", size = 1, linetype = "dashed") +
  geom_jitter(pch = 21, size = 3, show.legend = F) +
  ylim(0, 1.5) +
  labs(x = "") +
  scale_fill_manual(values = colsWth) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 16)))

summaryOC$all %>%
  select(SubjID, Cost, Gamma, Half) %>% 
  spread(Half, Gamma) %>%
  ggplot(aes(Half_1, Half_2, fill = Cost)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(pch = 21, color = "black", size = 3) +
    ylim(0, 2) +
    xlim(0, 2) +
    theme_minimal()
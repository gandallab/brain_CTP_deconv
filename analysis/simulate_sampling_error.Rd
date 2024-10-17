library(ggplot2)
library(cowplot)
library(gridGraphics)
library(tidyverse)

set.seed(1234)

# Create 2 groups (n=50 vs n=1000) with a true mean difference of 0.05
#size=50
size=1000
tab <- data.frame(case = rnorm(size, mean = 0.50, sd = 0.05), ctrl = rnorm(size, mean = 0.55, sd = 0.05))

t.test(tab$case, tab$ctrl) # t-test is significant

# Plots or not
plots <- "no"
# plots <- "yes"

# Simulate n=100 iterations

## Output data.frame
if (size == 50) {
	sim_out_50 <- data.frame(cor_case_sd05 = as.numeric(), cor_ctrl_sd05 = as.numeric(), tp_sd05 = as.numeric(), cor_case_sd10 = as.numeric(), cor_ctrl_sd10 = as.numeric(), tp_sd10 = as.numeric(), cor_case_sd20 = as.numeric(), cor_ctrl_sd20 = as.numeric(), tp_sd20 = as.numeric())
} else if (size == 1000) {
	sim_out_1000 <- data.frame(cor_case_sd05 = as.numeric(), cor_ctrl_sd05 = as.numeric(), tp_sd05 = as.numeric(), cor_case_sd10 = as.numeric(), cor_ctrl_sd10 = as.numeric(), tp_sd10 = as.numeric(), cor_case_sd20 = as.numeric(), cor_ctrl_sd20 = as.numeric(), tp_sd20 = as.numeric())
}


# Generate sampling variation

for (i in 1:100) {
	
print(i)
	
## sd=0.05
tab$case_s1_sd05 <- apply(tab, 1, function(x) rnorm(1, mean = x[1], sd = 0.05))
tab$case_s2_sd05 <- apply(tab, 1, function(x) rnorm(1, mean = x[1], sd = 0.05))
tab$ctrl_s1_sd05 <- apply(tab, 1, function(x) rnorm(1, mean = x[2], sd = 0.05))
tab$ctrl_s2_sd05 <- apply(tab, 1, function(x) rnorm(1, mean = x[2], sd = 0.05))

cor_case_sd05 <- cor.test(tab$case_s1_sd05, tab$case_s2_sd05)$estimate
cor_ctrl_sd05 <- cor.test(tab$ctrl_s1_sd05, tab$ctrl_s2_sd05)$estimate
tp_sd05 <- t.test(tab$case_s1_sd05, tab$ctrl_s1_sd05)$p.value

if (plots == "yes") {
	p1 <- ggplot(tab, aes(x = case_s1_sd05, y = case_s2_sd05)) + geom_point() + geom_smooth(method = "lm") + ggtitle(paste("r =", round(cor_case_sd05, 2)))
	p2 <- ggplot(tab, aes(x = ctrl_s1_sd05, y = ctrl_s2_sd05)) + geom_point() + geom_smooth(method = "lm") + ggtitle(paste("r =", round(cor_ctrl_sd05, 2)))
	p3 <- tab %>% tidyr::pivot_longer(cols = c("case_s1_sd05", "ctrl_s1_sd05"), names_to = "group", values_to = "prop") %>% ggplot(., aes(x = group, y = prop)) + geom_boxplot(colour = c("red", "blue")) + geom_hline(yintercept = 0.5, colour = "red", linetype = "dotted") + geom_hline(yintercept = 0.55, colour = "blue", linetype = "dotted") + ggtitle(paste("p =", signif(tp_sd05, 3)))
}

## sd=0.10
tab$case_s1_sd10 <- apply(tab, 1, function(x) rnorm(1, mean = x[1], sd = 0.10))
tab$case_s2_sd10 <- apply(tab, 1, function(x) rnorm(1, mean = x[1], sd = 0.10))
tab$ctrl_s1_sd10 <- apply(tab, 1, function(x) rnorm(1, mean = x[2], sd = 0.10))
tab$ctrl_s2_sd10 <- apply(tab, 1, function(x) rnorm(1, mean = x[2], sd = 0.10))

cor_case_sd10 <- cor.test(tab$case_s1_sd10, tab$case_s2_sd10)$estimate
cor_ctrl_sd10 <- cor.test(tab$ctrl_s1_sd10, tab$ctrl_s2_sd10)$estimate
tp_sd10 <- t.test(tab$case_s1_sd10, tab$ctrl_s1_sd10)$p.value

if (plots == "yes") {
	p4 <- ggplot(tab, aes(x = case_s1_sd10, y = case_s2_sd10)) + geom_point() + geom_smooth(method = "lm") + ggtitle(paste("r =", round(cor_case_sd10, 2)))
	p5 <- ggplot(tab, aes(x = ctrl_s1_sd10, y = ctrl_s2_sd10)) + geom_point() + geom_smooth(method = "lm") + ggtitle(paste("r =", round(cor_ctrl_sd10, 2)))
	p6 <- tab %>% tidyr::pivot_longer(cols = c("case_s1_sd10", "ctrl_s1_sd10"), names_to = "group", values_to = "prop") %>% ggplot(., aes(x = group, y = prop)) + geom_boxplot(colour = c("red", "blue")) + geom_hline(yintercept = 0.5, colour = "red", linetype = "dotted") + geom_hline(yintercept = 0.55, colour = "blue", linetype = "dotted") + ggtitle(paste("p =", signif(tp_sd10, 3)))
}

## sd=0.20
tab$case_s1_sd20 <- apply(tab, 1, function(x) rnorm(1, mean = x[1], sd = 0.20))
tab$case_s2_sd20 <- apply(tab, 1, function(x) rnorm(1, mean = x[1], sd = 0.20))
tab$ctrl_s1_sd20 <- apply(tab, 1, function(x) rnorm(1, mean = x[2], sd = 0.20))
tab$ctrl_s2_sd20 <- apply(tab, 1, function(x) rnorm(1, mean = x[2], sd = 0.20))

cor_case_sd20 <- cor.test(tab$case_s1_sd20, tab$case_s2_sd20)$estimate
cor_ctrl_sd20 <- cor.test(tab$ctrl_s1_sd20, tab$ctrl_s2_sd20)$estimate
tp_sd20 <- t.test(tab$case_s1_sd20, tab$ctrl_s1_sd20)$p.value

if (plots == "yes") {
	p7 <- ggplot(tab, aes(x = case_s1_sd20, y = case_s2_sd20)) + geom_point() + geom_smooth(method = "lm") + ggtitle(paste("r =", round(cor_case_sd20, 2)))
	p8 <- ggplot(tab, aes(x = ctrl_s1_sd20, y = ctrl_s2_sd20)) + geom_point() + geom_smooth(method = "lm") + ggtitle(paste("r =", round(cor_ctrl_sd20, 2)))
	p9 <- tab %>% tidyr::pivot_longer(cols = c("case_s1_sd20", "ctrl_s1_sd20"), names_to = "group", values_to = "prop") %>% ggplot(., aes(x = group, y = prop)) + geom_boxplot(colour = c("red", "blue")) + geom_hline(yintercept = 0.5, colour = "red", linetype = "dotted") + geom_hline(yintercept = 0.55, colour = "blue", linetype = "dotted") + ggtitle(paste("p =", signif(tp_sd20, 3)))
}

if (size == 50) {
	tmp <- data.frame(cor.case_sd05 = cor_case_sd05, cor.ctrl_sd05 = cor_ctrl_sd05, tp_sd05 = tp_sd05, cor.case_sd10 = cor_case_sd10, cor.ctrl_sd10 = cor_ctrl_sd10, tp_sd10 = tp_sd10, cor.case_sd20 = cor_case_sd20, cor.ctrl_sd20 = cor_ctrl_sd20, tp_sd20 = tp_sd20)
	sim_out_50 <- rbind(sim_out_50, tmp)
} else if (size == 1000) {
	tmp <- data.frame(cor.case_sd05 = cor_case_sd05, cor.ctrl_sd05 = cor_ctrl_sd05, tp_sd05 = tp_sd05, cor.case_sd10 = cor_case_sd10, cor.ctrl_sd10 = cor_ctrl_sd10, tp_sd10 = tp_sd10, cor.case_sd20 = cor_case_sd20, cor.ctrl_sd20 = cor_ctrl_sd20, tp_sd20 = tp_sd20)
	sim_out_1000 <- rbind(sim_out_1000, tmp)

}

}

#=====================================================================

# Generate plots

## Scatterplots
if (size == 50) {
	sim50_plot <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)
	save_plot("~/Downloads/simulate_sampling_error_n50.png", sim50_plot, base_height = 6)
} else if (size == 1000) {
	sim1000_plot <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)
	save_plot("~/Downloads/simulate_sampling_error_n1000.png", sim1000_plot, base_height = 6)
}

## Histograms

sim_out_50.long <- pivot_longer(sim_out_50, cols = everything(), names_to = "statistic", values_to = "values") %>% mutate(comparison = sub("_.*", "", statistic)) %>% mutate(sd = sub(".*_", "", statistic)) %>% mutate(comparison = gsub("tp", "p_ttest", comparison)) %>% mutate(sample_size = 50)
sim_out_1000.long <- pivot_longer(sim_out_1000, cols = everything(), names_to = "statistic", values_to = "values") %>% mutate(comparison = sub("_.*", "", statistic)) %>% mutate(sd = sub(".*_", "", statistic)) %>% mutate(comparison = gsub("tp", "p_ttest", comparison)) %>% mutate(sample_size = 1000)

sim_out_agg.long <- rbind(sim_out_50.long, sim_out_1000.long)

sim_out_cor_hist <- ggplot(sim_out_agg.long %>% filter(comparison != "p_ttest"), aes(x = values)) + geom_histogram() + facet_grid(sd ~ comparison ~ sample_size, scales = "free") + ggtitle("Correlation coefficients (Pearson's), Simulations = 100") + theme_bw()
ggsave("~/Downloads/simulate_sampling_error_n50_rep100_cor.png", sim_out_cor_hist, height = 6, width = 6)

sim_out_p_hist <- ggplot(sim_out_agg.long %>% filter(comparison == "p_ttest"), aes(x = values)) + geom_histogram() + facet_grid(sd ~ comparison ~ sample_size, scales = "free") + ggtitle("t-test p-values, Simulations = 100") + theme_bw()
ggsave("~/Downloads/simulate_sampling_error_n50_rep100_pttest.png", sim_out_p_hist, height = 6, width = 6)


sim_out_50_hist <- ggplot(sim_out_50.long, aes(x = values)) + geom_histogram() + facet_grid(sd ~ comparison, scales = "free") + ggtitle("Sample size = 50, Simulations = 100") + theme_bw()
ggsave("~/Downloads/simulate_sampling_error_n50_rep100.png", sim_out_50_hist, height = 6, width = 6)
sim_out_1000_hist <- ggplot(sim_out_1000.long, aes(x = values)) + geom_histogram() + facet_grid(sd ~ comparison, scales = "free") + ggtitle("Sample size = 1000, Simulations = 100") + theme_bw()
ggsave("~/Downloads/simulate_sampling_error_n1000_rep100.png", sim_out_1000_hist, height = 6, width = 6)


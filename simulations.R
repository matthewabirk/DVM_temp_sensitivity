library(dplyr); library(ggplot2)

# Simulations -------------------------------------------------------------

dg_raw = readxl::read_xlsx('Dosidicus_data_final_2.xlsx', sheet = 'values')

dg_raw = dg_raw %>% group_by(temp, metric) %>% summarise(mean_mr = mean(value), se_mr = birk::se(value))

sims = apply(dg_raw, 1, function(i) rnorm(n = 1e5, mean = as.numeric(i['mean_mr']), sd = as.numeric(i['se_mr'])))
colnames(sims) = paste(dg_raw$metric, dg_raw$temp, sep = '_')

sims = cbind(sims, FAS_10 = sims[, 'mmr_10'] / sims[, 'smr_10'])
sims = cbind(sims, FAS_20 = sims[, 'mmr_20'] / sims[, 'smr_20'])
sims = cbind(sims, FAS_25 = sims[, 'mmr_25'] / sims[, 'smr_25'])

sims = cbind(sims, alpha_10 = sims[, 'smr_10'] / sims[, 'pcrit_smr_10'])
sims = cbind(sims, alpha_20 = sims[, 'smr_20'] / sims[, 'pcrit_smr_20'])
sims = cbind(sims, alpha_25 = sims[, 'smr_25'] / sims[, 'pcrit_smr_25'])

sims = cbind(sims, pcrit_mmr_10 = sims[, 'mmr_10'] / sims[, 'alpha_10'])
sims = cbind(sims, pcrit_mmr_20 = sims[, 'mmr_20'] / sims[, 'alpha_20'])
sims = cbind(sims, pcrit_mmr_25 = sims[, 'mmr_25'] / sims[, 'alpha_25'])

sim_a = ggplot() +
	geom_histogram(data = data.frame(fas = sims[, 'FAS_10']), aes(fas), fill = 'blue', alpha = 0.3) +
	geom_histogram(data = data.frame(fas = sims[, 'FAS_20']), aes(fas), fill = 'green', alpha = 0.3) +
	geom_histogram(data = data.frame(fas = sims[, 'FAS_25']), aes(fas), fill = 'red', alpha = 0.3) +
	theme_bw() +
	labs(x = 'FAS') +
	scale_x_continuous(breaks = seq(0, 10, by = 1)) +
	theme(axis.title.y = element_blank())

sim_b = ggplot() +
	geom_histogram(data = data.frame(alpha = sims[, 'alpha_10']), aes(alpha), fill = 'blue', alpha = 0.3) +
	geom_histogram(data = data.frame(alpha = sims[, 'alpha_20']), aes(alpha), fill = 'green', alpha = 0.3) +
	geom_histogram(data = data.frame(alpha = sims[, 'alpha_25']), aes(alpha), fill = 'red', alpha = 0.3) +
	scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
	theme_bw() +
	labs(x = 'Alpha (umol/g/hr/kPa)') +
	theme(axis.title.y = element_blank())

sim_c = ggplot() +
	geom_histogram(data = data.frame(pcrit_mmr = sims[, 'pcrit_mmr_10']), aes(pcrit_mmr), fill = 'blue', alpha = 0.3) +
	geom_histogram(data = data.frame(pcrit_mmr = sims[, 'pcrit_mmr_20']), aes(pcrit_mmr), fill = 'green', alpha = 0.3) +
	geom_histogram(data = data.frame(pcrit_mmr = sims[, 'pcrit_mmr_25']), aes(pcrit_mmr), fill = 'red', alpha = 0.3) +
	theme_bw() +
	scale_x_continuous(breaks = seq(0, 40, by = 5)) +
	labs(x = 'Pcrit-MMR (kPa)') +
	theme(axis.title.y = element_blank())

ggsave('simulations.pdf', cowplot::plot_grid(sim_a, sim_b, sim_c), width = 8, height = 6.5)







# Monte Carlo mapping -----------------------------------------------------

d_roi = readRDS('WOA/d_roi.rds')

source_partial <- function(file, start, end, ...) {
	file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
	file.lines.collapsed <- paste(file.lines, collapse='\n')
	source(textConnection(file.lines.collapsed), ...)
}




######### wrap all of this in a loop through sims and append to each other

n_sims = 5000 # THIS WILL TAKE A LONG TIME TO RUN. I HAD TO RUN ON SERVER WITH A FEW HUNDRED GB RAM TO PROCESS IN REASONABLE TIME. CAN TEST OUT WITH n_sims = 10.

d_sims = matrix(nrow = nrow(d_transect) * n_sims, ncol = 15)

for(i in 1:n_sims){
	print(i)
	dg = data.frame(temp = c(10, 20, 25),
									smr = c(sims[i, 'smr_10'], sims[i, 'smr_20'], sims[i, 'smr_25']),
									mmr = c(sims[i, 'mmr_10'], sims[i, 'mmr_20'], sims[i, 'mmr_25']),
									pcrit_smr = c(sims[i, 'pcrit_smr_10'], sims[i, 'pcrit_smr_20'], sims[i, 'pcrit_smr_25']),
									row.names = NULL)
	
	source_partial('Dosidicus changing E_MMR.R', 12, 57)
	source_partial('run_transects.R', 1, 21)
	
	tmp = select(d_transect, -one_of('lon', 'lon0.5', 'lat0.5', 'o2', 'km_per_deg_lon', 'boundary', 'coast', 'roi_east', 'roi_west', 'MI', 'habitable', 'trust'))
	
	d_sims[(nrow(tmp) * (i - 1) + 1):(nrow(tmp) * (i)), ] = as.matrix(tmp)
}



df_sims = as.data.frame(d_sims)
colnames(df_sims) = colnames(tmp)
df_sims_2.5 = df_sims %>% group_by(lat, depth) %>% summarise_all(quantile, 0.025, na.rm = TRUE)
df_sims_97.5 = df_sims %>% group_by(lat, depth) %>% summarise_all(quantile, 0.975, na.rm = TRUE)
df_sims_diff = (df_sims_97.5 - df_sims_2.5)
df_sims_diff[, c('lat', 'depth', 'temp', 'po2', 'height')] = df_sims_2.5[, c('lat', 'depth', 'temp', 'po2', 'height')]

library(birk); library(metR)







# Make top portion of Fig S5 -------------------------------------------------------------


p2.5 = ggplot(df_sims_2.5, aes(lat, depth, height = height)) +
	labs(x = 'Latitude', y = 'Depth (m)') +
	scale_y_reverse(expand = c(0, 0)) +
	scale_x_latitude(ticks = 20) +
	geom_tile(aes(fill = FAS_min)) + 
	geom_contour(aes(z = FAS_min), col = 'white') + 
	geom_text_contour(aes(z = FAS_min), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + 
	labs(title = expression(paste('2.5% quantile of ', FAS[min]))) + 
	theme(legend.title = element_blank()) + 
	scale_fill_divergent(midpoint = 3)

p_real = ggplot(d_transect, aes(lat, depth, height = height)) +					# FROM make_Fig4.R
	labs(x = 'Latitude', y = 'Depth (m)') +
	scale_y_reverse(expand = c(0, 0)) +
	scale_x_latitude(ticks = 20) +
	geom_tile(aes(fill = FAS_min)) + 
	geom_contour(aes(z = FAS_min), col = 'white') + 
	geom_text_contour(aes(z = FAS_min), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + 
	labs(title = expression(paste('Modelled ', FAS[min]))) + 
	theme(legend.title = element_blank()) + 
	scale_fill_divergent(midpoint = 3) +
	theme(axis.title.y = element_blank())


p97.5 = ggplot(df_sims_97.5, aes(lat, depth, height = height)) +
	labs(x = 'Latitude', y = 'Depth (m)') +
	scale_y_reverse(expand = c(0, 0)) +
	scale_x_latitude(ticks = 20) +
	geom_tile(aes(fill = FAS_min)) + 
	geom_contour(aes(z = FAS_min), col = 'white') + 
	geom_text_contour(aes(z = FAS_min), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + 
	labs(title = expression(paste('97.5% quantile of ', FAS[min]))) + 
	theme(legend.title = element_blank()) + 
	scale_fill_divergent(midpoint = 3) +
	theme(axis.title.y = element_blank())


s5_top = cowplot::plot_grid(p2.5, p_real, p97.5, nrow = 1)


# Make bottom portion of Fig S5 -------------------------------------------




p = ggplot(df_sims_diff, aes(lat, depth, height = height)) +
	labs(x = 'Latitude', y = 'Depth (m)') +
	scale_y_reverse(expand = c(0, 0)) +
	scale_x_latitude(ticks = 20) +
	scale_alpha_discrete(range = c(0.6, 1), guide = 'none') +
	theme(axis.title.y = element_blank())


p1 = p + geom_tile(aes(fill = SMR)) + geom_contour(aes(z = SMR), col = 'white') + geom_text_contour(aes(z = SMR), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + labs(title = expression(paste('95% CI width of SMR (', mu, 'mol/g/hr)'))) + theme(legend.title = element_blank())
p2 = p + geom_tile(aes(fill = Pcrit)) + geom_contour(aes(z = Pcrit), col = 'white') + geom_text_contour(aes(z = Pcrit), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + labs(title = expression(paste('95% CI width of ', Pc[SMR], ' (kPa)'))) + theme(legend.title = element_blank())
p3 = p + geom_tile(aes(fill = alpha)) + geom_contour(aes(z = alpha), col = 'white') + geom_text_contour(aes(z = alpha), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + labs(title = expression(paste('95% CI width of ', alpha, ' (', mu, 'mol/g/hr/kPa)'))) + theme(legend.title = element_blank())
p4 = p + geom_tile(aes(fill = MMR_min)) + geom_contour(aes(z = MMR_min), col = 'white') + geom_text_contour(aes(z = MMR_min), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + labs(title = expression(paste('95% CI width of ', MMR[min], ' (', mu, 'mol/g/hr)'))) + theme(legend.title = element_blank())
p5 = p + geom_tile(aes(fill = Pcrit_MMR)) + geom_contour(aes(z = Pcrit_MMR), col = 'white') + geom_text_contour(aes(z = Pcrit_MMR), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + labs(title = expression(paste('95% CI width of ', Pc[max], ' (kPa)'))) + theme(legend.title = element_blank())
p6 = p + geom_tile(aes(fill = FAS_min)) + geom_contour(aes(z = FAS_min), col = 'white') + geom_text_contour(aes(z = FAS_min), stroke = 0.2, rotate = FALSE, check_overlap = TRUE) + labs(title = expression(paste('95% CI width of ', FAS[min]))) + theme(legend.title = element_blank())

cp = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)

ggsave('Dg_sims_95CI.pdf', cp, width = 7, height = 7)

saveRDS(df_sims_2.5, 'Map simulations/df_sims_2.5.rds')
saveRDS(df_sims_97.5, 'Map simulations/df_sims_97.5.rds')

s5 = cowplot::plot_grid(s5_top, cp, ncol = 1, rel_heights = c(0.25, 0.75))
ggsave('Fig S5.pdf', s5, width = 10, height = 10)

d_roi = readRDS('WOA/d_roi.rds')

# ALTERNATIVELY: source('Make d_roi.R')

# Add physiology data -----------------------------------------------------

tables1 = readxl::read_xlsx('table_S1 ETP only (updated 2021-04-09).xlsx')

taxon = filter(tables1, species == .GlobalEnv$species) %>% select(temp, smr, pcrit_smr, alpha)

taxonE = filter(tables1, species == .GlobalEnv$species) %>% select(E_SMR, E_PcSMR, E_alpha)
taxonE = as.numeric(as.data.frame(taxonE[1, ]))
names(taxonE) = c('smr', 'pcrit_smr', 'alpha')



d_roi$SMR = adj_by_temp(meas_temp = taxon$temp, meas_x = taxon$smr, temp_new = d_roi$temp, E = taxonE['smr'])
d_roi$Pcrit = adj_by_temp(meas_temp = taxon$temp, meas_x = taxon$pcrit_smr, temp_new = d_roi$temp, E = taxonE['pcrit_smr'])
d_roi$alpha = adj_by_temp(meas_temp = taxon$temp, meas_x = taxon$alpha, temp_new = d_roi$temp, E = taxonE['alpha'])


MMR25 = 21 * adj_by_temp(meas_temp = taxon$temp, meas_x = taxon$alpha, temp_new = 25, E = taxonE['alpha'])


d_roi$MMR_temp = adj_by_temp(meas_temp = 25, meas_x = MMR25, temp_new = d_roi$temp, E = E_MMR)

d_roi$MMR_o2 = d_roi$po2 / d_roi$Pcrit * d_roi$SMR
d_roi$FAS_o2 = d_roi$MMR_o2 / d_roi$SMR
d_roi$FAS_temp = d_roi$MMR_temp / d_roi$SMR
d_roi$MI = d_roi$po2 / d_roi$Pcrit
d_roi$MMR_min = pmin(d_roi$MMR_o2, d_roi$MMR_temp)

d_roi$Pcrit_MMR = d_roi$FAS_temp * d_roi$Pcrit

source('run_transects.R')

cowplot::save_plot(filename = paste0('FAS_plots/', species, ' transects.pdf'), plot = cp, base_height = 20, base_width = 24)




p_supplFAS = ggplot(d_transect, aes(lat, depth, height = height)) +
	labs(x = 'Latitude', y = 'Depth (m)') +
	scale_y_reverse(expand = c(0, 0)) +
	scale_x_latitude(ticks = 15) +
	scale_alpha_discrete(range = c(0.6, 1), guide = 'none') +
	geom_tile(aes(fill = FAS_min)) +
	scale_fill_divergent(midpoint = 3) +
	ggtitle(paste0(species, ' FAS (E_MMR = ', E_MMR, ')')) +
	theme(plot.title = element_text(hjust = 0.5),
				legend.title = element_blank())
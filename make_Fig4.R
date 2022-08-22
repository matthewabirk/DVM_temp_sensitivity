source('Make d_roi.R')

# ALTERNATIVELY: d_roi = readRDS('WOA/d_roi.rds')

# Add physiology data -----------------------------------------------------

dg = data.frame(
	temp = c(10L, 20L, 25L),
	smr = c(8.053662, 10.25983, 16.98664),
	mmr = c(17.35, 23.29, 69.77),
	pcrit_smr = c(1.62, 3.26, 4.63)
)

dg$alpha = dg$smr / dg$pcrit_smr


dgE = c(smr = calc_E(dg$smr, dg$temp), mmr = calc_E(dg$mmr, dg$temp), pcrit_smr = calc_E(dg$pcrit_smr, dg$temp), alpha = calc_E(dg$alpha, dg$temp))
dgE = c(dgE, 'pcrit_mmr' = unname(dgE['mmr'] - dgE['alpha']))


dg$pcrit_mmr = dg$mmr / dg$alpha



# Changing E --------------------------------------------------------------


t1_t2_E = function(x1, t1, t2, E){
	kb = 8.61733e-5
	x2 = exp(-E * (1 / (kb * (t2 + 273.15)) - 1 / (kb * (t1 + 273.15))) + log(x1))
	return(x2)
}

E_low_MMR = with(subset(dg, temp <= 20), calc_E(mmr, temp))
E_high_MMR = with(subset(dg, temp >= 20), calc_E(mmr, temp))
dgMMRE = function(temp){
	t1_t2_E(x1 = dg[dg$temp == 20, 'mmr'], t1 = dg[dg$temp == 20, 'temp'], t2 = temp, E = ifelse(temp < 20, E_low_MMR, E_high_MMR))
}

E_low_SMR = with(subset(dg, temp <= 20), calc_E(smr, temp))
E_high_SMR = with(subset(dg, temp >= 20), calc_E(smr, temp))
dgSMRE = function(temp){
	t1_t2_E(x1 = dg[dg$temp == 20, 'smr'], t1 = dg[dg$temp == 20, 'temp'], t2 = temp, E = ifelse(temp < 20, E_low_SMR, E_high_SMR))
}

library(respirometry)
d_roi$SMR = dgSMRE(d_roi$temp)
d_roi$Pcrit = adj_by_temp(meas_temp = dg$temp, meas_x = dg$pcrit_smr, temp_new = d_roi$temp, E = dgE['pcrit_smr'])
d_roi$alpha = adj_by_temp(meas_temp = dg$temp, meas_x = dg$alpha, temp_new = d_roi$temp, E = dgE['alpha'])

d_roi$MMR_temp = dgMMRE(d_roi$temp)

d_roi$MMR_o2 = d_roi$po2 / d_roi$Pcrit * d_roi$SMR
d_roi$FAS_o2 = d_roi$MMR_o2 / d_roi$SMR
d_roi$FAS_temp = d_roi$MMR_temp / d_roi$SMR
d_roi$MI = d_roi$po2 / d_roi$Pcrit
d_roi$MMR_min = pmin(d_roi$MMR_o2, d_roi$MMR_temp)

d_roi$Pcrit_MMR = d_roi$FAS_temp * d_roi$Pcrit

source('run_transects.R')


cowplot::save_plot(filename = 'dg min transects.pdf', plot = cp, base_height = 20, base_width = 24)


# Fig 4 -------------------------------------------------------------------

panel_a = p_FAS_min +
	ggtitle('Dosidicus gigas FAS') +
	theme(plot.title = element_text(hjust = 0.5),
				legend.title = element_blank())

panel_b = p_temp +
	scale_fill_viridis_c(breaks = c(ceiling(min(d_transect$temp)), floor(max(d_transect$temp)), 10, 15, 20, 25)) +
	ggtitle('Temperature (°C)') +
	theme(plot.title = element_text(hjust = 0.5),
				legend.title = element_blank(),
				axis.text.x = element_blank())

panel_c = p_po2 +
	scale_fill_viridis_c(breaks = c(ceiling(min(d_transect$po2)), floor(max(d_transect$po2)), 5, 10, 15, 20)) +
	ggtitle('Oxygen Pressure (kPa)') +
	theme(plot.title = element_text(hjust = 0.5),
				legend.title = element_blank())





species = 'Nyctiphanes simplex'

# E_MMR = 0.3
# source('all_species.R')
# panel_e = p_FAS_min +
# 	scale_fill_divergent(midpoint = 3, limits = c(0.5, 5)) +
# 	ggtitle('Nyctiphanes simplex FAS (E_MMR = 0.3)') +
# 	theme(axis.title.x = element_blank(),
# 				axis.text.x = element_blank(),
# 				plot.title = element_text(hjust = 0.5),
# 				legend.title = element_blank())

E_MMR = 1
source('all_species.R')
panel_f = p_FAS_min +
	scale_fill_divergent(midpoint = 3, limits = c(0.5, 5)) +
	ggtitle('Nyctiphanes simplex FAS (E_MMR = 1.0)') +
	theme(plot.title = element_text(hjust = 0.5),
				legend.title = element_blank())

library(cowplot)

fig4 = plot_grid(panel_b, panel_c, panel_a, panel_f, ncol = 1, labels = 'AUTO')

ggsave('fig4.pdf', fig4, height = 7, width = 4.5)

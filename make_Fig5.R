source('/Users/matthewbirk/Documents/Career/Writing/Manuscripts_in_prep/Dosidicus/final_scripts/import_d_roi.R')

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













# 97/98 ENSO event --------------------------------------------------------

library(ncdf4)

make_enso = function(wd){
	setwd(wd)
	
	enso = data.frame()
	
	# data from https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/
	
	for(i in list.files(pattern = '.nc$')){
		sst = nc_open(i)
		tmp1 = list()
		tmp1$lon = ncvar_get(sst,'lon') # extract the longitude variable
		tmp1$lat = ncvar_get(sst,'lat')
		tmp1$temp = ncvar_get(sst,'sst')
		nc_close(sst)
		sst_data = reshape2::melt(tmp1$temp)
		sst_data = cbind(expand.grid(tmp1$lon, tmp1$lat), sst_data)
		sst_data = sst_data[, -(3:4)]
		colnames(sst_data) = c('lon', 'lat', 'sst')
		sst_data[which(sst_data$lon > 180), 'lon'] = sst_data[which(sst_data$lon > 180), 'lon'] - 360
		sst_data = filter(sst_data, lon < -60, lat < 55, lat > -55)
		
		sst_data$km_per_deg_lon = cos(measurements::conv_unit(sst_data$lat, 'degree', 'radian')) * 111.321 # add column with conversion from 1* lon to km
		
		sst_data$boundary = ifelse(abs(sst_data$lat) > 5, sst_data$lon > (sst_data$lat + 276) / -2.4, sst_data$lon > -85) # ignore Galapagos
		
		coast = sst_data[is.na(sst_data$sst) & sst_data$boundary, ] %>% group_by(lat) %>% summarise(coast = min(lon)) # choose the westernmost NA value that is east of the boundary
		
		sst_data = left_join(sst_data, coast, by = 'lat')
		sst_data$roi_east = sst_data$coast - 500 / sst_data$km_per_deg_lon # should be 500 km
		sst_data$roi_west = sst_data$coast - 2000 / sst_data$km_per_deg_lon # should be 2000 km
		
		sst_roi = filter(sst_data, lon > roi_west, lon < roi_east, !is.na(sst)) %>% select(lon, lat, sst)
		
		sst_roi$date = gsub('oisst-avhrr-v02r01.(\\d+).nc', '\\1', i)
		
		enso = rbind(enso, sst_roi)
	}
	
	enso$date = lubridate::ymd(enso$date)
	
	enso_transect = enso %>% group_by(lat, date) %>% summarise_all(mean, na.rm = TRUE)
	enso_transect$trust = enso_transect$sst >= 10 & enso_transect$sst <= 25
	return(enso_transect)
}






#############################



enso9798 = make_enso('/Users/matthewbirk/Documents/Career/Writing/Manuscripts_in_prep/Dosidicus/97-98 ENSO event')

tmp = enso9798 %>% group_by(date) %>% slice(cumsum(rle(.data$trust)$lengths))
tmp$FAS_temp = NA
tmp = bind_rows(enso9798, tmp)

tmp$PcSMR = adj_by_temp(meas_temp = dg$temp, meas_x = dg$pcrit_smr, temp_new = tmp$sst, E = dgE['pcrit_smr'])
tmp$PcMMR_byfns = dgMMRE(tmp$sst) / dgSMRE(tmp$sst) * tmp$PcSMR
tmp$PcMMR_byfns = pmin(tmp$PcMMR_byfns, 21)
tmp$FAS_byPcSMR = tmp$PcMMR_byfns / tmp$PcSMR
tmp[cumsum(rle(tmp$trust)$lengths), 'FAS_byPcSMR'] = NA # make NA to break lines







d_transect = d_roi %>% group_by(lat, depth) %>% summarise_all(mean, na.rm = TRUE)
inter2 = data.frame(depth = unique(d_roi$depth))
inter2$height = c(5, diff(unique(d_roi$depth)))
inter2[22, 'height'] = 50 # fill in gap between 100 and 125 m depth
d_transect = left_join(d_transect, inter2, by = 'depth')
d_transect = filter(d_transect, depth <= 500)
d_transect$trust = d_transect$temp >= 10 & d_transect$temp <= 25
d_surface = filter(d_transect, depth == 0)
d_surface$PcSMR = adj_by_temp(meas_temp = dg$temp, meas_x = dg$pcrit_smr, temp_new = d_surface$temp, E = dgE['pcrit_smr'])
d_surface$PcMMR_byfns = dgMMRE(d_surface$temp) / dgSMRE(d_surface$temp) * d_surface$PcSMR
d_surface$PcMMR_byfns = pmin(d_surface$PcMMR_byfns, 21)
d_surface$FAS_byPcSMR = d_surface$PcMMR_byfns / d_surface$PcSMR
d_surface[cumsum(rle(d_surface$trust)$lengths), 'FAS_byPcSMR'] = NA # make NA to break lines

d_cmip = data.frame()

for(i in list.files('../CMIP6/', pattern = '.rds', full.names = TRUE)){
	
	cmip = readRDS(i)
	
	cmip_transect = cmip %>% group_by(lat) %>% summarise_all(mean, na.rm = TRUE)
	cmip_transect$trust = cmip_transect$temp >= 10 & cmip_transect$temp <= 25
	cmip_transect$PcSMR = adj_by_temp(meas_temp = dg$temp, meas_x = dg$pcrit_smr, temp_new = cmip_transect$temp, E = dgE['pcrit_smr'])
	cmip_transect$PcMMR_byfns = dgMMRE(cmip_transect$temp) / dgSMRE(cmip_transect$temp) * cmip_transect$PcSMR
	cmip_transect$PcMMR_byfns = pmin(cmip_transect$PcMMR_byfns, 21)
	cmip_transect$FAS_byPcSMR = cmip_transect$PcMMR_byfns / cmip_transect$PcSMR
	cmip_transect[cumsum(rle(cmip_transect$trust)$lengths), 'FAS_byPcSMR'] = NA # make NA to break lines
	cmip_transect$model = gsub('cmip_roi_(.*) SSP5-8.5\\.rds', '\\1', basename(i))
	
	d_cmip = rbind(d_cmip, cmip_transect)
	
}


pan_a = ggplot() +
	scale_y_continuous(breaks = seq(1, 6, by = 0.5)) +
	scale_x_latitude(ticks = 15) +
	scale_color_date(date_labels = '%b %Y', date_breaks = '3 months', low = 'green', high = 'orange') +
	labs(y = 'Dosidicus gigas FAS', color = 'Month during 1997-1998   \nEl Niño') +
	geom_line(data = tmp, aes(lat, FAS_byPcSMR, group = interaction(date, !trust), color = date, linetype = !trust)) +
	geom_line(data = d_surface, aes(lat, FAS_byPcSMR, linetype = !trust), lwd = 0.75) +
	theme_bw() +
	scale_linetype_discrete(guide = 'none')

pan_b = ggplot() +
	scale_y_continuous(breaks = seq(1, 6, by = 0.5)) +
	scale_x_latitude(ticks = 15) +
	labs(y = 'Dosidicus gigas FAS', color = 'CMIP6 SSP5-8.5') +
	geom_line(data = d_surface, aes(lat, FAS_byPcSMR, linetype = !trust), lwd = 0.75) +
	geom_line(data = d_cmip, aes(lat, FAS_byPcSMR, color = model, linetype = !trust), lwd = 0.75) +
	theme_bw() +
	scale_linetype_discrete(guide = 'none')


cowplot::plot_grid(pan_a, pan_b, ncol = 1)
ggsave('../dg_surface_FAS.pdf', height = 4, width = 6)

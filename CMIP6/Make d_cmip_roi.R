########## Create profile data frames from netcdf files #########################

setwd('/Users/matthewbirk/Documents/Career/Writing/Manuscripts_in_prep/Dosidicus/CMIP6')


library(ncdf4);library(lubridate);library(marelac);library(measurements); library(dplyr); library(ggplot2)

for(i in list.files(pattern = '.nc')){
	
	model_id = gsub('CMIP6 - Sea Surface Temperature \\(SST\\) deg C - (.*) - Annual \\(27 models\\)\\.nc', '\\1', i)
	
	cmip_sst = nc_open(i)
	inter1 = list()
	inter1$lat = ncvar_get(cmip_sst, 'lat')
	inter1$lon = ncvar_get(cmip_sst, 'lon')
	inter1$temp = ncvar_get(cmip_sst, 'tos')
	nc_close(cmip_sst)
	
	temp = expand.grid(lon = as.vector(inter1$lon), lat = as.vector(inter1$lat))
	temp$temp = as.vector(inter1$temp)
	attr(temp, 'out.attrs') = NULL
	
	temp = filter(temp, lon < -60)
	temp = filter(temp, lat < 55)
	temp = filter(temp, lat > -55)
	
	temp$lon0.5 = trunc(temp$lon) + ifelse(temp$lon < 0, -0.5, 0.5)
	temp$lat0.5 = trunc(temp$lat) + ifelse(temp$lat < 0, -0.5, 0.5)
	
	d = temp
	
	# add column with conversion from 1* lon to km
	d$km_per_deg_lon = cos(conv_unit(d$lat, 'degree', 'radian')) * 111.321
	
	
	d$boundary = ifelse(abs(d$lat) > 5, d$lon > (d$lat + 276) / -2.4, d$lon > -85) # ignore Galapagos
	
	coast = d[is.na(d$temp) & d$boundary, ] %>% group_by(lat) %>% summarise(coast = min(lon)) # choose the westernmost NA value that is east of the boundary
	
	
	d = left_join(d, coast, by = 'lat')
	d$roi_east = d$coast - 500 / d$km_per_deg_lon
	d$roi_west = d$coast - 2000 / d$km_per_deg_lon
	
	d_roi = filter(d, lon > roi_west, lon < roi_east, !is.na(temp))
	
	library(metR)
	
	ggplot() + 
		geom_raster(data = filter(d), aes(lon, lat, fill = temp)) +
		geom_raster(data = filter(d_roi), aes(lon, lat), fill = 'yellow') +
		#	geom_line(data = data.frame(x = -140:-90, y = -140:-90 * -2.4 - 276), aes(x, y), col = 'red') +
		#	geom_path(data = coast, aes(x = coast, y = lat)) +
		coord_cartesian() +
		scale_x_longitude() +
		scale_y_latitude(breaks = seq(-50, 50, by = 10)) +
		scale_fill_continuous(breaks = seq(0, 30, by = 5)) +
		labs(fill = 'SST (°C)')
	
	ggsave(paste0('cmip_map ', model_id, '.pdf'), device = 'pdf', width = 4, height = 4)
	
	
	
	saveRDS(d_roi, paste0('cmip_roi_', model_id, '.rds'))
	
}

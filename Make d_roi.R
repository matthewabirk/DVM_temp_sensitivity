########## Create profile data frames from netcdf files #########################

setwd('WOA')


library(ncdf4);library(lubridate);library(marelac);library(birk); library(dplyr); library(respirometry)

temp0.25 = nc_open('woa18_decav_t00_04.nc')
inter1 = list()
inter1$lat = ncvar_get(temp0.25, 'lat')
inter1$lon = ncvar_get(temp0.25, 'lon')
inter1$depth = ncvar_get(temp0.25, 'depth')
inter1$temp = ncvar_get(temp0.25, 't_an')
nc_close(temp0.25)

temp = expand.grid(lon = as.vector(inter1$lon), lat = as.vector(inter1$lat), depth = as.vector(inter1$depth))
temp$temp = as.vector(inter1$temp)
attr(temp, 'out.attrs') = NULL

temp = filter(temp, lon < -60)
temp = filter(temp, lat < 55)
temp = filter(temp, lat > -55)

temp$lon0.5 = trunc(temp$lon) + ifelse(temp$lon < 0, -0.5, 0.5)
temp$lat0.5 = trunc(temp$lat) + ifelse(temp$lat < 0, -0.5, 0.5)


oxygen1 = nc_open('woa18_all_O00_01.nc')
inter1 = list()
inter1$lat = ncvar_get(oxygen1, 'lat')
inter1$lon = ncvar_get(oxygen1, 'lon')
inter1$depth = ncvar_get(oxygen1, 'depth')
inter1$o2 = ncvar_get(oxygen1, 'O_an')
nc_close(oxygen1)

o2 = expand.grid(lon = as.vector(inter1$lon), lat = as.vector(inter1$lat), depth = as.vector(inter1$depth))
o2$o2 = as.vector(inter1$o2)
attr(o2, 'out.attrs') = NULL

d = left_join(temp, o2, by = c('lon0.5' = 'lon', 'lat0.5' = 'lat', 'depth'))

# add column with conversion from 1* lon to km
d$km_per_deg_lon = cos(conv_unit(d$lat, 'degree', 'radian')) * 111.321


d$boundary = ifelse(abs(d$lat) > 5, d$lon > (d$lat + 276) / -2.4, d$lon > -85) # ignore Galapagos

coast = d[is.na(d$temp) & d$boundary & d$depth == 0, ] %>% group_by(lat) %>% summarise(coast = min(lon)) # choose the westernmost NA value that is east of the boundary


d = left_join(d, coast, by = 'lat')
d$roi_east = d$coast - 500 / d$km_per_deg_lon
d$roi_west = d$coast - 2000 / d$km_per_deg_lon

d_roi = filter(d, lon > roi_west, lon < roi_east, depth <= 1000, !is.na(temp), !is.na(o2))
d_roi$po2 = respirometry::conv_o2(d_roi$o2, 'percent_a.s.', 'kPa', temp = d_roi$temp)

library(metR)

ggplot() + 
	geom_raster(data = filter(d, depth == 0), aes(lon, lat, fill = temp)) +
	geom_raster(data = filter(d_roi, depth == 0), aes(lon, lat), fill = 'yellow') +
	#	geom_line(data = data.frame(x = -140:-90, y = -140:-90 * -2.4 - 276), aes(x, y), col = 'red') +
	#	geom_path(data = coast, aes(x = coast, y = lat)) +
	coord_cartesian() +
	scale_x_longitude() +
	scale_y_latitude(breaks = seq(-50, 50, by = 10)) +
	scale_fill_continuous(breaks = seq(0, 30, by = 5)) +
	labs(fill = 'SST (°C)')

ggsave('map.pdf', device = 'pdf', width = 4, height = 4)



saveRDS(d_roi, 'd_roi.rds')

setwd('..')

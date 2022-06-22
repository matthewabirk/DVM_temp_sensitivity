# Transects ---------------------------------------------------------------


d_transect = d_roi %>% group_by(lat, depth) %>% summarise_all(mean, na.rm = TRUE)

FAS_o2_threshold = 3
d_transect$habitable = d_transect$FAS_o2 >= FAS_o2_threshold

d_transect[which(d_transect$po2 < d_transect$Pcrit), 'SMR'] = NA

d_transect$FAS_min = d_transect$MMR_min / d_transect$SMR


inter2 = data.frame(depth = unique(d_roi$depth))
inter2$height = c(5, diff(unique(d_roi$depth)))
inter2[22, 'height'] = 50 # fill in gap between 100 and 125 m depth
d_transect = left_join(d_transect, inter2, by = 'depth')

d_transect = filter(d_transect, depth <= 500)
d_transect$trust = d_transect$temp >= 10 & d_transect$temp <= 25
d_transect[which(d_transect$po2 < d_transect$Pcrit), 'trust'] = TRUE

library(birk)

p = ggplot(d_transect, aes(lat, depth, height = height)) +
	labs(x = 'Latitude', y = 'Depth (m)') +
	scale_y_reverse(expand = c(0, 0)) +
	scale_x_latitude(ticks = 15) +
	scale_alpha_discrete(range = c(0.6, 1), guide = 'none')

p_temp = p +
	geom_tile(aes(fill = temp)) +
	scale_fill_gradientn(colours = rev(rainbow(7)), breaks = round(range_seq(d_transect$temp, length.out = 5), 1) + c(0.1, 0, 0, 0, -0.1))


p_po2 = p +
	geom_tile(aes(fill = po2)) +
	scale_fill_gradientn(colours = rev(rainbow(7)), breaks = round(range_seq(d_transect$po2, length.out = 5), 1) + c(0.1, 0, 0, 0, -0.1))



p_Pcrit = p +
	geom_tile(aes(fill = Pcrit, alpha = trust)) +
	scale_fill_gradientn(colours = rev(rainbow(7)), breaks = round(range_seq(d_transect$Pcrit, length.out = 5), 1) + c(0.1, 0, 0, 0, -0.1))




p_Pcrit_MMR = p +
	geom_tile(aes(fill = Pcrit_MMR, alpha = trust)) +
	scale_fill_gradientn(colours = rev(rainbow(7)), limits = c(0, 21), breaks = seq(0, 21, by = 3))


p_SMR = p +
	geom_tile(aes(fill = SMR, alpha = trust)) +
	scale_fill_gradientn(colours = rev(rainbow(7)), breaks = round(range_seq(d_transect$SMR, length.out = 5), 1) + c(0.1, 0, 0, 0, -0.1))


p_MMR_min = p +
	geom_tile(aes(fill = MMR_min, alpha = trust)) +
	scale_fill_gradientn(colours = rev(rainbow(7)))


p_alpha = p +
	geom_tile(aes(fill = alpha, alpha = trust)) +
	scale_fill_gradientn(colours = rev(rainbow(7)))

p_FAS_min = p +
	geom_tile(aes(fill = FAS_min)) +
	scale_fill_divergent(midpoint = 3)

cp = cowplot::plot_grid(p_temp,
												p_po2,
												p_Pcrit,
												p_SMR,
												p_Pcrit_MMR,
												p_MMR_min,
												p_FAS_min,
												p_alpha,
												labels = 'AUTO', nrow = 4, ncol = 2, align = 'hv')

Required R packages:
	birk
	cowplot
	dplyr
	ggplot2
	lubridate
	marelac
	measurements
	metR
	ncdf4
	reshape2
	respirometry

Before anything, download the WOA temperature data that were too big to push to GitHub:
	https://www.ncei.noaa.gov/thredds-ocean/catalog/ncei/woa/temperature/decav/0.25/catalog.html
	Place in the WOA directory.
Also, gunzip WOA/woa18_all_O00_01.nc.gz.

To generate figure 4:
	Run "make_Fig4.R"
	
To generate figure 5:
	Run "make_Fig5.R"
	
To generate figures S4 and S5:
	Run "simulations.R"
	
To generate figure S6:
	Run "all_species.R"
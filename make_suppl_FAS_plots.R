species_options = c('Nyctiphanes simplex', 'Euphausia mucronata', 'Tesserobrachion occulatum', 'Nematobrachion flexipes')
E_MMR_options = c(0.3, 1)
tmp = expand.grid(species_options, E_MMR_options)

apply(tmp, 1, function(i){
	species <<- i['Var1']
	E_MMR <<- as.numeric(i['Var2'])
	source('all_species.R')
})

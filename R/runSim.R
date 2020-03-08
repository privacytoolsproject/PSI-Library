source('functions-unbiasedprivacy.R')
source('utilities-unbiasedprivacy.R')

param_combos <- read.csv('ellen_param_combos.csv')

param_row <- param_combos[5,]

param_row

# run UDP sim
save_path = '.'
udpSim(param_row, save_path = save_path)


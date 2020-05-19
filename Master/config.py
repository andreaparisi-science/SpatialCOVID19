simList = {}
simList['Default'] = {
	'cores': 12,			# Number of cores used by replica
	'gridres': 5,			# Grid resolution in km
	'mobility': 0.2,		# Value for the mobility parameter
	'simLen': 365,			# Duration of simulations for multi-runs
	'model': 				# Model type; options include 'seir', 'full', 'consensus'
		'consensus',
	't0': 21,				# Number of days before Day 0: this is the day for the start of the fitting and intial programmed infections
	'R0': 2.5,				# Value of the basic reproductive number
	'sigma': 1/3.,			# Rate of transition from latent to infective classes
	'gamma': 1/4.,			# Rate of removal
	'tauleap': 'no',		# Tau leaping? NOT OPERATIVE, LEAVE TO 'no'
	'households': 'yes', 	# WHETHER TO INCLUDE HOUSEHOLDS OR NOT
	'betaMul': 1.0,			# CORRECTIVE TERM TO BETA
	'familyAttackMul': 1.0,  # CORRECTIVE TERM TO BETA IN HOUSEHOLD
	'nParticles': 140,  	# Number of particles for fitting procedure
	'fitRestart': 1,			# Whether fitting should restart from scratch
}

simList['Kenya'] = {'cores': 12, 'gridres': 5, 'mobility': 0.870, 'simLen': 365, 'R0': 2.5, 'model': 'consensus', 't0': 21, 'nParticles': 6, 'households': 'yes', 'betaMul': 1.172}
simList['Italy'] = {'cores': 12, 'gridres': 5, 'mobility': 0.395, 'simLen': 60, 'R0': 2.5, 'model': 'consensus', 't0':  21, 'nParticles': 16, 'households': 'yes', 'betaMul': 1.2368}
#simList['Italy'] = {'cores': 12, 'gridres': 5, 'mobility': 0.395, 'simLen': 60, 'R0': 2.5, 'model': 'seir', 't0':  21, 'nParticles': 16, 'households': 'no', 'tauleap': 'no'}
simList['Test'] = {'cores': 1, 'gridres': 5000, 'mobility': 0.870, 'simLen': 365, 'R0': 2.5, 'model': 'consensus', 't0':  0, 'nParticles': 16}


simList = {}
simList['Default'] = {
	# General parameters
	'cores': 12,			# Number of cores used by replica
	'gridres': 5,			# Grid resolution in km
	'mobility': 0.2,		# Value for the mobility parameter
	'simLen': 365,			# Duration of simulations for multi-runs
	'model': 				# Model type; options include 'seir', 'full', 'consensus'
		'consensus',
	'tauleap': 'no',		# Tau leaping? NOT OPERATIVE, LEAVE TO 'no'
	't0': 21,				# Number of days before Day 0: this is the day for the start of the fitting and intial programmed infections

	# Epidemiological parameters
	'R0': 2.5,				# Value of the basic reproductive number
	'sigma': 1/3.,			# Rate of transition from latent to infective classes
	'gamma': 1/4.,			# Rate of removal
	'omega': 0,				# Rate of imports per day
	'households': 'yes', 	# WHETHER TO INCLUDE HOUSEHOLDS OR NOT
	'betaMul': 1.0,			# CORRECTIVE TERM TO BETA

	# Other model parameters
	'rates': 'traditional', # Possible values: 'spanish', 'traditional'
	'stayQuadratic': 'no', # Whether work reduction means stopping people and reducing interaction for the remaining
	'youngestLimit': 19,    # All indivs with age <= this will be affected by young stay-at-home
	'eldestLimit': 70,      # All indivs with age >= this will be affected by old stay-at-home
	'minAgeFirstcase': 30,  # Minimal age to consider for initial cases
	'weeklyWorkHours': 40,  # Number of working hours per week

	# Household parameters
	'avgFamilySize': 3.9,   # Average family size
	'familyAttackMul': 1.0, # CORRECTIVE TERM TO BETA IN HOUSEHOLD

	# Fitting parameters
	'nParticles': 140,  	# Number of particles for fitting procedure
	'fitRestart': 1,		# Whether fitting should restart from scratch (1), continue fitting (0), reset fitting from last (2), forecasting (3); (reserved 4)
	'fitRate': 0.50,        # Fraction of discarded particles
	'fitThreshold': 3.00,   # Target acceptance threshold (in %)
	'allowZeroFits': 'yes', # Whether runs that have no cases should be included or automatically rejected during fitting
}

#simList['Kenya'] = {'mobility': 0.870*0.337, 'simLen': 365, 'R0': 2.5, 'model': 'consensus', 't0': 0, 'nParticles': 16, 'households': 'yes', 'betaMul': 1.172, 'omega': 0}
simList['Kenya'] = {'mobility': 0.870*0.337, 'simLen': 365, 'R0': 2.5, 'model': 'consensus', 't0': 0, 'nParticles':  6, 'households': 'no',  'betaMul': 1.0, 
					'omega': 0, 'avgFamilySize': 3.9, 'weeklyWorkHours' : 45}
simList['Italy'] = {'mobility': 0.395, 'simLen': 365, 'R0': 2.5, 'model': 'consensus', 't0': 21, 'nParticles': 16, 'households': 'yes', 'betaMul': 1.23852, 
					'avgFamilySize': 2.31, 'weeklyWorkHours': 37.5}
#simList['Italy'] = {'mobility': 0.395, 'simLen': 365, 'R0': 2.5, 'model': 'consensus', 't0':  21, 'nParticles': 4, 'households': 'no', 'tauleap': 'no', 'betaMul': 1.00, 'omega': 0.0}
#simList['Italy'] = {'mobility': 0.395, 'simLen': 60, 'R0': 2.5, 'model': 'seir', 't0':  21, 'nParticles': 16, 'households': 'no', 'tauleap': 'no'}
simList['Test'] = {'cores': 12, 'gridres': 20, 'mobility': 0.0750*0.337, 'simLen': 365, 'R0': 2.5, 'model': 'consensus', 't0':  0, 'nParticles': 16, 'households': 'no'}


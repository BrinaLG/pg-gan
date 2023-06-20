#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 15:49:46 2023

@author: blopezgf
"""

import msprime 

T_OOA = 920
demography = msprime.Demography()
demography.add_population(name="AFR", description="African", initial_size=14474)
demography.add_population(
	name="EUR",
	description="European",
	initial_size=34039,
	growth_rate=0.0038,
)
demography.add_population(
	name="EAS",
	description="East Asian",
	initial_size=45852,
	growth_rate=0.0048,
)
demography.add_population(
	name="ADMIX",
	description="Admixed America",
	initial_size=54664,
	growth_rate=0.05,
)
demography.add_population(
	name="OOA",
	description="Bottleneck out-of-Africa",
	initial_size=1861,
)
demography.add_population(
	name="AMH", description="Anatomically modern humans", initial_size=14474
)
demography.add_population(
	name="ANC",
	description="Ancestral equilibrium",
	initial_size=7310,
)
demography.set_symmetric_migration_rate(["AFR", "EUR"], 2.5e-5)
demography.set_symmetric_migration_rate(["AFR", "EAS"], 0.78e-5)
demography.set_symmetric_migration_rate(["EUR", "EAS"], 3.11e-5)

demography.add_admixture(
	12,
	derived="ADMIX",
	ancestral=["AFR", "EUR", "EAS"],
	proportions=[1 / 6, 2 / 6, 3 / 6],
);

demography.add_population_split(T_OOA, derived=["EUR", "EAS"], ancestral="OOA")
demography.add_symmetric_migration_rate_change(
	time=T_OOA, populations=["AFR", "OOA"], rate=15e-5
)
demography.add_population_split(2040, derived=["OOA", "AFR"], ancestral="AMH")
demography.add_population_split(5920, derived=["AMH"], ancestral="ANC")


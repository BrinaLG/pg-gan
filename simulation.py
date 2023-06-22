"""
Simulate data for training or testing using msprime.
Author: Sara Matheison, Zhanpeng Wang, Jiaping Wang, Rebecca Riley
Date: 9/27/22
"""

# python imports
import math
import msprime
import numpy as np
# from stdpopsim
import sps.engines
import sps.species
import sps.HomSap

# our imports
import global_vars
import util

################################################################################
# SIMULATION
################################################################################

def simulate_im(params, sample_sizes, seed, reco):
    """Note this is a 2 population model"""
    assert len(sample_sizes) == 2

    # condense params
    N1 = params.N1.value
    N2 = params.N2.value
    T_split = params.T_split.value
    N_anc = params.N_anc.value
    mig = params.mig.value

    population_configurations = [
        msprime.PopulationConfiguration(sample_size=sample_sizes[0],
            initial_size = N1),
        msprime.PopulationConfiguration(sample_size=sample_sizes[1],
            initial_size = N2)]

    # no migration initially
    mig_time = T_split/2

    # directional (pulse)
    if mig >= 0:
        # migration from pop 1 into pop 0 (back in time)
        mig_event = msprime.MassMigration(time = mig_time, source = 1,
            destination = 0, proportion = abs(mig))
    else:
        # migration from pop 0 into pop 1 (back in time)
        mig_event = msprime.MassMigration(time = mig_time, source = 0,
            destination = 1, proportion = abs(mig))

    demographic_events = [
        mig_event,
		# move all in deme 1 to deme 0
		msprime.MassMigration(
			time = T_split, source = 1, destination = 0, proportion = 1.0),
        # change to ancestral size
        msprime.PopulationParametersChange(time=T_split, initial_size=N_anc,
            population_id=0)
	]

    # simulate tree sequence
    ts = msprime.simulate(
		population_configurations = population_configurations,
		demographic_events = demographic_events,
		mutation_rate = params.mut.value,
		length = global_vars.L,
		recombination_rate = reco,
        random_seed = seed)

    return ts

def simulate_ooa2(params, sample_sizes,seed, reco):
    """Note this is a 2 population model"""
    assert len(sample_sizes) == 2

    # condense params
    T1 = params.T1.value
    T2 = params.T2.value
    mig = params.mig.value

    population_configurations = [
        msprime.PopulationConfiguration(sample_size=sample_sizes[0],
            initial_size = params.N3.value), # YRI is first
        msprime.PopulationConfiguration(sample_size=sample_sizes[1],
            initial_size = params.N2.value)] # CEU/CHB is second

    # directional (pulse)
    if mig >= 0:
        # migration from pop 1 into pop 0 (back in time)
        mig_event = msprime.MassMigration(time = T2, source = 1,
            destination = 0, proportion = abs(mig))
    else:
        # migration from pop 0 into pop 1 (back in time)
        mig_event = msprime.MassMigration(time = T2, source = 0,
            destination = 1, proportion = abs(mig))

    demographic_events = [
        mig_event,
        # change size of EUR
        msprime.PopulationParametersChange(time=T2,
            initial_size=params.N1.value, population_id=1),
		# move all in deme 1 to deme 0
		msprime.MassMigration(time = T1, source = 1, destination = 0,
            proportion = 1.0),
        # change to ancestral size
        msprime.PopulationParametersChange(time=T1,
            initial_size=params.N_anc.value, population_id=0)
	]

    ts = msprime.simulate(
		population_configurations = population_configurations,
		demographic_events = demographic_events,
		mutation_rate = params.mut.value,
		length =  global_vars.L,
		recombination_rate = reco,
        random_seed = seed)

    return ts

def simulate_postOOA(params, sample_sizes, seed, reco):
    """Note this is a 2 population model for CEU/CHB split"""
    assert len(sample_sizes) == 2

    # condense params
    T1 = params.T1.value
    T2 = params.T2.value
    mig = params.mig.value
    #m_EU_AS = params.m_EU_AS.value

    population_configurations = [
        msprime.PopulationConfiguration(sample_size=sample_sizes[0],
            initial_size = params.N3.value), # CEU is first
        msprime.PopulationConfiguration(sample_size=sample_sizes[1],
            initial_size = params.N2.value)] # CHB is second

    # symmetric migration
    #migration_matrix=[[0, m_EU_AS],
    #                  [m_EU_AS, 0]]

    # directional (pulse)
    if mig >= 0:
        # migration from pop 1 into pop 0 (back in time)
        mig_event = msprime.MassMigration(time = T2/2, source = 1,
            destination = 0, proportion = abs(mig))
    else:
        # migration from pop 0 into pop 1 (back in time)
        mig_event = msprime.MassMigration(time = T2/2, source = 0,
            destination = 1, proportion = abs(mig))

    demographic_events = [
        mig_event,
		# move all in deme 1 to deme 0
		msprime.MassMigration(time = T2, source = 1, destination = 0,
            proportion = 1.0),
        # set mig rate to zero (need if using migration_matrix)
        #msprime.MigrationRateChange(time=T2, rate=0),
        # ancestral bottleneck
        msprime.PopulationParametersChange(time=T2,
            initial_size=params.N1.value, population_id=0),
        # ancestral size
        msprime.PopulationParametersChange(time=T1,
            initial_size=params.N_anc.value, population_id=0)
	]

    ts = msprime.simulate(
		population_configurations = population_configurations,
		demographic_events = demographic_events,
        #migration_matrix = migration_matrix,
		mutation_rate = params.mut.value,
		length =  global_vars.L,
		recombination_rate = reco,
        random_seed = seed)

    return ts

def simulate_exp(params, sample_sizes, seed, reco):
    """Note this is a 1 population model"""
    assert len(sample_sizes) == 1

    T2 = params.T2.value
    N2 = params.N2.value

    N0 = N2 / math.exp(-params.growth.value * T2)

    demographic_events = [
        msprime.PopulationParametersChange(time=0, initial_size=N0,
            growth_rate=params.growth.value),
        msprime.PopulationParametersChange(time=T2, initial_size=N2,
            growth_rate=0),
		msprime.PopulationParametersChange(time=params.T1.value,
            initial_size=params.N1.value)
	]

    ts = msprime.simulate(sample_size = sum(sample_sizes),
		demographic_events = demographic_events,
		mutation_rate = params.mut.value,
		length =  global_vars.L,
		recombination_rate = reco,
        random_seed = seed)

    return ts

def simulate_admix(params, sample_sizes, seed, reco):
    """Note this is a 4 population model for Admixed Individuals"""
    mu=1.25e-8 # mutation rate per bp
    rho=1e-8 # recombination rate per bp
    nbp = 1e8 # generate 100 Mb
    N0=7310 # initial population size
    Thum=5920 # time (gens) of advent of modern humans 
    Naf=14474 # size of african population
    Tooa=2040 # number of generations back to Out of Africa 
    Nb=1861 # size of out of Africa population
    mafb=1.5e-4 # migration rate Africa and Out-of-Africa 
    Teu=920 # number generations back to Asia-Europe split 
    Neu=1032; Nas=554 # bottleneck population sizes 
    mafeu=2.5e-5; mafas=7.8e-6; meuas=3.11e-5 # mig. rates 
    reu=0.0038 # growth rate per generation in Europe 
    ras=0.0048 # growth rate per generation in Asia 
    Tadmix=12 # time of admixture
    Nadmix=30000 # initial size of admixed population 
    radmix=.05 # growth rate of admixed population
    NeADMIX= params.NeADMIX.value
    Tadmix = params.Tadmix.value
    propAFR = params.propJOINT.value[0]
    propEUR = params.propJOINT.value[1]
    propEAS = 1 - (propAFR + propEUR)
    radmix = params.radmix.value
    # pop0 is Africa, pop1 is Europe, pop2 is Asia, pop3 is admixed
    pop_config = [
    msprime.PopulationConfiguration(sample_size=sample_sizes[0],initial_size=Naf,growth_rate=0.0), msprime.PopulationConfiguration(sample_size=sample_sizes[1],initial_size=Neu*math.exp(reu*Teu),growth_rate=reu), msprime.PopulationConfiguration(sample_size=sample_sizes[2],initial_size=Nas*math.exp(ras*Teu),growth_rate=ras), msprime.PopulationConfiguration(sample_size=sample_sizes[3],initial_size=Nadmix*math.exp(radmix*Tadmix),growth_rate=radmix)] 
    mig_mat = [[0,mafeu,mafas,0],[mafeu,0,meuas,0], [mafas,meuas,0,0],[0,0,0,0]]
    # Admixture event, 1/6 Africa, 2/6 Europe, 3/6 Asia
    admixture_event = [ msprime.MassMigration(time=Tadmix,source=3,destination=0,proportion=1.0/6.0), msprime.MassMigration(time=Tadmix+0.0001,source=3,destination=1,proportion=2.0/5.0), msprime.MassMigration(time=Tadmix+0.0002,source=3,destination=2,proportion=1.0)]

    # Asia and Europe split
    eu_event = [
    msprime.MigrationRateChange(time=Teu,rate=0.0), msprime.MassMigration(time=Teu+0.0001,source=2,destination=1,proportion=1.0), msprime.PopulationParametersChange(time=Teu+0.0002,initial_size=Nb,growth_rate=0.0,population_id=1), msprime.MigrationRateChange(time=Teu+0.0003,rate=mafb,matrix_index=(0,1)), msprime.MigrationRateChange(time=Teu+0.0003,rate=mafb,matrix_index=(1,0))]
    # Out of Africa event
    ooa_event = [
    msprime.MigrationRateChange(time=Tooa,rate=0.0), msprime.MassMigration(time=Tooa+0.0001,source=1,destination=0,proportion=1.0)]

    # initial population size
    init_event = [ msprime.PopulationParametersChange(time=Thum,initial_size=N0,population_id=0)]

    events = admixture_event + eu_event + ooa_event + init_event

    ts = msprime.simulate(
    population_configurations=pop_config,
    migration_matrix=mig_mat,
    demographic_events=events,
    length=global_vars.L,
    mutation_rate=mu, # Should be mu!
    recombination_rate=reco)
 
    return ts

def simulate_const(params, sample_sizes, seed, reco):
    assert len(sample_sizes) == 1

    # simulate data
    ts = msprime.simulate(sample_size=sum(sample_sizes), Ne=params.Ne.value,
        length=global_vars.L, mutation_rate=params.mut.value,
        recombination_rate=reco, random_seed = seed)

    return ts

def simulate_ooa3(params, sample_sizes, seed, reco):
    """From OOA3 as implemented in stdpopsim"""
    assert len(sample_sizes) == 3

    sp = sps.species.get_species("HomSap")

    mult = global_vars.L/141213431 # chr9
    contig = sp.get_contig("chr9",length_multiplier=mult) # TODO vary the chrom

    # 14 params
    N_A = params.N_A.value
    N_B = params.N_B.value
    N_AF = params.N_AF.value
    N_EU0 = params.N_EU0.value
    N_AS0 = params.N_AS0.value
    r_EU = params.r_EU.value
    r_AS = params.r_AS.value
    T_AF = params.T_AF.value
    T_B = params.T_B.value
    T_EU_AS = params.T_EU_AS.value
    m_AF_B = params.m_AF_B .value
    m_AF_EU = params.m_AF_EU.value
    m_AF_AS = params.m_AF_AS.value
    m_EU_AS = params.m_EU_AS.value

    model = sps.HomSap.ooa_3(N_A, N_B, N_AF, N_EU0, N_AS0, r_EU, r_AS, T_AF,
        T_B, T_EU_AS, m_AF_B, m_AF_EU, m_AF_AS, m_EU_AS)
    samples = model.get_samples(sample_sizes[0], sample_sizes[1],
        sample_sizes[2]) #['YRI', 'CEU', 'CHB']
    engine = sps.engines.get_engine('msprime')
    ts = engine.simulate(model, contig, samples)

    return ts



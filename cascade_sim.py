from crpropa import *
import numpy as np

N = input("Injection particles: ")


#---------------------------------
# SIMULATION
#---------------------------------

sim = ModuleList()
Emin = 1 * GeV
Emax = 100 * TeV
dmax = 585 * Mpc
    #Propagation
sim.add( SimplePropagation() )
    #Interactions
sim.add( EMInverseComptonScattering( CMB(), True ) )
sim.add( EMInverseComptonScattering( IRB_Franceschini08(), True ) )
sim.add( EMPairProduction( CMB(), True ) )
sim.add( EMPairProduction( IRB_Franceschini08(), True ) )
sim.add( EMDoublePairProduction( CMB(), True ) )
sim.add( EMDoublePairProduction( IRB_Franceschini08(), True ) )
sim.add( EMTripletPairProduction( CMB(), True ) )
sim.add( EMTripletPairProduction( IRB_Franceschini08(), True ) )
    #Limits
sim.add( MaximumTrajectoryLength( dmax ) )
sim.add( MinimumEnergy( Emin ) )

#---------------------------------
# PLASMA INSTABILITIES
#---------------------------------

alpha = 0.0
tildeE = 1 * TeV
lambda_0 = 1.2*kpc
sim.add( PlasmaInstability( alpha, tildeE, lambda_0 ) )

#---------------------------------
# OBSERVERS
#---------------------------------

    #First observer sphere around source (for initial spectrum)
obs1 = Observer()
obs1.add( ObserverSurface( Sphere( Vector3d(0.), 2 * pc ) ) )
obs1.setDeactivateOnDetection(False)
t1 = TextOutput( "blazar_initial_spectrum.txt", Output.Event3D )
obs1.onDetection( t1 )
sim.add( obs1 )

    #Second observer sphere around source (for final spectrum)
obs2 = Observer()
obs2.add( ObserverSurface( Sphere( Vector3d(0.), 580 * Mpc ) ) )
obs2.setDeactivateOnDetection(True)
t2 = TextOutput( "blazar_cascade_spectrum.txt", Output.Event3D )
obs2.onDetection( t2 )
sim.add( obs2 )
print(obs2)

#---------------------------------
# SOURCE
#---------------------------------

source = Source()
source.add( SourcePosition( 0 * Mpc ) )
    #Target emission of photons
angle = 0.17453
P = 0.90
kappa = np.log( 1 - P ) / ( np.cos( angle ) - 1 )
source.add( SourceDirectedEmission( Vector3d(1,0,0), kappa ) )
source.add( SourceParticleType( 22 ) )
    #Spectrum
source.add( SourcePowerLawSpectrum( Emin, Emax, -1 ) )
print(source)

#---------------------------------
# RUNNING
#---------------------------------

sim.setShowProgress(True)
sim.run( source, int(N) )
t1.close()
t2.close()

import numpy as np

#---------------------------------

# INITIAL VALUES
num_bins_E_0 = 50   #Number og initial energy bins
num_bins_E = 50     #Number of energy bins
num_bins_t_0 = 50   #Number of initial time bins
num_bins_t = 50     #Number of time bins

#---------------------------------

# COMPUTATION OF CASCADE SIGNAL G

    #Loading data and initial setting
data = np.loadtxt( "blazar_cascade_spectrum.txt" )
data0 = np.loadtxt( "blazar_initial_spectrum.txt" )
E_0_norm = np.log10( data0[:,10] * ( 10**6 ) )
E = np.log10( data[:,2] * ( 10**6 ) )
E_0 = np.log10( data[:,10] * ( 10**6 ) )
t = ( ( data[:,0] * ( 3.086 * ( 10**22 ) ) ) / ( 3 * ( 10**8 ) ) ) - ( ( 580 * ( 3.086 * ( 10**22 ) ) ) / ( 3 * ( 10**8 ) ) )
index = np.where( np.logical_and( data[:,1] == 22, data[:,1] == 22 ) )
E_0 = E_0[index]
E = E[index]
t = t[index]
N = float(len(E))
E_0_max = 2
E_0_min = -3
E_max = E_0_max
E_min = E_0_min

    #Artificial t_0 dependence function
def gaussian(x, mu, sig):
            return np.exp( -np.power( x - mu, 2. ) / ( 2 * np.power(sig,2.) ) ) * 1000
x_values = np.linspace( -3, 3, num_bins_t_0 )
mu = 0
sig = 1
    #Single-time emission
t_0 = np.ones( num_bins_t_0 )
    #Delay time limits
timestep = 6308001
delta_t = 1
t_min = 0
t_max = 3.154e+8

    #Assigning weights for t_0
weight = np.floor( t_0 )
w = weight.astype(int)
counts_E_0, edges_E_0 = np.histogram( E_0, bins=num_bins_E_0, range=( E_0_min, E_0_max ), density=False )
counts_t, edges_t = np.histogram( t, bins=num_bins_t, range=( t_min, t_max ), density=False )
G = np.zeros((num_bins_t, num_bins_E))
    #Time steps (moving along t_0)
for i in range(0,num_bins_t_0):
    time = t + (i*timestep)
        #Selecting energies at source (moving along E_0)
    for j in range(0,num_bins_E_0):
        index1 = np.where(np.logical_and(E_0>=edges_E_0[j], E_0<=edges_E_0[j+1]))
        EE = E[index1]
        tt = time[index1]
            #Correct normalization
        index_norm = np.where(np.logical_and(E_0_norm>=edges_E_0[j], E_0_norm<=edges_E_0[j+1]))
        E_norm = E_0_norm[index_norm]
        N = float(len(E_norm))
            #Going row by row, sorting by energy
        for k in range(0,num_bins_t):
            index2 = np.where(np.logical_and(tt>=edges_t[k], tt<=edges_t[k+1]))
            EEE = EE[index2]
            a = np.full(len(EEE), w[i])
            counts_E, edges_E = np.histogram(EEE, bins=num_bins_E, range=(E_min, E_max), weights=a, density=False)
            #if w[i] == 0:
                #G[num_bins_t-1-k,:] = 0
            #elif N == 0:
            if N == 0:
        G[num_bins_t-1-k,:] = 0
            else:
                G[num_bins_t-1-k,:] = counts_E/N
            if k == num_bins_t-1:
                mat = np.matrix(G)
                with open('MATRICES/matrix_powerlaw_int_{t_0}_{E_0}.txt'.format(t_0 = np.around(edges_t[i], decimals=3), E_0 = np.around(edges_E_0[j], decimals=3)),'wb') as f:
                    for line in mat:
                        np.savetxt(f, line, fmt='%.2f')
            else:
                pass

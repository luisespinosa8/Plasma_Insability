import numpy as np

# INPUT BLAZAR PARAMETERS

A = 1   #Normalization
gamma = 0.96    #Spectral index
Ecut = 7.7  #Cut-off

# INITIAL VALUES. MUST BE THE SAME AS IN "blazar_eval.py"
num_bins_E_0 = 50
num_bins_t_0 = 50
num_bins_E = num_bins_E_0
num_bins_t = num_bins_t_0
logE0_max = 2
logE0_min = -3

# LOADING DATA AND CONSTRUCTION OF J(logE_0)
data0 = np.genfromtxt("blazar_initial_spectrum.txt")
data = data0
logE0 = np.log10(data0[:,10]*(10**6))
index = np.where(np.logical_and(data0[:,1] == 22,data0[:,1] == 22))
logE0 = logE0[index]
t = ((data[:,0]*(3.086*(10**22)))/(3*(10**8))) - ((580*(3.086*(10**22)))/(3*(10**8)))
t_min = 0
t_max = 3.154e+8
counts_t, edges_t = np.histogram(t, bins=num_bins_t, range=(t_min, t_max), density=False)
counts_E_0, edges_E_0 = np.histogram(logE0, bins=num_bins_E_0, range=(logE0_min,logE0_max), density=False)
index_range = np.where(np.logical_and(counts_E_0 >0, counts_E_0 > 0))
dlogE = edges_E_0[1:]-edges_E_0[:-1]
E = 10**(edges_E_0[0:num_bins_E_0]+(dlogE/2))
J_E0 = np.zeros((num_bins_E_0,num_bins_E_0))
J_E0[num_bins_t-1,index_range] = ((A*(E[index_range]/1)**(-gamma))*np.exp(-E[index_range]/Ecut))
J_logE0 = J_E0*E

# SAVING THE MATRICES IN A DICTIONARY OF A DICTIONARY. THEN G IS FOUR-DIMENSIONAL AND HAS 4 INDICES
G = {-1:{}}
for a in range(0,num_bins_t_0):
    G[a] = {}
    for b in range(0,num_bins_E_0):
        G[a]['{E_0}'.format(E_0=b)] = np.loadtxt('MATRICES/matrix_powerlaw_int_{t_0}_{E_0}.txt'.format(t_0 = np.around(edges_t[a], decimals=3), E_0 = np.around(edges_E_0[b], decimals=3)))

# COMPUTATION OF J(log(E)). FOUR FOR LOOPS FOR THE FOUR DIMENSIONS OF G
J_logE = np.zeros((len(J_logE0),len(J_logE0)))
J_E = np.zeros((len(J_logE0),len(J_logE0)))
for i in range(0,num_bins_E):
    for j in range(0,num_bins_t):
        for k in range(i,num_bins_E_0):
            for w in range(0,num_bins_t_0):            
                J_logE[j,i] += G[w]['{E_0}'.format(E_0=k)][j,i]*J_logE0[num_bins_t-1-w,k]
        J_E[j,i] = J_logE[j,i]/E[i]

np.savetxt('spectrum.txt',J_E[num_bins_t-1,:])

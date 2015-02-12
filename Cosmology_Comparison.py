import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize

plt.rc('font', family='serif') #Changes all plotting fonts.

#Parameters to change
suffix = 'log.png'              #Plot tags
Omega_R = 5.04*10**(-5)      #EQ blah Dodelson (COBE?)

#Constants
H0  = 70/(3.08567758*10**19) #Hubble constant today in s^-1 (70 km/s/Mpc)
gyr = 3.15569*10**16         #Number of seconds in a Gyr

#quantities = ['a','t','H','rho_M','rho_R','rho_DE','rho_tot']


CosmologyDict = {'EdS':(0,1,0,-1),'Closed':(0,2,0,-1),'Open':(0,0.3,0,-1),
                 'LCDM':(0,0.3,0.7,-1),'Quint':(0,0.3,0.7,-2/3.),
                 'Phant':(0,0.3,0.7,-4/3.),
                 'EdSR':(Omega_R,1-Omega_R,0,-1),
                 'ClosedR':(Omega_R,2-Omega_R,0,-1),
                 'OpenR':(Omega_R,0.3-Omega_R,0,-1),
                 'LCDMR':(Omega_R,0.3-Omega_R,0.7,-1),
                 'QuintR':(Omega_R,0.3-Omega_R,0.7,-2/3.),
                 'PhantR':(Omega_R,0.3-Omega_R,0.7,-4/3.)
                }

def cosmo_dict(cosmo):
    '''
    Returns an array in the form (Omega_R,Omega_M,Omega_DE,w).
    '''
    #if cosmo in CosmologyDict.keys():
    try:
        return CosmologyDict[cosmo]
    except KeyError:
        print 'Unrecognized Cosmology, please choose from:'
        print CosmologyDict.keys()   
    return ''

def dadt(a,cosmo):
    Omega_R,Omega_M,Omega_DE,w = cosmo_dict(cosmo)
    Omega_k = Omega_R + Omega_M + Omega_DE
    return H0*np.sqrt( Omega_R/a**2 + Omega_M/a + Omega_DE/a**(1+3*w)+(1-Omega_k) )

def lookback_time(a,cosmo):
    Omega_R,Omega_M,Omega_DE,w = cosmo_dict(cosmo)
    Omega_k = 1. - Omega_R + Omega_M + Omega_DE
    dt = lambda x: 1./dadt(x,cosmo)
    return quad(dt,1,a)[0]

def cosmic_time(a,cosmo):
    t_0 = lookback_time(0,cosmo)
    t_a = lookback_time(a,cosmo)
    return -t_0 + t_a

def Hubble(a,cosmo):
    return dadt(a,cosmo)/a

def density(a,cosmo):
    '''
    Returns [rho_R, rho_M, rho_DE, rho_total] in units of rho_c_0.
    density(a,cosmo)/(H(a,cosmo)/H0)**2 is in units of rho_c(a)
    '''
    Omega_R,Omega_M,Omega_DE,w = cosmo_dict(cosmo)
    rho_R  = Omega_R/a**4
    rho_M  = Omega_M/a**3
    rho_DE = Omega_DE/a**(3+3*w)
#    return np.array([rho_R, rho_M, rho_DE, rho_R+rho_M+rho_DE])
    return rho_R, rho_M, rho_DE, rho_R+rho_M+rho_DE

def check_closure(cosmo):
    Omega_R,Omega_M,Omega_DE,_ = cosmo_dict(cosmo)
    if Omega_R+Omega_M+Omega_DE > 1:
        return True
    else:
        return False

def find_max_size(cosmo):
    if check_closure(cosmo):
        max_a = minimize(lambda x: dadt(x,cosmo),1,method='Nelder-Mead',
                         options={'xtol': 10**-10}
                        ).x[0]
        return max_a
    else:
        print 'Cosmology is not closed.'
        return None

def make_quantities(cosmo,a_start = 0, a_end = 1, npoints =1000000):
    a_sample  = np.linspace(a_start,a_end,npoints)
    t_sample  = np.array([cosmic_time(a,cosmo) for a in a_sample])
    H_sample  = np.array(Hubble(a_sample,cosmo))
    rho_all   = np.array(density(a_sample,cosmo))
    rho_R   = rho_all[0,:]
    rho_M   = rho_all[1,:]
    rho_DE  = rho_all[2,:]
    rho_tot = rho_all[3,:]

    np.savez(cosmo+'a_'+str(a_end),
             a=a_sample,t=t_sample,H=H_sample,
             rho_M=rho_M,
             rho_R=rho_R,
             rho_DE=rho_DE,
             rho_tot=rho_tot)

def load(quantity, cosmo,a_end):
    quantities = ['a','t','H','rho_M','rho_R','rho_DE','rho_tot']
    if quantity in quantities:
        quant = np.load(cosmo+'a_'+str(a_end)+'.npz')
        return quant[quantity][1::]
    else:
        print 'Quantity not recognized. Choose from:'
        print quantities
    
def quad_plot(cosmo, a_start = 0, a_end = 1, npoints =10000,
              axscale='linear', ayscale='linear', 
              rho_x_scale='linear', rho_y_scale='linear',
              f=None, return_axes=False,suffix=suffix
             ):


    Omega_R,Omega_M,Omega_DE,_ = cosmo_dict(cosmo)

    a_sample  = np.linspace(a_start,a_end,npoints)
    t_sample  = np.array([cosmic_time(a,cosmo) for a in a_sample])
    H_sample  = np.array(Hubble(a_sample,cosmo))
    rho_all= np.array(density(a_sample,cosmo))

    rho_R   = rho_all[0,:]
    rho_M   = rho_all[1,:]
    rho_DE  = rho_all[2,:]
    rho_tot = rho_all[3,:]

    rho_all = rho_all/(H_sample/H0)**2

    if f == None:
        fig = plt.figure()
        ax1 = plt.subplot2grid((2,2), (0,0))
        ax2 = plt.subplot2grid((2,2), (0,1))
        ax3 = plt.subplot2grid((2,2), (1,0))
        ax4 = plt.subplot2grid((2,2), (1,1))
    else:
        fig = f[0]
        ax1 = f[1]
        ax2 = f[2]
        ax3 = f[3]
        ax4 = f[4]

    #Plot a(t)
    ax1.plot(t_sample/gyr,a_sample,label=cosmo,color='black')
    #ax1.scatter(0,1,label='Today',color='black')#Today Value
    ax1.set_title('Scale Factor')
    ax1.set_xlabel(r'$t$ [Gyr]')
    ax1.set_ylabel(r'$a$')
#        ax1.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
#        ax1.set_ylim([0,np.max(As)])
    ax1.set_xscale(axscale)
    ax1.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)

    #Plot H(t)
    ax2.plot(t_sample/gyr,np.arcsinh(H_sample*70/H0),label=cosmo,color='black')
    #ax2.scatter(0,np.arcsinh(70),label='Today',color='black')#Today Value
    ax2.set_title('Hubble "Constant"')
    ax2.set_xlabel(r'$t$ [Gyr]')
    ax2.set_ylabel(r'$H$ [km/s/Mpc]')
    newticks = np.sinh(ax2.get_yticks())
    newticks = ["%.1e" % tick for tick in newticks]
    for j in range(len(newticks)):
        [mant,exp] = newticks[j].split('e')
        exp = str(int(exp))
        newticks[j] = r'${!s}\times 10^{!s}$'.format(mant,exp)
    ax2.set_yticklabels(newticks)
#        ax2.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
    ax2.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)

    #Plot rho(t) [rho_c_0]
    rho_R   = rho_all[0,:]
    rho_M   = rho_all[1,:]
    rho_DE  = rho_all[2,:]
    rho_tot = rho_all[3,:]
    ax3.plot(t_sample/gyr, rho_tot, label=cosmo+':Total',color='black')
    #Add Radiation if applicable
    if Omega_R > 0:
        ax3.plot(t_sample/gyr, rho_R, label=cosmo+':Radiation',color='red')
    #Add Matter
    if Omega_M > 0:
        ax3.plot(t_sample/gyr, rho_M, label=cosmo+':Matter',color='blue')
    #Add Dark Energy if applicable
    if Omega_DE> 0:
        ax3.plot(t_sample/gyr, rho_DE,label=cosmo+':Dark Energy',color='green')
    ax3.set_title('Density')
    ax3.set_xlabel(r'$t$ [Gyr]')
    ax3.set_ylabel(r'$\rho/\rho_{c,0}$')
    ax3.set_xscale(rho_x_scale)
    ax3.set_yscale(rho_y_scale)
#        ax3.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
    ax3.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)

    #Plot rho(t) [rho_c(t)]
    rho_all = rho_all/(H_sample/H0)**2 #Change Units
    rho_R   = rho_all[0,:]
    rho_M   = rho_all[1,:]
    rho_DE  = rho_all[2,:]
    rho_tot = rho_all[3,:]
    ax4.plot(t_sample/gyr, rho_tot, label=cosmo+':Total',color='black')
    #Add Radiation if applicable
    if Omega_R > 0:
        ax4.plot(t_sample/gyr, rho_R, label=cosmo+':Radiation',color='red')
    #Add Matter
    if Omega_M > 0:
        ax4.plot(t_sample/gyr, rho_M, label=cosmo+':Matter',color='blue')
    #Add Dark Energy if applicable
    if Omega_DE> 0:
        ax4.plot(t_sample/gyr, rho_DE,label=cosmo+':Dark Energy',color='green')
    ax4.set_title('Density')
    ax4.set_xlabel(r'$t$ [Gyr]')
    ax4.set_ylabel(r'$\rho/\rho_c(t)$')
    ax4.set_xscale(rho_x_scale)
    ax4.set_yscale(rho_y_scale)
#        ax4.set_xlim([np.min(times)/gyr,np.max(times)/gyr])
    ax4.legend(loc='best', scatterpoints=1,prop={'size':8}, framealpha=0)

    plt.tight_layout()
    plt.savefig(cosmo+suffix)
    if return_axes:
        return fig,ax1,ax2,ax3,ax4
    else:
        plt.close()

def plot_all_in_one(y,x='a', a_end=1):
    colors    = ['r','g','b','c','m','k']
    linetypes = ['-','--']
    color_count=0
    for cosmo in CosmologyDict.keys():
        if 'R' in cosmo:
            pass    
        else:
            color = colors[color_count]
            color_count=color_count+1
            a = load(x, cosmo,a_end)
            if 'rho_tot' in y:
                q = load('rho_tot', cosmo,a_end)
            else:
                q = load(y, cosmo,a_end)

            if y=='H':
                q = q*70/H0
            elif y=='rho_tot_t':
                h = np.array(load('H',cosmo,a_end))
                q = q/(h/H0)**2
            plt.plot(a,q,label=cosmo,color=color, ls='-')

            a = load(x, cosmo+'R',a_end)
            if 'rho_tot' in y:
                q = load('rho_tot', cosmo+'R',a_end)
            else:
                q = load(y, cosmo+'R',a_end)
            if y=='H':
                q = q*70/H0
            elif y=='rho_tot_t':
                h = load('H',cosmo+'R',a_end)
                q = q/(h/H0)**2
            plt.plot(a,q,label=cosmo,color=color, ls='--')

    plt.yscale('log')
    plt.xlim([0,10])      #Dynamically Adjust?
#    plt.ylim([10,10000000]) #Dynamically Adjust?
    plt.ylabel(y)# [km/s/Mpc]')
    plt.xlabel(x)
    plt.legend(loc='best')
    plt.show()

def EdS_Analytic_check():
    t = np.linspace(0,9*gyr,10000)
    t_0 = lookback_time(0,'EdS')
    t = t+t_0
    a=lambda t: np.power(3*H0*t/2. +1,2/3.)
    H=lambda t: H0/(3*H0*t/2. +1)
    rho=lambda t: np.power(3*H0*t/2. +1,-2)
    
    t_num =load('t', 'EdS',1)
    a_num = load('a', 'EdS',1)
    H_num = load('H', 'EdS',1)
    rho_num = load('rho_tot', 'EdS',1)
    
    t_num=t_num/gyr#now in Gyrs
    plt.plot(t/gyr-t_0/gyr,a(t), color='black')
    plt.plot(t_num,a_num, color='red', ls='--')
    plt.title('a')
    plt.show()

    plt.plot(t/gyr-t_0/gyr,H(t), color='black')
    plt.plot(t_num,H_num, color='red', ls='--')
    plt.yscale('log')
    plt.title('H')
    plt.show()

    plt.plot(t/gyr-t_0/gyr,rho(t), color='black')
    plt.plot(t_num,rho_num, color='red', ls='--')
    plt.yscale('log')
    plt.title('rho')
    plt.show()
if __name__ == '__main__':

    EdS_Analytic_check()



#    quantities = ['a','t','H','rho_M','rho_R','rho_DE','rho_tot']
    #for cosmo in CosmologyDict.keys():
    #    make_quantities(cosmo,a_start = 0, a_end = 10, npoints =1000)

#    plot_all_in_one('H',a_end=10)
#    plot_all_in_one('rho_M', a_end=10)
#    plot_all_in_one('rho_R', a_end=10)
#    plot_all_in_one('rho_DE', a_end=10)
#    plot_all_in_one('rho_tot', a_end=10)
#    plot_all_in_one('rho_tot_t', a_end=10)

    make_H_a       = False
    make_rho_tot_a = False

    if make_H_a:
        colors    = ['r','g','b','c','m','k']
        linetypes = ['-','--']
        color_count=0
        for cosmo in CosmologyDict.keys():
            if 'R' in cosmo:
                pass    
            else:
                color = colors[color_count]
                color_count=color_count+1
                a = load('a', cosmo,1)
                H = load('H', cosmo,1)*70/H0
                plt.plot(a,H,label=cosmo,color=color, ls='-')

                a = load('a', cosmo+'R',1)
                H = load('H', cosmo+'R',1)*70/H0
                plt.plot(a,H,label=cosmo,color=color, ls='--')
        plt.yscale('log')
        plt.xlim([0,0.01])
        plt.ylim([10,10000000])
        plt.ylabel('H(a) [km/s/Mpc]')
        plt.xlabel('a')
        plt.legend(loc='best')
        plt.show()

    if make_rho_tot_a:
        colors    = ['r','g','b','c','m','k']
        linetypes = ['-','--']
        color_count=0
        for cosmo in CosmologyDict.keys():
            if 'R' in cosmo:
                pass    
            else:
                color = colors[color_count]
                color_count=color_count+1
                a = load('a', cosmo,1)
                rho_tot = load('rho_tot', cosmo,1)
                plt.plot(a,rho_tot,label=cosmo,color=color, ls='-')

                a = load('a', cosmo,1)
                rho_tot = load('rho_tot', cosmo+'R',1)
                plt.plot(a,rho_tot,label=cosmo,color=color, ls='--')
        plt.yscale('log')
        plt.xlim([0,1])
        plt.ylim([0.01,10000])
        plt.ylabel('rho_tot(a) [?]')
        plt.xlabel('a')
        plt.legend(loc='best')
        plt.show()

#        print cosmic_time(1,'Phant')/gyr
#        print cosmic_time(1,'PhantR')/gyr

#Run the below once to generate the data out to your scale factor of interest.
#Then just use load(quantity,cosmo,a_end), ie 1 in the below.
#    for cosmo in CosmologyDict.keys():
#        make_quantities(cosmo,a_start = 0, a_end = 1, npoints =1000)

#        quad_plot(cosmo, a_start = 0, a_end = 1, npoints =10000,
#                  axscale='linear', ayscale='linear', 
#                  rho_x_scale='linear', rho_y_scale='linear',
#                  f=None, return_axes=False,suffix = suffix
#                 )
#        quad_plot(cosmo, a_start = 0, a_end = 1, npoints =10000,
#                  axscale='linear', ayscale='linear', 
#                  rho_x_scale='log', rho_y_scale='log',
#                  f=None, return_axes=False,suffix = 'log'+suffix
#                 )










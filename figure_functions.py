import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


plt.rc('font', size = 15)
plt.rc('legend', fontsize = 12)


#based on Perets & Naoz 09, v1 on arxiv supplement D
#g aka aop aka w aka argument of pericenter
def find_e_max(i0, e0, g0):
    X = (5*e0*e0*np.sin(g0)**2 + 2*(1-e0*e0)) *np.sin(i0)**2
    Y = (1-e0*e0)*np.cos(i0)**2
    
    #solution for w=0
    imax_w0 = np.arctan(np.sqrt(X/Y/2.))
    emin_w0 = np.sqrt(1-Y/np.cos(imax_w0)**2)
    
    #solutions for w=pi/2                            
    A = -3.
    B = (1.-3.*Y+X)
    C = 2.-2.*Y-X
    D = B*B - 4.*A*C

#            e2 = (-B +- np.sqrt(D))/2/A            
    emin_1 = np.sqrt((-B + np.sqrt(D))/2/A)             
    emin_2 = -1*np.sqrt((-B + np.sqrt(D))/2/A)             
    emin_3 = np.sqrt((-B - np.sqrt(D))/2/A)             
    emin_4 = -1*np.sqrt((-B - np.sqrt(D))/2/A)             
    
    emin = np.array([e0, emin_w0, emin_1, emin_2, emin_3, emin_4])
#            print emin
#            emax = max(emin)
    emax = np.nanmax(emin, axis=0)
    imin = np.arccos(np.sqrt((1-e0*e0)/(1-emax*emax) * np.cos(i0)**2))
    return emax





def plot_figure(m1,m2,m3,e_in_max,incl,oct_limit,list_of_objects,PLOT_PERIOD_RATIO, SAVE_FIG=False, fig_name='TRES_diagnostic.pdf'):
    
    #Object plotting stuff
    e_out_objects = []
    a_ratio_objects = []
    errors =[]
    labels = []
    for i in list_of_objects:
        a_in = i[0]
        a_out = i[1]
        e_out = i[2]
        error = i[3]
        label = i[4]
        e_out_objects.append(e_out)
        a_ratio_objects.append(a_out/a_in)
        errors.append(error)
        labels.append(label)
    



    q_out   = m3/(m1+m2) 
    m_in    = m1+m2
    m_tot   = m1+m2+m3
    

    
    lines = np.array(['-', '--', '-.', ':', '--'])
    colors = np.array(['b', 'g', 'r','c', 'm'])
    e_out_vec = np.arange(1000)/1000.

    ar_stable = 2.8/(1-e_out_vec) * ( (1+q_out)*(1+e_out_vec)/np.sqrt(1-e_out_vec) )**0.4 *(1-0.3*incl/np.pi)
    pr_stable = ar_stable**1.5 * np.sqrt(m_in/m_tot)

    ar_oct = (m1-m2)/(m1+m2) *  e_out_vec/(1-e_out_vec**2) / oct_limit
    pr_oct = ar_oct**1.5 * np.sqrt(m_in/m_tot)


    ar_ss = 1. / (1-e_out_vec) * (5*np.pi * m3/(m1+m2)/np.sqrt(1-e_in_max))**(1./3.) 
    pr_ss = ar_ss**1.5 * np.sqrt(m_in/m_tot)
    
    if PLOT_PERIOD_RATIO  == False:
    
        plt.figure(figsize=(20,10))
        #to fix legend  duplicating red line
        #red
        plt.semilogy(e_out_vec, ar_stable, label='Stability limit', color=colors[2], ls=lines[0], lw=2)
        
        #green 
        w_oct = np.arange(len(ar_oct))[ar_oct > ar_ss]
        plt.semilogy(e_out_vec[w_oct], ar_oct[w_oct], label='Octupole limit', color=colors[1], ls=lines[2], lw=3)
        
        #blue
        w_ss = np.arange(len(ar_ss))[ar_ss > ar_stable]
        plt.semilogy(e_out_vec[w_ss], ar_ss[w_ss], label='Semisecular limit', color=colors[0], ls=lines[1], lw=3)
        
        #red
        plt.semilogy(e_out_vec, ar_stable, color=colors[2], ls=lines[0], lw=2)
        
        plt.fill_between(e_out_vec, 1, ar_stable, alpha=0.7, color='k')
        plt.fill_between(e_out_vec[w_ss], ar_stable[w_ss], ar_ss[w_ss], alpha=0.45, color='k')
        max_ss_oct = np.maximum(ar_ss,ar_stable)
        plt.fill_between(e_out_vec[w_oct], max_ss_oct[w_oct], ar_oct[w_oct], alpha=0.1, color='k')
    
        
        
        alpha_kozai = 1
        
        #    ar_kl_100 = (100**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.)
        ar_kl_1000 = (1000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.)
        ar_kl_10000 = (10000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.)
        ar_kl_100000 = (100000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.)
        ar_kl_1000000 = (1000000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.)
        
        w_kl_1000 = np.arange(len(e_out_vec))[ar_kl_1000>ar_oct]
        w_kl_10000 = np.arange(len(e_out_vec))[ar_kl_10000>ar_oct]
        w_kl_100000 = np.arange(len(e_out_vec))[ar_kl_100000>ar_oct]
        w_kl_1000000 = np.arange(len(e_out_vec))[ar_kl_1000000>ar_oct]
        plt.semilogy(e_out_vec[w_kl_1000], ar_kl_1000[w_kl_1000], color='0.6', ls=':', lw=2)
        plt.annotate('$10^{4}P_{out}$',(e_out_vec[w_kl_1000][0],ar_kl_1000[w_kl_1000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
      
        plt.semilogy(e_out_vec[w_kl_10000], ar_kl_10000[w_kl_10000], color='0.6', ls=':', lw=2)
        plt.annotate('$10^{5}P_{out}$',(e_out_vec[w_kl_10000][0],ar_kl_10000[w_kl_10000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
        
        plt.semilogy(e_out_vec[w_kl_100000], ar_kl_100000[w_kl_100000], color='0.6', ls=':', lw=2)
        plt.annotate('$10^{6}P_{out}$',(e_out_vec[w_kl_100000][0],ar_kl_100000[w_kl_100000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
        
        plt.semilogy(e_out_vec[w_kl_1000000], ar_kl_1000000[w_kl_1000000], color='0.6', ls=':', lw=2, label=r'LK timescale [$10^4,10^5,10^6,10^7]P_{\rmout}$')
        plt.annotate('$10^{7}P_{out}$',(e_out_vec[w_kl_1000000][0],ar_kl_1000000[w_kl_1000000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
    
        
        plt.errorbar(e_out_objects,a_ratio_objects,xerr=errors,fmt='o',ecolor='black',mfc='black',mec='black',capsize=5)
        
        for i in range(len(labels)):
            plt.annotate(labels[i],(e_out_objects[i],a_ratio_objects[i]),fontsize=12,textcoords='offset pixels',xytext=(0,10))
        
        ax = plt.gca()
        majorLocatory   = MultipleLocator(0.2)
        minorLocatory   = MultipleLocator(0.05)
        ax.xaxis.set_major_locator(majorLocatory)
        ax.xaxis.set_minor_locator(minorLocatory)
        plt.xlim((0,1))
        
        ax.set_xlabel(r'Eccentricity of outer orbit ($e_{\rm out}$)')
        ax.xaxis.set_label_coords(.5, -0.08)
        plt.ylabel(r'Ratio of semimajor axes ($a_{\rm out}/a_{\rm in}$)')
        plt.ylim((2, 1e4))
        
        plt.legend(bbox_to_anchor=[0.05, 1.],loc='upper left')
        if SAVE_FIG:
            plt.savefig(fig_name)
        plt.show()
        
    else:
        
            
        plt.figure(figsize=(20,10))
        
        #to fix legend  duplicating red line
        #red
        plt.semilogy(e_out_vec, pr_stable, label='Stability limit', color=colors[2], ls=lines[0], lw=2)
        
        #green 
        w_oct = np.arange(len(pr_oct))[pr_oct > pr_ss]
        plt.semilogy(e_out_vec[w_oct], pr_oct[w_oct], label='Octupole limit', color=colors[1], ls=lines[2], lw=3)
        
        #blue
        w_ss = np.arange(len(pr_ss))[pr_ss > pr_stable]
        plt.semilogy(e_out_vec[w_ss], pr_ss[w_ss], label='Semisecular limit', color=colors[0], ls=lines[1], lw=3)
        
        #red
        plt.semilogy(e_out_vec, pr_stable, color=colors[2], ls=lines[0], lw=2)
        
        plt.fill_between(e_out_vec, 1, pr_stable, alpha=0.7, color='k')
        plt.fill_between(e_out_vec[w_ss], pr_stable[w_ss], pr_ss[w_ss], alpha=0.45, color='k')
        max_ss_oct = np.maximum(pr_ss,pr_stable)
        plt.fill_between(e_out_vec[w_oct], max_ss_oct[w_oct], pr_oct[w_oct], alpha=0.1, color='k')
    
        
        
        alpha_kozai = 1
        
        pr_kl_1000 = ((1000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.))**1.5*np.sqrt(m_in/m_tot)
        pr_kl_10000 = ((10000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.))**1.5*np.sqrt(m_in/m_tot)
        pr_kl_100000 = ((100000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.))**1.5*np.sqrt(m_in/m_tot)
        pr_kl_1000000 = ((1000000**2 /alpha_kozai**2 * m3*m3/m_in/m_tot * (1-e_out_vec**2)**-3)**(1./3.))**1.5*np.sqrt(m_in/m_tot)
        
        w_kl_1000 = np.arange(len(e_out_vec))[pr_kl_1000>pr_oct]
        w_kl_10000 = np.arange(len(e_out_vec))[pr_kl_10000>pr_oct]
        w_kl_100000 = np.arange(len(e_out_vec))[pr_kl_100000>pr_oct]
        w_kl_1000000 = np.arange(len(e_out_vec))[pr_kl_1000000>pr_oct]
        
        plt.semilogy(e_out_vec[w_kl_1000], pr_kl_1000[w_kl_1000], color='0.6', ls=':', lw=2)
        plt.annotate('$10^{4}P_{out}$',(e_out_vec[w_kl_1000][0],pr_kl_1000[w_kl_1000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
      
        plt.semilogy(e_out_vec[w_kl_10000], pr_kl_10000[w_kl_10000], color='0.6', ls=':', lw=2)
        plt.annotate('$10^{5}P_{out}$',(e_out_vec[w_kl_10000][0],pr_kl_10000[w_kl_10000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
        
        plt.semilogy(e_out_vec[w_kl_100000],pr_kl_100000[w_kl_100000], color='0.6', ls=':', lw=2)
        plt.annotate('$10^{6}P_{out}$',(e_out_vec[w_kl_100000][0],pr_kl_100000[w_kl_100000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
        
        plt.semilogy(e_out_vec[w_kl_1000000], pr_kl_1000000[w_kl_1000000], color='0.6', ls=':', lw=2, label=r'LK timescale [$10^4,10^5,10^6,10^7]P_{\rmout}$')
        plt.annotate('$10^{7}P_{out}$',(e_out_vec[w_kl_1000000][0],pr_kl_1000000[w_kl_1000000][0]),fontsize=12,color='gray',textcoords='offset pixels',xytext=(0,5))
    
        a_ratio_objects = np.array(a_ratio_objects)
        p_ratio_objects = a_ratio_objects ** 1.5*np.sqrt(m_in/m_tot)
        
        plt.errorbar(e_out_objects,p_ratio_objects,xerr=errors,fmt='o',ecolor='black',mfc='black',mec='black',capsize=5)
        
        for i in range(len(labels)):
            plt.annotate(labels[i],(e_out_objects[i],p_ratio_objects[i]),fontsize=12,textcoords='offset pixels',xytext=(0,10))
        
        ax = plt.gca()
        majorLocatory   = MultipleLocator(0.2)
        minorLocatory   = MultipleLocator(0.05)
        ax.xaxis.set_major_locator(majorLocatory)
        ax.xaxis.set_minor_locator(minorLocatory)
        plt.xlim((0,1))
        
        ax.set_xlabel(r'Eccentricity of outer orbit ($e_{\rm out}$)')
        ax.xaxis.set_label_coords(.5, -0.08)
        plt.ylabel(r'Ratio of orbital periods ($P_{\rm out}/P_{\rm in}$)')
        plt.ylim((2, 1e6))
        
        plt.legend(bbox_to_anchor=[0.05, 1.],loc='upper left')
        if SAVE_FIG:
            plt.savefig(fig_name)
        plt.show()
        






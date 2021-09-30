from __future__ import print_function
import numpy as np
import matplotlib, sys
sys.path.append("/home/vgupta/resources")
from rcparams import rc_params
matplotlib.rcParams.update(rc_params)
import matplotlib.pyplot as plt
plt.rcParams.update(rc_params)

from matplotlib.patches import Patch
import argparse
from Survey import Survey

from mpl_toolkits.mplot3d import axes3d
from astropy.cosmology import WMAP9


def dms_width(dm, survey):                                                                      
    #tau (usec) = 8.3 * BW (MHz) * DM * f**-3 (GHz)             from Cosmos webpage
    return 8.3 * (survey.CHW*1e-6) * dm * (survey.CFREQ*1e-9)**-3 * 1e-6
                       

def dm_to_dist(dm):
  z = (dm)/1000.
  dist = WMAP9.luminosity_distance(z).si.value
  return dist

def dm_to_z(dm):
  return dm/1000.

def fluence_limit_for_telescope(t, iw_axis, dm_axis, snr_limit=10, idealize=False,tsamp = None):
    survey = Survey(t)
    if idealize:
        survey.idealize()
    iw_axis = 10**iw_axis.copy()
    
    if tsamp is None:
        tsamp = survey.TSAMP
        
    ofl_dm_plane = []
    for dm in dm_axis:
        dms_width_axis = dms_width(dm, survey)
        samp_axis = tsamp * np.ones_like(iw_axis)
        '''
        ow_axis = np.where(dms_width_axis > iw_axis, dms_width_axis, iw_axis)
        
        ow_axis = np.where(ow_axis > samp_axis, ow_axis, samp_axis)
        '''

        ow_axis = np.sqrt(iw_axis**2 + samp_axis**2 + dms_width_axis**2)
        
        #ifl_axis = np.log10(snr_limit * 10**iw_axis * survey.SEFD / np.sqrt(10**iw_axis * survey.BW * survey.NPOL))
        ofl_axis = np.log10(snr_limit * ow_axis*1e3 * survey.SEFD / np.sqrt(ow_axis * survey.BW * survey.NPOL))  #Jy ms
        #ofl_axis = np.log10(snr_limit * ow_axis * survey.SEFD / np.sqrt(ow_axis * survey.BW * survey.NPOL) * 1e-26)  #Joules / m^2 / Hz
        ofl_dm_plane.append(ofl_axis)

    return np.array(ofl_dm_plane)

def fix_mpl_bug(s):
  s._facecolors2d = s._facecolors3d
  s._edgecolors2d = s._edgecolors3d

def main():

    iw_axis = np.arange(args.wr[0], args.wr[1], 0.01)
    dm_axis = np.arange(args.dmr[0], args.dmr[1], 2)
    ofl_dm_plane_actual = fluence_limit_for_telescope(args.t, iw_axis, dm_axis) 
    ofl_dm_plane_ideal = fluence_limit_for_telescope(args.t, iw_axis, dm_axis, idealize=True) 

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(iw_axis, dm_axis)
    si = ax.plot_surface(X, Y, ofl_dm_plane_ideal, cmap='viridis') 
    sc = ax.plot_surface(X, Y, ofl_dm_plane_actual, cmap='plasma') 
    label_Jy = [Patch(facecolor='pink', edgecolor='r', label="Actual"), Patch(facecolor='#00846b', edgecolor='k', label='Ideal')]
  
    #fix_mpl_bug(si)
    #fix_mpl_bug(sc)
    
    ax.set_xlabel("log(sky width) [s]", labelpad=30)
    ax.set_ylabel("DM [pc/cc]", labelpad=30)
    ax.set_zlabel("log(Ob. fluence) [Jy ms]", labelpad=20)
    ax.legend(handles = label_Jy, fontsize=25)
    ax.view_init(azim=-116, elev=30)
    ax.tick_params(axis='x', labelsize=22)
    ax.tick_params(axis='y', labelsize=22)
    ax.tick_params(axis='z', labelsize=22)
    ax.xaxis.label.set_size(30)
    ax.yaxis.label.set_size(30)
    ax.zaxis.label.set_size(30)
    #fig.legend()
    fig.savefig("{}_Ob_fl_dm_w_threshold3D.png".format(args.t))
   
    survey = Survey(args.t)
    dist = dm_to_dist(dm_axis)
    z = dm_to_z(dm_axis)
    Joules_ofl_dm_plane_actual = ofl_dm_plane_actual -29    #(-3 for 1ms and -26 for Jy to Joules)
    Joules_ofl_dm_plane_ideal = ofl_dm_plane_ideal -29    #(-3 for 1ms and -26 for Jy to Joules)
    ifl_dm_plane_actual = np.log10(10**Joules_ofl_dm_plane_actual * (4 * np.pi * dist**2 / (1 + z) * survey.CFREQ)[:, None] )   #Joules/Hz
    ifl_dm_plane_ideal = np.log10(10**Joules_ofl_dm_plane_ideal * (4 * np.pi * dist**2 / (1 + z) * survey.CFREQ)[:, None] )   #Joules/Hz
   
    fig2 = plt.figure(figsize=(10,8))
    ax2 = fig2.add_subplot(111, projection='3d')
    isi= ax2.plot_surface(X, Y, ifl_dm_plane_ideal, cmap='viridis') 
    isc= ax2.plot_surface(X, Y, ifl_dm_plane_actual, cmap='plasma') 
    
    label_J = [Patch(facecolor='orange', edgecolor='r', label="Current"), Patch(facecolor='#00846b', edgecolor='k', label='Ideal')]
    #fix_mpl_bug(isi)
    #fix_mpl_bug(isc)

    ax2.set_xlabel("log(sky width) [s]", labelpad=30)
    ax2.set_ylabel("DM [pc/cc]", labelpad=30)
    ax2.set_zlabel("log(In. energy) [Joules]", labelpad=10)
    ax2.legend(handles=label_J, fontsize=25)
    ax2.view_init(azim=-146, elev=20)
    ax2.tick_params(axis='x', labelsize=22)
    ax2.tick_params(axis='y', labelsize=22)
    ax2.tick_params(axis='z', labelsize=22)
    ax2.xaxis.label.set_size(30)
    ax2.yaxis.label.set_size(30)
    ax2.zaxis.label.set_size(30)
    fig2.savefig("{}_In_fl_dm_w_threshold3D.png".format(args.t))
    
    DM_frb = 715. #pc/cc
    frb_sig_835 = 11.35e-6 #s
    frb_tau_835 = 21.43e-6 #s
    frb_tau_cfreq = frb_tau_835 * (835.5e6 / survey.CFREQ)**4
    sky_w_frb = np.sqrt(frb_tau_cfreq**2 + frb_sig_835**2)
    #sky_w_frb = 24.3e-6 #s
    pl_index = args.pli - 1
    print("DM_frb, sky_w_frb, pli == {0}, {1}, {2}".format(DM_frb, sky_w_frb, pl_index))

    
    det_fraction = np.zeros_like(ofl_dm_plane_ideal) 
    for i_dm in range(len(dm_axis)):                                                   
        for i_iw in range(len(iw_axis)):                                    
            fluence_cutoff_ideal = ofl_dm_plane_ideal[i_dm, i_iw]
            fluence_cutoff_actual = ofl_dm_plane_actual[i_dm, i_iw]
            area_until_ideal = -1. * ((10**fluence_cutoff_ideal)**(1 + pl_index) / (1 + pl_index) )
            area_until_actual = -1. * ((10**fluence_cutoff_actual)**(1 + pl_index) / (1 + pl_index) )
            ratio_of_areas = area_until_actual  /area_until_ideal
            det_fraction[i_dm, i_iw] = ratio_of_areas
            if (DM_frb - 1.5 < dm_axis[i_dm] < DM_frb + 1.5) and (np.log10(sky_w_frb) - 0.005 < iw_axis[i_iw] < np.log10(sky_w_frb) + 0.005):
              print("i_dm = {0}, i_w = {1}, fluence_cutoff_actual = {2}, fluence_cutoff_ideal = {3}, det_fraction = {4}".format(dm_axis[i_dm], iw_axis[i_iw], fluence_cutoff_actual, fluence_cutoff_ideal, ratio_of_areas)) 

    plt.figure(figsize=(15,12))
    plt.imshow(det_fraction, aspect='auto', interpolation=None, origin = 'lower', extent=[iw_axis[0], iw_axis[-1], dm_axis[0], dm_axis[-1]])

    iw_frb_index = np.argmin(np.abs(10**iw_axis - sky_w_frb))
    dm_frb_index = np.argmin(np.abs(dm_axis - DM_frb))
    print("det_fraction.shape = {}, iw_frb_index = {}, dm_frb_index = {}".format(det_fraction.shape, iw_frb_index, dm_frb_index))
    print("{} survey would miss {}% of FRBs at the DM and sky width of FRB191107B".format(args.t, (1-det_fraction[dm_frb_index, iw_frb_index])*100.))
    print("{} survey would miss {}% of FRBs at the DM and sky width of FRB191107B".format(args.t, (1-det_fraction[iw_frb_index, dm_frb_index])*100.))


    telescope_names = np.array(['UTMOST', 'SUPERB', 'ASKAP_FE', 'ASKAP_ICS', 'CHIME'])
    telescope_indices = ['A', 'B', 'E', 'D', 'C']
    idx = np.where(telescope_names == args.t)[0][0]
    print(idx)
    telescope_idx = telescope_indices[idx]

    plt.plot(np.log10(sky_w_frb), DM_frb, '*', c='#FC76EC', ms=25, label="FRB191107B at {0:.2f} GHz".format(survey.CFREQ/1e9))
    cbar= plt.colorbar()
    cbar.set_label("Fraction of FRBs detected")
    plt.xlabel(r"log($w_{sky}$) [s]")
    plt.ylabel("Total DM [pc/cc]")
    #plt.title("Power law index: {}".format(args.pli))
    plt.title("({0}) {1}".format(telescope_idx.lower(), str(args.t).replace('_', ' ')), fontsize=40)
    plt.legend(fontsize=30)
    plt.tight_layout()
    plt.savefig("Det_fraction_{}_pli_{}.png".format(args.t, args.pli))

    plt.show()

if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument("-t", type=str, help="Telescope/Survey to simulate [UTMOST, SUPERB, ASKAP_FE, ASKAP_ICS, CHIME](def:UTMOST)", default="UTMOST")
    a.add_argument("-pli", type=float, help="Power law index to simulate (def = -1.8)", default=-1.8)
    a.add_argument("-dmr", type=float, nargs=2, help="DM range to simulate - edge values delimetered by comma (def: 0.1, 3000)", default=[0.1, 3000])
    a.add_argument("-wr", type=float, nargs=2, help="log10(Width) range to simulate (in sec) - edge values delimetered by comma (def: -6, -2)", default=[-6, -2])
    args = a.parse_args()
    main()


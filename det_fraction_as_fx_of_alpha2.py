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
    
    if tsamp is None:
        tsamp = survey.TSAMP
    
    dms_w = dms_width(dm_axis, survey)
    ow = np.sqrt(iw_axis**2 + tsamp**2 + dms_w**2)
    print("ow = {0}, iw_axis = {1}, tsamp = {2}, dms_w = {3}".format(ow, iw_axis, tsamp, dms_w))
    o_fl_limit = np.log10(snr_limit * ow*1e3 * survey.SEFD / np.sqrt(ow * survey.BW * survey.NPOL)) # Jy ms

    return o_fl_limit


def main():
  #plt.figure(figsize=(10, 8))
  for tel in [r"UTMOST", r"SUPERB", r"CHIME", r"ASKAP_ICS"]:

    survey = Survey(tel)
    DM_frb = 715. #pc/cc
    DM_EG = DM_frb
    #DM_EG = 537. #pc/cc
    
    z_FRB = dm_to_z(DM_EG)
    dist_FRB = dm_to_dist(DM_EG)

    frb_sig_835 = 11.35e-6 #s
    frb_tau_835 = 21.43e-6 #s
    frb_tau_cfreq = frb_tau_835 * (835.5e6 / survey.CFREQ)**4
    sky_w_frb = np.sqrt(frb_tau_cfreq**2 + frb_sig_835**2)

    print("DM_frb, sky_w_frb == {0}, {1}".format(DM_frb, sky_w_frb))
    ofl_dm_plane_actual = fluence_limit_for_telescope(tel, sky_w_frb, DM_frb) 
    ofl_dm_plane_ideal = fluence_limit_for_telescope(tel, sky_w_frb, DM_frb, idealize = True) 

    print("ofl_dm_plane_actual_limit_at_FRB191107B: {0}".format(ofl_dm_plane_actual))
    print("ofl_dm_plane_ideal_limit_at_FRB191107B: {0}".format(ofl_dm_plane_ideal))
    Joules_ofl_dm_plane_actual = ofl_dm_plane_actual -29
    Joules_ofl_dm_plane_ideal = ofl_dm_plane_ideal -29

    ifl_dm_plane_actual = np.log10(10**Joules_ofl_dm_plane_actual * (4 * np.pi * dist_FRB**2 / (1 + z_FRB) * survey.CFREQ) )   #Joules/Hz
    ifl_dm_plane_ideal = np.log10(10**Joules_ofl_dm_plane_ideal * (4 * np.pi * dist_FRB**2 / (1 + z_FRB) * survey.CFREQ) )   #Joules/Hz 
      
    det_fractions = []
    alphas = np.arange(-1.5, -3.5,-0.1)
    method = "old"
    if method == "old":
      cutoff_plane_actual = ofl_dm_plane_actual
      cutoff_plane_ideal = ofl_dm_plane_ideal
    elif method == "new":
      cutoff_plane_actual = ifl_dm_plane_actual
      cutoff_plane_ideal = ifl_dm_plane_ideal

    for pl_index in alphas:
      area_until_ideal = -1. * ((10**cutoff_plane_ideal)**(1 + pl_index) / (1 + pl_index) )
      area_until_actual = -1. * ((10**cutoff_plane_actual)**(1 + pl_index) / (1 + pl_index) )
      det_fraction = area_until_actual / area_until_ideal
      det_fractions.append(det_fraction)

    tname = tel
    if tname == r"ASKAP_ICS":
      tname = r"ASKAP ICS"
    plt.plot(alphas+1, det_fractions, label=tname)

  #plt.yscale('log')
  #plt.ylim(0.01, 1.1)
  plt.xlabel(r"Power law index ($\gamma$)")
  plt.ylabel("Detected fraction")
  plt.axvline(-1.8, ls='--', c='k', label=r"$\gamma$ for FRB121102A")
  plt.legend(fontsize=30)
  plt.tight_layout()
  if args.s:
    plt.savefig("Det_fraction_fx_alpha.png")
     
  plt.show()

if __name__ == '__main__':
    a = argparse.ArgumentParser()
    #a.add_argument("-t", type=str, help="Telescope/Survey to simulate [UTMOST, SUPERB, ASKAP_FE, ASKAP_ICS, CHIME](def:UTMOST)", default="UTMOST")
    a.add_argument("-s", action='store_true', default=False, help="Save fig (def: False)")
    args = a.parse_args()
    main()


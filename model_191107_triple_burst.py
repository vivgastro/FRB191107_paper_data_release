import numpy as np
import bilby
import sys, argparse, glob
from matplotlib import rcParams
import matplotlib.pyplot as plt
rcParams['text.usetex'] = False
#rcParams['text.latex.preamble'] = r'\newcommand{\mathdefault}[1][]{}'

ts = 0

def model(x, **kwargs): 
    A1 = kwargs["A1"]
    A2 = kwargs["A2"]
    A3 = kwargs['A3']

    sig1 = kwargs["sig1"]
    sig2 = kwargs["sig2"]
    sig3 = kwargs['sig3']

    tau = kwargs["tau"]
    c1 = kwargs["c1"]
    c2 = kwargs["c2"]
    c3 = kwargs['c3']

    y = A1 * np.exp(-1*(x - c1)**2 / (2*sig1**2)) + A2 * np.exp(-1*(x-c2)**2 / (2*sig2**2)) + A3 * np.exp(-1*(x-c3)**2 / (2*sig3**2))
    peak_loc = np.argmax(y) 
    peak_val = y[peak_loc] 
    y = convolve_exp(y, tau) + sys.float_info[3]
    return norm(y, peak_val, peak_loc) 
 
def convolve_exp(x, tau): 
    t = np.arange(tau/ts * 10) *ts
    kernel = np.exp(-t / tau) 
    y = np.convolve(x, kernel, mode='full') 
    return y 
     
def norm(d, peak_val, peak_loc): 
    cpeak_loc = np.argmax(d) 
    cpeak_val = d[cpeak_loc] 
    loc_diff = peak_loc - cpeak_loc 
    y = d.copy()
    #y = np.roll(d, loc_diff) 
    return y * peak_val / cpeak_val

class VGLikelihood(bilby.Likelihood):
    def __init__(self, x, data, model_func):
        self.x = x
        self.N = len(x)
        self.data = data
        self.model_func = model_func
        self.params = {
            'A1' : None,
            'A2' : None,
            'A3' : None,
            'c1' : None,
            'c2' : None,
            'c3' : None,
            'sig1' : None,
            'sig2' : None,
            'sig3' : None,
            'tau' : None,
            'noise': None}

        super().__init__(parameters = self.params)

    def log_likelihood(self):
        noise = self.params["noise"]
        res = self.data - self.model_func(self.x, **self.params)[:len(self.data)]
        return  -0.5 * ( np.sum((res / noise)**2 + np.log(2*np.pi * noise**2) ) )

def interp(data, ix):
   Np = len(data)
   new_data = np.zeros((Np-1)*(ix-1) + Np)
   for ii in range(Np):
     if ii == Np-1:
         break
     new_data[ii*ix:(ii+1)*ix]  = np.linspace(data[ii], data[ii+1], ix, endpoint=False)
   new_data[-1] = data[-1]
   return new_data

def get_tsamp(frb):
  tsamp_file = "FRB_tsamps.txt"
  frbs = np.loadtxt(tsamp_file, usecols=(0), dtype=str)
  tsamps = np.loadtxt(tsamp_file, usecols=(1))

  frb_id = np.where(frbs == frb)[0]
  tsamp = tsamps[frb_id]
  return tsamp[0]

def main():
  fname = glob.glob("{0}/{0}*npy".format(args.FRB))
  global ts
  ts = get_tsamp(args.FRB)
  print("Processing {}".format(fname))
  data = np.load(fname[0])
  data = data - np.median(data)
  frb_location = np.argmax(data)
  ns = len(data)
  #data = data[int(frb_location - ns/2): int(frb_location + ns/2)]
  #data = interp(data, 4)
  
  Nsamples = len(data)
  print(data.shape, ns, frb_location)
  print(ts)
  x = np.arange(Nsamples)*ts
  
  #noise = 1.5
  #data = model(x, 10, Nsamples/2, 5, 5) + np.random.normal(0, noise, Nsamples)
  
  priors = dict()
  priors["A1"] = bilby.core.prior.Uniform(0.000001, 0.3, name=r"$Amp_1$")
  priors["c1"] = bilby.core.prior.Uniform((ns/2+1 - 100)*ts, (ns/2+1 - 40)*ts,  name=r"$t_{1}$")
  priors["sig1"] = bilby.core.prior.Uniform((0.00001)*ts, 100*ts, name=r"$\sigma_1$")
  
  priors["A2"] = bilby.core.prior.Uniform(0.000001, 0.4, name=r"$Amp_2$")
  priors["c2"] = bilby.core.prior.Uniform((ns/2+1 - 26)*ts, (ns/2+1 - 19)*ts,  name=r"$t_{2}$")
  priors["sig2"] = bilby.core.prior.Uniform((0.00001)*ts, 10*ts, name=r"$\sigma_2$")
  
  priors["A3"] = bilby.core.prior.Uniform(0.4, 1.00, name=r"$Amp_3$")
  priors["c3"] = bilby.core.prior.Uniform((ns/2+1 - 3)*ts, (ns/2+1 + 2)*ts,  name=r"$t_{3}$")
  priors["sig3"] = bilby.core.prior.Uniform((0.00001)*ts, 4*ts, name=r"$\sigma_3$")

  priors["tau"] = bilby.core.prior.Uniform((0.00001)*ts, 5*ts, name=r"$\tau$")
  priors["noise"] = bilby.core.prior.Uniform(0, 0.067, name=r"$noise$")
  
  outdir = "{}/triple_peak_modelled_nsteps5000/".format(args.FRB)
  #outdir = "{}/triple_peak_modelled_dynesty/".format(args.FRB)
  likelihood_obj = VGLikelihood(x, data, model)
  
  result = bilby.run_sampler(likelihood=likelihood_obj, priors = priors,
      sampler = "emcee", nwalkers = 100, nsteps = 8000, outdir = outdir, nburn = 3000)
  
  #result =bilby.sampler.run_sampler(likelihood=likelihood_obj, priors=priors,
  #    sampler="dynesty", npoints = 1000, outdir=outdir, label=None, resume=True)
  fg = plt.figure(figsize=(8,9)) 
  fg.suptitle(r"${0}\ tsamp={1:.2f}\mu s$".format(args.FRB, ts))
  result.plot_corner(filename="./{0}/posterior_corner.png".format(outdir), fig=fg, dpi=100)

if __name__ == '__main__':
  a = argparse.ArgumentParser()
  a.add_argument("-FRB", type=str, help="Name of the FRB to model", default="FRB191107")
  args = a.parse_args()
  main()

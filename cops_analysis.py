'''
During manual assignment, extracts Cbeta and other fit parameters for a given peak using pyruvate-GRADCOPs lineshape fitting.
Also contains a module for simulation of GRADCOPs lineshapes. 

May, 2022
'''

__author__ = "Harrison Wang"
__email__ = "hwang6@fas.harvard.edu"

import numpy as np
import matplotlib.pyplot as plt
import scipy
from matplotlib.widgets import Slider, Button
import scipy.optimize
import scipy.io
from scipy.interpolate import interp1d
import nmrglue as ng
import os


###To do: add a function to analyze non-pyruvate experiments.

class cops_analyze():
    
    def __init__(self,data_strs, mode='HNCA', pyruvate_on=True, cop_num=[1,2,3,4,5]):
        
        '''        
        Import requirements 
        ___________
        
        HNCA_nocop.ft3: a 3D HNCA spectrum, processed with NMRpipe and pipetoccp.com, without applying GRADCOPs decoupling pulses.
        Can be omitted in the case of uniformly-labeled spectra.
        
        HNCA_cop%i.ft3: a series of 3D HNCA spectra processed in the same manner as HNCA_nocop.ft3. 
        See definition of cop_num and definition of self.cop_dics for format of input string.
        
        dec_profiles_named.mat: A MATLAB .mat file with 6 rows and 274 columns. The 13C chemical shift ranges over the columns,
        and the rows are: 
        row 1: 13C chemical shift. 
        rows 2-6: decoupling magnetization inversion profile of GRADCOPs 1, 3-6, in order. 
        This needs to be converted to decoupled fraction before use.
        
        
        Parameters
        __________
        cop_num: int or list
        If int, number of COPs spectra to import. Default value is 5, for GRADCOPs 1, 3-6.
        If list, the named COPs spectra to import. Default value is [1,3,4,5,6], for GRADCOPs 1, 3-6. 
        pyruvate_on: boolean, default True
        indicates if pyruvate was used for nonuniform labeling. 
        
        '''
        ###
        #initialize spectra
        ###
        self.pyr_on = pyruvate_on
        self.cop_nums = cop_num
        self.copnames = data_strs
        
        #allows us to squash cops analysis variables into an array without gaps
        self.cop_num = len(self.cop_nums)
        self.mode = mode
        self.init_spectra(self.mode)

        ###
        #initialize simulated decoupling profiles
        ###
        
        self.mat = np.loadtxt('./files/dec_profiles.csv')
        #self.mat.get('dec_profiles')[0]: frequency axis. self.mat.get('dec_profiles')[1]-[5]: decoupling profiles for gradcop 1, 3, 4, 5, 6.  
        self.cops = self.mat[0].reshape(1,-1)
        for i in self.cop_nums:
            self.cops=np.concatenate((self.cops, self.mat[i].reshape(1,-1)), axis=0)

        
        
        #trim cops to remove CO decoupling profile
        self.cops=self.cops[:,:171]
        
        #generates array of interpolation functions to determine value of decoupling for a particular Cb. 
        self.dec_interpolation = [interp1d(self.cops[0], self.cops[i+1], kind='cubic') for i in range(self.cop_num)]
        
        #for Calc_CB_old. 
        #number of interpolation points for decoupling profile fitting
        #self.interp_points=10000
        
        #initialize interpolated decoupling profiles.
        #self.gen_decoupling_profiles()
    
    ###############################################
    #SECTION 0. EXPERIMENTAL LINESHAPE EXTRACTION # 
    ###############################################
    def import_spectrum(self, filestr, dimension):
        filetype = os.path.splitext(filestr)[1]
        if filetype=='.ucsf':
            if dimension==3:
                dic, dat = ng.sparky.read_3D(filestr)
            elif dimension==2:
                dic, dat = ng.sparky.read_2D(filestr)
            unit_conv = [ng.sparky.make_uc(dic, dat, i) for i in range(dimension)] 
        elif filetype=='.ft3' or filetype=='.dat':
            if dimension==3:
                dic, dat = ng.pipe.read_3D(filestr)
            elif dimension==2:
                dic, dat = ng.pipe.read_2D(filestr)
            unit_conv = [ng.pipe.make_uc(dic, dat, i) for i in range(dimension)]
        else:
            raise ValueError('Check file format. Currently supports .ft3, .ucsf, .dat')
        return dic, dat, unit_conv
    
    def init_spectra(self, mode):
        ###
        #initialize spectra
        ###
        #import nocop HNCA data
        if self.mode == 'HNCA':
            dims = 3
        elif self.mode == 'HCA':
            dims = 2
        if self.pyr_on:
            # make unit conversion object for each axis of the nocop HNCA spectrum. 
            #Indexes correspond to dimensions as follows:
            #i=0: 15N, i=1: 13CA, i=2: 1H
            self.nocop_dic, self.nocop_dat, self.nocop_unit_convs = self.import_spectrum(self.copnames[0], dims)
        
        ######
        ######currently inefficient!
        ######
        self.cop_dics = []
        self.cop_dats = []
        self.cop_unit_convs = []
        for i in range(self.cop_num):
            imported = self.import_spectrum(self.copnames[i+self.pyr_on], dims)
            ##array of n dictionaries, 1 per COP spectrum
            self.cop_dics.append(imported[0])
            ##array of n datasets, 1 per COP spectrum. shape of array: [number of COPs, #points in w1 (15N), #pts in w2 (13C), #pts in w3 (1H)]
            self.cop_dats.append(imported[1])
            ##cop_num by 3 array, 1 nmrglue unit converter object per dimension per COP spectrum
            ##shape: [number of COPs, 3]
            self.cop_unit_convs.append(imported[2])
        self.cop_dics = np.array(self.cop_dics)
        self.cop_unit_convs = np.array(self.cop_unit_convs)

        return None
   
    
    #given a spectrum and unit conversion triple, extract 1D trace from center of peak at data_pt_ppm. tw: trace width, Hertz. 
    def extract1D(self, data_pt, spectrum, uc, sw=70, C_offset = 0, normalize=False):
        '''
        DEFINITION
        __________
        /given/ a peak center in ppm and a spectrum, 
        /returns/ the 13C peak profile. 
        
        PARAMETERS
        __________
        data_pt_ppm: list of length 3
        peak center (ppm).
        
        spectrum: list of dimension 3
        spectral data from nmrglue.
        
        uc: list of length 3 containing 3 nmrglue unit converter objects
        one nmrglue UC per dimension. #i=0: 15N, i=1: 13CA, i=2: 1H
        
        sw: float
        spectral width (Hz)
        
        normalize: boolean
        if True, the 13C peak profile is normalized by volume. 
        
        OUTPUT
        ______
        hz_vals: list
        x axis (Hz)
        
        trace: list
        intensity of 13C peak profile (arb. units)
         
        '''

        #convert data_pt_ppm to index

        if self.mode=='HNCA':
             
            idx = np.array([uc[0](data_pt[0], "ppm"), uc[1](data_pt[1]+C_offset, "ppm"), uc[2](data_pt[2], "ppm")])
            #calculate indices for trace boundary, based on tw. 
            hz_bounds = self.hz_to_idx(uc[1], sw)
            hz_vals = np.linspace(-hz_bounds, hz_bounds, num=2*hz_bounds+1)*sw/hz_bounds

            #1D slice through peak center, weight-added by tensor product to 1D 13C slices nearby (in the HN, N dimensions)
            weights = np.array([[0.6, 0.9, 1, 0.9, 0.6]]) #weights vector to compute weighted sum
            slices = np.array(spectrum[idx[0]-2:idx[0]+3, idx[1]-hz_bounds:idx[1]+hz_bounds+1,idx[2]-2:idx[2]+3])
            trace = np.tensordot(slices,weights.T@weights, axes=([0,2],[0,1]))
            trace = np.array(trace)
        elif self.mode=='HCA':
            
            idx = np.array([uc[i](data_pt[i], "ppm") for i in range(len(uc))])
                
            #calculate indices for trace boundary, based on tw. 
            hz_bounds = self.hz_to_idx(uc[1], sw)
            hz_vals = np.linspace(-hz_bounds, hz_bounds, num=2*hz_bounds+1)*sw/hz_bounds

            #1D slice through peak center, weight-added by tensor product to 1D 13C slices nearby (in the HN, N dimensions)
            weights = np.array([[0.4, 0.6, 1, 0.6, 0.4]]) #weights vector to compute weighted sum
            slices = spectrum[idx[0]-hz_bounds:idx[0]+hz_bounds+1, idx[1]-2:idx[1]+3]
            trace = slices@weights.T
            trace = np.array(trace)
        
        normalizer=np.sum(np.abs(trace))
        
        if normalize:
            #normalize by peak volume
            trace=trace/normalizer*(2*hz_bounds+1)
        return hz_vals, trace
    
    
    def hz_to_idx(self, uconv, Hz):
        return uconv(0, "Hz")-uconv(Hz, "Hz")
    
    ##################################
    #SECTION 1. LINESHAPE SIMULATION # 
    ##################################
    
    #### should be vectorizable.
    def lineshape(self, w, k_abs, fraction, c, j, lwid):
        '''
        DEFINITION
        __________
        /given/ an input frequency axis w and peak profile parameters:
        /returns/ a simulated 13C peak profile. 
        
        
        PARAMETERS
        __________
        
        w: list
        frequency axis of the 13C peak profile, in Hz, centered at 0. 
        
        k_abs: float
        absolute height of fully-decoupled gaussian, in arbitrary and normalized units. 
        
        fraction: float between 0 and 1
        Cb recoupling fraction. 
        
        c: float
        offset from center c of the the 13C peak profile, in Hz. 
        
        j: float
        J coupling between CA and CB, in Hz. 
        
        lwid: float
        13C linewidth, in Hz. 
        
        OUTPUT
        _______
        
        C_profile: list
        13C peak profile, in arbitrary intensity units. 1 intensity value per point on the frequency axis (w). 
        
        '''
        C_profile = k_abs/2*fraction*np.exp(-(w-c-j)**2/(2*lwid**2))+k_abs*(1-fraction)*np.exp(-(w-c)**2/(2*lwid**2))+k_abs/2*fraction*np.exp(-(w-c+j)**2/(2*lwid**2))
        return C_profile
    
    
    def lineshape_Cb(self, w, k_abs, pyr_fraction, Cb, c, j, lwid):
        '''
        DEFINITION
        __________
        /given/ a value of CB chemical shift, pyruvate Cb labeling fraction, and a frequency axis w, along with peak profile parameters:
        /returns/ a self.cop_num*len(w) list containing predicted lineshapes for the no-cop spectrum (index 0 of the lineshapes variable
        pre-reshaping) and each GRADCOP (index 1-self.cop_num). This is used for fitting experimental lineshapes. 
        
        PARAMETERS
        __________
        
        peak profile parameters (w, k_abs, c, j, lwid) defined as in the lineshape function. 
        
        Cb: float
        CB chemical shift in ppm. 
        
        pyr_fraction: float between 0 and 1
        CB labeling fraction for the given resonance.
        
        OUTPUT
        ______
        lineshapes: list
        a list of the predicted 13C peak profile for each GRADCOP spectrum, appended end-to-end.
        
        '''
        
        #calculates the decoupling fraction for each GRADCOP given the Cb. 
        decoupling_fractions = [(1+self.dec_interpolation[i](Cb))/2 for i in range(self.cop_num)]
        #generates lineshape
        lineshapes = np.array([self.lineshape(w, k_abs, pyr_fraction*decoupling_fractions[i], c, j, lwid) for i in range(self.cop_num)])
        #reshapes to provide a target list for fitting experimental lineshapes. 
        lineshapes=lineshapes.reshape(-1)
        return lineshapes
        

    ###########################
    #SECTION 2. MAIN FUNCTION # 
    ###########################
    
    #takes in data point coordinates in ppm, and outputs the triangulated Cb, directly optimized over lineshapes. 
    def CalcCB(self, data_pt, fit_tol=10, simple_output=True):
        '''
        DEFINITION
        __________
        /given/ a single peak center in ppm:
        /returns/ the triangulated Cb value.
        
        PARAMETERS
        __________
        data_pt: list of length 3
        15N, 13C, 1H chemical shift of the peak center (ppm).
        
        fit_tol: float, default value 3
        tolerance in arbitrary units for fitting COPs spectral parameters given the nocop paramter fit. 
        
        simple_output: boolean, default value True
        if True, outputs only CB value. if False, outputs entire fit parameters as well as a credence value.
        
        OUTPUT
        ______
        best_params: list
        fit parameters. See description of lineshape_Cb for variable definitions.  
        index 0: best fit of the k_abs variable
        index 1: best fit of the pyr_fraction variable
        index 2: best fit of CB (ppm)
        index 3: best fit of the c variable
        index 4: best fit of the j variable
        index 5: best fit of the linewidth variable
        
        1/min_sq: float
        fit credence. a large value (>50) indicates the fit is very good. A small value (~5) indicates the fit is poor. 
        
        '''
        
        if (self.mode=='HNCA' and len(data_pt)!=3) or self.mode=='HCA' and len(data_pt)!=2:
            raise ValueError('Format of peak shift list is incorrect.')
        
        if self.pyr_on:
            #strategy: optimize Cb coupling fraction in parallel to extract Cb directly
            
            ####
            ##have to remove C_offset before alpha
            ####
            hz, nocop_trace = self.extract1D(data_pt, self.nocop_dat, self.nocop_unit_convs, C_offset=-0.09, normalize=True)
            nocop_params = self.lineshape_fit(hz, nocop_trace)
            #unpack some nocop lineshape fit parameters
            pyr_fraction = nocop_params[1]
            j_ab = nocop_params[3]
            lwid=nocop_params[4]

        else: #for uniformly-labeled Cb. 
            j_ab = 40
            lwid = 10

            

        #reshapes experimental lineshape
        hz = self.extract1D(data_pt, self.cop_dats[1], self.cop_unit_convs[1], normalize=True)[0]
        cop_1Ds = np.array([self.extract1D(data_pt, self.cop_dats[i], self.cop_unit_convs[i], normalize=True)[1] for i in range(self.cop_num)])
        cop_1Ds = cop_1Ds.reshape(-1)

        '''bounds: [k_abs, pyr_fraction, Cb, c, j, lwid]'''

        for i in range(7):
            
            #initializes CB to a ppm value between 10 and 45, inclusive. 
            prior_cb = 15+i*5
            
            #initialize prior and bounds of fitting.
            if self.pyr_on:
                p = [15,pyr_fraction,prior_cb,0,j_ab,lwid]
                bounds = ([0,pyr_fraction-fit_tol/200,9,-5,j_ab-fit_tol/4,lwid-fit_tol/4],[40,pyr_fraction+fit_tol/200,46,5,j_ab+fit_tol/4,lwid+fit_tol/4])
                params = self.lineshape_Cb_fit(hz, cop_1Ds, prior=p, bounding=bounds)
            else:
                p = [15,1,prior_cb,0,j_ab,lwid]
                bounds = ([0,0.999,9,-30,j_ab-20,lwid-5],[40,1,46,30,j_ab+20,lwid+5])
                params = self.lineshape_Cb_fit(hz, cop_1Ds, prior=p, bounding=bounds)
                
            
            #measures the squared error between experimental and simulated peak profile
            sq_error = np.sum((self.lineshape_Cb(hz, *params)-cop_1Ds)**2)
            
            #determines which of the initial CB values produces the smallest error.  
            if i ==0:
                min_sq = sq_error
                best_params = params
            if sq_error<min_sq:
                best_params = params
                min_sq = sq_error
                
        if simple_output: 
            return best_params[2] #returns CB shift value in ppm
        else:
            return best_params, 1/min_sq #returns every fit parameter: CB shift (ppm), linewidth, J coupling; as well as 1/error of the CB estimate

    '''old function that uses the strategy of fitting individual lineshapes separately.'''
    #fit_tol: tolerance (Hz) for COP trace fitting
    def CalcCB_old(self, data_pt, fit_tol = 2):
        #strategy: fit the no COP trace, then use fit parameters to bound COP fitting
        hz, nocop_trace = self.extract1D(data_pt, self.nocop_dat, self.nocop_unit_convs, normalize=True)
        t = self.lineshape_fit(hz, nocop_trace)

        
        #COP trace fitting bounded by determined fit parameters
        bounds=([0,0,-5,t[3]-fit_tol,t[4]-fit_tol],[40,1,5,t[3]+fit_tol,t[4]+fit_tol])

        frac = []
        for i in range(len(self.cop_dats)):
            hz, line,_ = self.extract1D(data_pt, self.cop_dats[i], self.cop_unit_convs[i], normalize=True)
            t1 = self.lineshape_fit(hz, line, bounding = bounds)
            frac = np.append(frac, t1[1])
        
        #determine fraction of decoupled intensity
        frac = frac/t[1]
        errors, pt = self.min_error_calc(frac)
        #if 1/np.min(errors) is low, then fit is not good. 
        return pt, 1/np.min(errors)
    
    def min_error_calc(self, fracs_vector):
        error=np.array([])
        
        fracs_vectors = np.reshape(fracs_vector,[-1,1])@np.reshape(np.ones(self.interp_points),[1,-1])
        errors = np.sum((self.predicts-fracs_vectors)**2,axis=0)

        return errors, self.freqs[np.argmin(errors)]
    
    ##generate decoupling profile interpolated function and function values.   
    def gen_decoupling_profiles(self, Print=False):
        
        self.dec_interpolation = [interp1d(self.cops[0], self.cops[i+1], kind='cubic') for i in range(self.cop_num)]

        self.freqs = np.linspace(self.cops[0][0], self.cops[0][-1],self.interp_points)
        self.predicts =[(1+self.dec_interpolation[i](self.freqs))/2 for i in range(self.cop_num)]

        if Print:
            return self.predicts
        else:
            return None

    ###############################
    #SECTION 3. FITTING UTILITIES # 
    ###############################
    
    '''bounding: [k_abs, fraction, c, j, lwid]'''
    def lineshape_fit(self, x, y, bounding=([0,0,-10,0,0],[40,1,10,40,10])):
        '''
        DEFINITION
        /given/ a list of frequency values and a list of experimental lineshape intensity values, 
        /returns/ the lineshape function's parameters that produce the best fit. 
        
        PARAMETERS
        __________
        x: list
        list of frequency values (Hz). 
        
        y: list, size of x
        lineshape intensity values (arbitrary units). 
        
        bounding: tuple of two lists of length 5, default: ([0,0,-2,0,0],[30,1,2,40,10])
        bounds for arguments of the lineshape function; see the lineshape definition. 
        
        OUTPUT
        ______
        param_best: list
        fit parameters. See description of lineshape for variable definitions.  
        index 0: best fit of the k_abs variable
        index 1: best fit of the fraction variable
        index 3: best fit of the c variable
        index 4: best fit of the j variable
        index 5: best fit of the linewidth variable
        
        '''
        
        param_best, _=scipy.optimize.curve_fit(self.lineshape, x, y,p0=np.add(bounding[0], bounding[1])/2, bounds=bounding)
        return param_best
    
    '''bounding: [k_abs, pyr_fraction, Cb, c, j, lwid]'''
    def lineshape_Cb_fit(self, x, y, prior=None, bounding=([0,0,9,-20,0,0],[10,1,46,20,70,15])):
        '''
        DEFINITION
        /given/ a list of frequency values and a list of experimental lineshape intensity values for the COPs, 
        /returns/ the lineshape function's parameters that produce the best fit. 
        
        PARAMETERS
        __________
        x: list
        list of frequency values (Hz). 
        
        y: list, size of x
        lineshape intensity values of COPs spectra, appended end-to-end (arbitrary units). 
        
        prior: list, default None
        fitting initial parameters for the lineshape_Cb function parameters.
        
        bounding: tuple of two lists of length 5, default: ([0,0,9,-5,0,0],[20,1,46,5,40,10])
        fitting bounds for the lineshape_Cb function parameters; see the lineshape_Cb definition. 
        
        OUTPUT
        ______
        param_best: list
        fit parameters. See description of lineshape_Cb for variable definitions.  
        index 0: best fit of the k_abs variable
        index 1: best fit of the pyr_fraction variable
        index 2: best fit of the Cb variable
        index 3: best fit of the c variable
        index 4: best fit of the j variable
        index 5: best fit of the linewidth variable
        
        '''
        
        if prior==None:
            prior = np.add(bounding[0], bounding[1])/2
        param_best,_ = scipy.optimize.curve_fit(self.lineshape_Cb, x, y, p0 = prior, bounds=bounding)
        return param_best
    
    ######################################
    #SECTION 4. Live lineshape plotting. #
    ######################################
    #for development. 
    
    #generate noisy simulated data. currently takes fixed parameters. 
    def lineshape_gen(self, snr=10, sampling = 150):

        noisewindow = 8
        k_abs_init = 1
        noise_level = k_abs_init/snr
        noise=np.convolve(np.random.normal(0,noise_level,sampling+noisewindow-1), np.ones(noisewindow)/noisewindow, mode='valid')

        w = np.linspace(-100, 100, sampling)
        fraction_init = 0.35

        c_init = 0
        j_init = 20
        lw_init=7
        return lineshape(w, k_abs_init, fraction_init, c_init, j_init, lw_init)+noise
    
    #compute coupling fraction based on Cb. 
    def cop_frac(self, Cbeta, copnum=1):
        
        if copnum > self.copnum:
            print("invalid COP number.")
            return None
        
        a = self.interp1d(self.cops[0], self.cops[copnum], kind='cubic')
        return (1+a(Cbeta))/2
    
    #interactive plot of lineshape depending on Cb.  
    def lineshape_inter_plot(self):
        
        w = np.linspace(-100, 100, 200)
        fraction_init = 0.35
        k_abs_init = 10
        c_init = 0
        j_init = 30
        lw_init=7
        Cb_init=30


        fig, ax = plt.subplots()
        line1, = plt.plot(w, self.lineshape(w, k_abs_init, fraction_init*self.cop_frac(Cb_init), c_init, j_init,lw_init), lw=2, label='grad1')
        line2, = plt.plot(w, self.lineshape(w, k_abs_init, fraction_init*self.cop_frac(Cb_init, copnum=2), c_init, j_init,lw_init), lw=2, label='grad3')
        line3, = plt.plot(w, self.lineshape(w, k_abs_init, fraction_init*self.cop_frac(Cb_init, copnum=3), c_init, j_init,lw_init), lw=2, label='grad4')
        line4, = plt.plot(w, self.lineshape(w, k_abs_init, fraction_init*self.cop_frac(Cb_init, copnum=4), c_init, j_init,lw_init), lw=2, label='grad5')
        line5, = plt.plot(w, self.lineshape(w, k_abs_init, fraction_init*self.cop_frac(Cb_init, copnum=5), c_init, j_init,lw_init), lw=2, label='grad6')
        ax.set_xlabel('frequency (Hz)')
        plt.subplots_adjust(bottom=0.3)
        plt.xlim([-100,100])
        plt.ylim([0, 15])
        plt.legend()

        ax1 = plt.axes([0.2, 0.17, 0.65, 0.03])
        fraction_slider = Slider(
            ax=ax1,
            label='doublet fraction',
            valmin=0,
            valmax=1,
            valinit=fraction_init,
        )


        ax2 = plt.axes([0.2, 0.13, 0.65, 0.03])
        j_slider = Slider(
            ax=ax2,
            label="$J_{ab}$ (Hz)",
            valmin=0,
            valmax=50,
            valinit=j_init,
        )

        ax3 = plt.axes([0.2, 0.09, 0.65, 0.03])
        lw_slider = Slider(
            ax=ax3,
            label="linewidth (Hz)",
            valmin=0,
            valmax=50,
            valinit=lw_init,
        )


        ax4 = plt.axes([0.2, 0.05, 0.65, 0.03])
        Cb_slider = Slider(
            ax=ax4,
            label="C$b$ (ppm)",
            valmin=9,
            valmax=46,
            valinit=Cb_init,
        )




        # The function to be called anytime a slider's value changes
        def update(val):
            line1.set_ydata(self.lineshape(w, k_abs_init, fraction_slider.val*self.cop_frac(Cb_slider.val, copnum=1), c_init, j_slider.val,lw_slider.val))
            line2.set_ydata(self.lineshape(w, k_abs_init, fraction_slider.val*self.cop_frac(Cb_slider.val, copnum=2), c_init, j_slider.val,lw_slider.val))
            line3.set_ydata(self.lineshape(w, k_abs_init, fraction_slider.val*self.cop_frac(Cb_slider.val, copnum=3), c_init, j_slider.val,lw_slider.val))
            line4.set_ydata(self.lineshape(w, k_abs_init, fraction_slider.val*self.cop_frac(Cb_slider.val, copnum=4), c_init, j_slider.val,lw_slider.val))
            line5.set_ydata(self.lineshape(w, k_abs_init, fraction_slider.val*self.cop_frac(Cb_slider.val, copnum=5), c_init, j_slider.val,lw_slider.val))
            fig.canvas.draw_idle()


        # register the update function with each slider
        fraction_slider.on_changed(update)
        j_slider.on_changed(update)
        lw_slider.on_changed(update)
        Cb_slider.on_changed(update)

        # Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
        resetax = plt.axes([0.8, 0.21, 0.1, 0.04])
        button = Button(resetax, 'Reset', hovercolor='0.975')


        def reset(event):
            fraction_slider.reset()
            j_slider.reset()
            lw_slider.reset()
        button.on_clicked(reset)

        plt.show()
        return None
    

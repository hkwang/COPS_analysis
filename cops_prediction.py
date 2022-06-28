#COPs analysis toolbox containing functions for
# 1) amino acid type distributions given CA, CB, and pyruvate
# 2) unsupervised and supervised lineshape matching

from cops_analysis import cops_analyze
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nmrglue as ng


try:
    import pluqin
except:
    print("Check location of folder pluq and pluqin.py.")
    
#read table with pyruvate likelihood distributions.
df = pd.read_excel('./files/pyruvate_likelihood.xlsx')
df.index = df['Residue']


#predict type assignment by CA and CB
def print_probabilities(analyzer, CA, shifts, TMS, pyruvate_on, verbose=False):
    '''
    DEFINITION
    __________
    given a 13C alpha chemical shift and a peak location, this function computes the probability of each amino acid type.
    
    PARAMETERS
    __________
    analyzer: cops_analyze object
    cops_analyze object to conduct lineshape analysis and CB extraction.
    
    CA: float
    13C alpha chemical shift
    
    shifts: 1D array
    chemical shifts of active peak
    
    TMS: boolean
    indicate if chemical shifts are referenced to TMS. If so, converts to DSS. 
    
    pyruvate_on: boolean
    indicate if sample is pyruvate labelled. 
    
    OUTPUTS
    __________
    pandas DataFrame containing amino acid type distribution. 
    
    Also prints out the predicted Cb value. 
    
    '''
    global df
    try:     
        params, error = analyzer.CalcCB(shifts, simple_output=False)
    except: 
        print("No optimum found. invalid or glycine peak selection; or check your initialization.")
    CB_predict = params[2]

    print("\n peak shifts:", shifts)
    print("predicted Cb:", str(round(CB_predict, 2))+"ppm")
    if CA < 50 and error >80:
        prediction = ['G',1]
    elif CA > 55 and CB_predict > 43.5:
        prediction = ['T',1]
    else:
        pluq_out = pluqin.main([[CA+2.6*TMS, CB_predict+2.6*TMS]], 'cc')
        if pluq_out !=None:
            prediction = [[pluq_out[i][0], pluq_out[i][2]/100] for i in range(len(pluq_out))]
        else:
            prediction=['X',1]
    prediction=np.array(prediction).reshape(-1,2)
    if pyruvate_on:
        print('pyruvate fraction:', np.around(params[1],decimals=2))
        prediction = pyruvate_posterior(prediction, params[1],df, verbose=verbose)
    return pd.DataFrame({"type":prediction[:,0], "probability":prediction[:,1]})
    
    
#alter prediction using pyruvate lineshape
def pyruvate_posterior(input_prob, pyruvate_frac, df, verbose=True):
    '''
    DEFINITION
    __________
    given an input amino acid type distribution and a pyruvate fraction, 
    this function computes the posterior distribution conditioned on pyruvate.
    
    PARAMETERS
    __________
    input_prob: 2D array
    prior probability of amino acid type based on Ca and Cb alone.
    format: column 0: amino acid type. column 1: probability of that type. 
    
    pyruvate_frac: float
    pyruvate Cb labeling fraction. 
    
    df: Pandas DataFrame
    dataframe containing probability distributions of pyruvate fraction for each residue type. 
    
    verbose: boolean, default True
    if verbose, prints out prior and conditional distributions. 
    
    OUTPUTS
    __________
    posterior: 2D array
    posterior probability of amino acid type based on Ca, Cb, and pyruvate. 
    format: column 0: amino acid type. column 1: probability of that type. 
    
    '''
    if verbose:
        print("\n prior probability from pluq:\n", pd.DataFrame({"type":input_prob[:,0], "probability":input_prob[:,1]}))
    
    input_prob = np.array(input_prob)
    df['conditional_prob']=gaussian(100*pyruvate_frac, df['mean'].to_numpy(), df['stdev'].to_numpy())
    if verbose:
        print("\n likelihood of each AA given pyruvate fraction:")
        print(df.loc[list(input_prob[:,0])])
    posterior = input_prob[:,1].astype('float64')*df.loc[input_prob[:,0]]['conditional_prob']
    if np.sum(df.loc[input_prob[:,0]]['conditional_prob'])<0.05: 
        print("credence is low. amino acid type may be unlisted.")
    posterior = posterior/np.sum(posterior)
    input_prob[:,1] = np.around(posterior, decimals=2)
    posterior = input_prob
    posterior = posterior[np.argsort(-posterior[:,1].astype('float64'))]

    return posterior

def gaussian(x, mu, sig):
    '''
    DEFINITION
    __________
    this function returns the height at point x of a gaussian with mean mu and standard deviation sig.
    
    '''
    return 1/np.sqrt(2*np.pi*sig)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    

#predict sequential peaks
class int_seq_match():
    def __init__(self, cops_analyzer, cops_mode='HCA', peak_table_dir=None):
        self.cops_mode = cops_mode
        self.cops_analyzer = cops_analyzer
        
        try:
            self.load_peak_table_dir(peak_table_dir)
            self.load_1Ds()
            self.list_mode = True
        except:
            self.list_mode = False
        
    def load_peak_table_dir(self, tb):
        self.tb = tb
        self.shifts_array = self.tb[['CA', 'N', 'HN']].to_numpy(dtype=np.float32)
        self.shifts_array[:,[0,1]]=self.shifts_array[:,[1,0]]
        
    def load_peak_table(self, peak_table_dir):
        '''
        DEFINITION
        __________
        Loads a peak table containing labeled peaks. reads chemical shifts and whether or not peak is sequential or internal.

        PARAMETERS
        __________
        peak_table_dir: string
        location of SPARKY peak table

        OUTPUTS
        __________
        None. This fucntion updates the shifts_array instance and the tb_sequential instance with the peaks from the table.
        '''

        
        if self.cops_mode == 'HCA':
            self.tb = pd.read_fwf(peak_table_dir, infer_nrows=300)
            self.tb = self.tb.rename(columns={'w1':'CA','w2':'HN'})
            self.shifts_array = self.tb[['CA', 'HN']].to_numpy(dtype=np.float32)
            
        elif self.cops_mode == 'HNCA':
            self.tb = pd.read_fwf(peak_table_dir, infer_nrows=300)
            self.tb = self.tb.rename(columns={'w1':'CA','w2':'N','w3':'HN'})
            self.shifts_array = self.tb[['CA', 'N', 'HN']].to_numpy(dtype=np.float32)
            self.shifts_array[:,[0,1]]=self.shifts_array[:,[1,0]]
        
        try:
            self.tb['Assignment']
        except:
            self.tb['is_sequential']=np.append([False], [len(self.tb['Assignment'][i+1]) > len(self.tb['Assignment'][i]) for i in range(len(self.tb)-1)])
    
    def conv_to_sparky_table(self, peak_table_dir):
        return None
    
    
    def load_1Ds(self, C_center=None):
        '''
        DEFINITION
        __________
        Loads 1D 13C slices of peaks in the peak table, centered at C_center (ppm).  
        '''
        ca = self.cops_analyzer
        #extract 1Ds from absolute center
        temp_shifts = np.copy(self.shifts_array)
        if C_center != None:
            if self.cops_mode == 'HNCA':
                temp_shifts[:,1] = C_center
            elif self.cops_mode == 'HCA':
                temp_shifts[:,0] = C_center

        try:
            #this would be faster with a vectorizable nmrglue unit conversion object
            for j in range(len(temp_shifts)):
                slice_1D= np.array([ca.extract1D(temp_shifts[j], ca.cop_dats[i], ca.cop_unit_convs[i],sw=40, normalize=True)[1] 
                                             for i in range(ca.cop_num)]).reshape(-1)
                slice_1D_wide = np.array([ca.extract1D(temp_shifts[j], ca.cop_dats[i], ca.cop_unit_convs[i],sw=150, normalize=True)[1] 
                                             for i in range(ca.cop_num)]).reshape(-1)
                if j == 0:
                    dslice_cop = ([slice_1D])
                    dslice_cop_wide = ([slice_1D_wide])
                else:
                    dslice_cop = dslice_cop+[slice_1D]
                    dslice_cop_wide = dslice_cop_wide + [slice_1D_wide]
                    
            self.slices = np.vstack(dslice_cop) 
            self.slices_wide = np.vstack(dslice_cop_wide)
        except:
            print("Error with loading 1D slices. Check peak list!")
    
    #cleans shifts array of copies of the same peak as the active one
    def remove_picked_copies(self, shifts, peak):
        peaks = np.ones(len(shifts)).reshape([-1,1])@peak.reshape([1,-1])
        diff = np.sum(np.abs(peaks-shifts)**2, axis = 1)
        return shifts[diff>0.1]
        
        
        
    #works for HNCA only
    def pick_slice_peaks(self, peak, blockwidth=15, snr = 5):
        '''
        DEFINITION
        __________
        
        picks peaks in a slice around the active peak.

        PARAMETERS
        __________
        
        peak: 1D array 
        chemical shifts (ppm) of active peak
        
        blockwidth: float
        2*blockwidth is the CA range (in Hz) of the spectrum block in which peaks are picked
        
        snr: float
        size of peak to search for. this number multiplies the standard deviation of intensity of the block 
        

        OUTPUTS
        __________
        None. This function updates the shifts_array instance with nearby peaks.

        '''
        ca = self.cops_analyzer
        
        if self.cops_mode == 'HNCA':
            CA_center = peak[1]
            peak_ind = ca.cop_unit_convs[0][1](CA_center, "ppm")
            bound = ca.hz_to_idx(ca.cop_unit_convs[0][1],blockwidth)
            datablock = ca.cop_dats[0][5:-5,peak_ind-bound:peak_ind+bound,5:-5]
        elif self.cops_mode == 'HCA':
            CA_center = peak[0]
            peak_ind = ca.cop_unit_convs[0][0](CA_center, "ppm")
            bound = ca.hz_to_idx(ca.cop_unit_convs[0][0],blockwidth)
            datablock = ca.cop_dats[0][peak_ind-bound:peak_ind+bound,5:-5]
        
        approx_noise = np.std(datablock)
        results,_, amps = ng.analysis.peakpick.pick(datablock, snr*approx_noise, cluster=False, table=False, est_params=True)
        results_shift = np.array(results)
        if self.cops_mode == 'HNCA':
            results_shift[:,1] += peak_ind - bound
            results_shift[:,0] +=5
            results_shift[:,2] +=5
        elif self.cops_mode == 'HCA':
            results_shift[:,0] += peak_ind - bound
            results_shift[:,1] +=5
            
        #the below is why i would like a vectorizable unit conversion.
        results_shift_ppm = [[ca.cop_unit_convs[0][i].ppm(results_shift[j][i]) for i in range(len(results_shift[0]))] for j in range(len(results_shift))]
        results_shift_ppm = np.array(results_shift_ppm)
        results_shift_ppm = self.remove_picked_copies(results_shift_ppm, peak)
        self.shifts_array = results_shift_ppm
        self.load_1Ds(C_center=CA_center)
        return amps/approx_noise
        
        
    
    def find_best_matches(self, peak, num_best=7, gen_plot=False, label='current peak', snr=None, verbose=False, sequential_mode=True):
        '''
        DEFINITION
        __________
        ranks the highest-correlated peaks to the active peak.

        PARAMETERS
        __________
        peak: 1D array
        chemical shifts (ppm) of active peak
        
        num_best: int, default 7
        number of ranked candidate peaks
        
        gen_plot: Boolean, default False
        if true, returns a plot of the candidate peak lineshapes
        
        label: string, default 'current peak'
        peak label on the plot

        OUTPUTS
        __________
        If gen_plot, outputs a plot of the candidate peak lineshapes. Otherwise, prints table of best matches.

        '''
        
        #if needed, pick peaks and extract amplitudes. #to do: set amplitudes for verbose mode.  
        if self.list_mode: 
            amps =[]
        else:
            amps = self.pick_slice_peaks(peak, snr=snr)

        #compute CA difference likelihood.
        if self.cops_mode == 'HCA':
            CA_diff = peak[0]-self.shifts_array[:,0]
        elif self.cops_mode == 'HNCA':
            CA_diff = peak[1]-self.shifts_array[:,1]

        CA_likelihood = gaussian(CA_diff, 0, 0.05)
        ca = self.cops_analyzer
        data_current = np.array([ca.extract1D(peak, ca.cop_dats[i], ca.cop_unit_convs[i],sw=40, normalize=True)[1] for i in range(ca.cop_num)]).reshape(-1)

        #set indices to correlate to the relevant ones, if in list mode. 
        if self.list_mode:
            if sequential_mode:
                relevant_indices_bool=self.tb['is_sequential']*CA_likelihood>0.01
            else:
                relevant_indices_bool=~self.tb['is_sequential']*CA_likelihood>0.01
        else:
            relevant_indices_bool=range(len(self.shifts_array))

        relevant_indices = relevant_indices_bool.to_numpy().nonzero()[0]
        print(relevant_indices)
        
        #reset num_best if not in list mode and number of peaks picked is small. 
        try:
            assert(len(relevant_indices)>0)
            if len(relevant_indices)<num_best:
                num_best = len(relevant_indices)
        except:
            return 'no peaks picked or read within CA range.', None

        #compute correlation. 
        correlations = np.corrcoef(data_current, self.slices[relevant_indices_bool])[1:,0]
        correlations = np.multiply(correlations, CA_likelihood[relevant_indices_bool])

        #pull out the top num_best matches and their data.     
        index_bestmatch = np.argsort(-correlations, axis = 0)[0:num_best]


        #cleaner reindexing
        if self.list_mode:
            df_index = relevant_indices[index_bestmatch]
        else:
            df_index = index_bestmatch

        correlations_output = correlations[index_bestmatch]
        correlations_output = np.around(correlations_output, decimals=2)


        if self.list_mode:
            labels_output = self.tb['Assignment'][df_index]
            labels_output=np.array(labels_output)
        else:
            labels_output = np.around(self.shifts_array[index_bestmatch], decimals=2).tolist()


        if verbose:


            #amplitude lookup
            if self.list_mode:
                amps = np.zeros(len(df_index))
            else:
                amps = np.around(amps[index_bestmatch], decimals=0)

            #CB calculation
            pred_CB = []
            for i in df_index:
                try:
                    CB = ca.CalcCB(self.shifts_array[i])
                except:
                    CB = 0.0
                pred_CB.append(CB)    


            out_df = pd.DataFrame({"peak": labels_output, "likelihood":np.around(correlations_output, decimals=2), "CB (ppm)": np.around(pred_CB, decimals=2), "amplitude (arb)": amps})
        else:
            out_df = pd.DataFrame({"peak": labels_output, "likelihood":np.around(correlations_output, decimals=2)})

        if gen_plot:
            data_current_wide = np.array([ca.extract1D(peak, ca.cop_dats[i], ca.cop_unit_convs[i],sw=150, normalize=True)[1] for i in range(ca.cop_num)]).reshape(-1)

            cmap = ['b','r','m','y']

            hz,_ = ca.extract1D(peak, ca.cop_dats[0], ca.cop_unit_convs[0], sw=150, normalize=True)
            hz_long = np.array([hz+400*i for i in range(ca.cop_num)]).reshape(-1)

            fig = plt.figure(figsize=(13,6))

            for k in range(ca.cop_num):
                if k ==0:
                    plt.plot(hz+400*k, data_current_wide[0+len(hz)*k:len(hz)*(k+1)],'k', figure = fig, label =label)
                else:
                    plt.plot(hz+400*k, data_current_wide[0+len(hz)*k:len(hz)*(k+1)],'k', figure = fig)

            for i in range(num_best):
                for k in range(ca.cop_num):
                    if k==0:
                        label_str = str(labels_output[i])+', p:'+str(np.around(correlations_output, decimals=2)[i])
                        plt.plot(hz+400*k, self.slices_wide[df_index[i]][0+len(hz)*k:len(hz)*(k+1)]-15*i,cmap[i%4],label = label_str, figure = fig)
                    else:
                        plt.plot(hz+400*k, self.slices_wide[df_index[i]][0+len(hz)*k:len(hz)*(k+1)]-15*i,cmap[i%4], figure = fig)

            fig.legend()

            return out_df, fig
        else:
            return out_df, None

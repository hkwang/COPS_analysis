from cops_analysis import cops_analyze
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


try:
    import pluqin
except:
    print("Check location of folder pluq and pluqin.py.")
 

df = pd.read_excel('./files/pyruvate_likelihood.xlsx')
df.index = df['Residue']


#predict type assignment by CA and CB
def print_probabilities(analyzer, CA, shifts, TMS, pyruvate_on):
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
        pluq_out = pluqin.main([[CA+2*TMS, CB_predict+2*TMS]], 'cc')
        if pluq_out !=None:
            prediction = [[pluq_out[i][0], pluq_out[i][2]/100] for i in range(len(pluq_out))]
        else:
            prediction=['X',1]
    prediction=np.array(prediction).reshape(-1,2)
    if pyruvate_on:
        print('pyruvate fraction:', np.around(params[1],decimals=2))
        prediction = pyruvate_posterior(prediction, params[1],df)
    print(prediction)
    
    
#alter prediction using pyruvate lineshape
def pyruvate_posterior(input_prob, pyruvate_frac, df, verbose=True):
    if verbose:
        print("\n prior probability from pluq:\n", input_prob)
    
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
    return 1/np.sqrt(2*np.pi*sig)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    

#predict sequential peaks
class int_seq_match():
    def __init__(self, peak_table_dir, cops_analyzer, cops_mode='HCA'):
        
        self.cops_mode = cops_mode
        self.cops_analyzer = cops_analyzer
        
        if self.cops_mode == 'HCA':
            self.tb = pd.read_fwf(peak_table_dir, infer_nrows=300)
            self.tb = self.tb.rename(columns={'w1':'CA','w2':'HN'})
            try: 
                self.tb = self.tb.set_index(self.tb['Assignment'])
            except: 
                raise ValueError("SPARKY table assignment column should be titled: 'Assignment'")
            self.tb['is_sequential']=np.append([False], [len(self.tb['Assignment'][i+1]) > len(self.tb['Assignment'][i]) for i in range(len(self.tb)-1)])
            self.shifts_array = self.tb[['CA', 'HN']].to_numpy(dtype=np.float32)
            
        elif self.cops_mode == 'HNCA':
            self.tb = pd.read_fwf(peak_table_dir, infer_nrows=300)
            self.tb = self.tb.rename(columns={'w1':'CA','w2':'N','w3':'HN'})
            try: 
                self.tb = self.tb.set_index(self.tb['Assignment'])
            except: 
                raise ValueError("SPARKY table assignment column should be titled: 'Assignment'")
            self.tb['is_sequential']=np.append([False], [len(self.tb['Assignment'][i+1]) > len(self.tb['Assignment'][i]) for i in range(len(self.tb)-1)])
            self.shifts_array = self.tb[['CA', 'N', 'HN']].to_numpy(dtype=np.float32)
            self.shifts_array[:,[0,1]]=self.shifts_array[:,[1,0]]
        
        self.shifts_array = self.shifts_array[self.tb['is_sequential']]
        self.tb_sequential = self.tb[self.tb['is_sequential']]
        self.load_1Ds()
    
    def load_1Ds(self):
        ca = self.cops_analyzer
        self.seq_slices = np.array([np.array([ca.extract1D(self.shifts_array[j], ca.cop_dats[i], ca.cop_unit_convs[i],sw=90, normalize=True)[1] for i in range(ca.cop_num)]).reshape(-1) for j in range(len(self.shifts_array))])
        self.seq_slices_wide = np.array([np.array([ca.extract1D(self.shifts_array[j], ca.cop_dats[i], ca.cop_unit_convs[i],sw=150, normalize=True)[1] for i in range(ca.cop_num)]).reshape(-1) for j in range(len(self.shifts_array))])
    
    def find_best_matches(self, peak, num_best=3, gen_plot=False, label='current peak'):
        if self.cops_mode == 'HCA':
            CA_diff = peak[0]-self.shifts_array[:,0]
        elif self.cops_mode == 'HNCA':
            CA_diff = peak[1]-self.shifts_array[:,1]
        CA_likelihood = gaussian(CA_diff, 0, 0.05)
        ca = self.cops_analyzer
        data_internal = np.array([ca.extract1D(peak, ca.cop_dats[i], ca.cop_unit_convs[i],sw=90, normalize=True)[1] for i in range(ca.cop_num)]).reshape(-1)
            
        correlations = np.corrcoef(data_internal, self.seq_slices)[1:,0]
        correlations = np.multiply(correlations, CA_likelihood)
        
        index_bestmatch = np.argsort(-correlations, axis = 0)
        
        correlations_output = correlations[index_bestmatch[0:num_best]]/np.sum(correlations[index_bestmatch[0:num_best]])
        labels_output = self.tb_sequential.index[index_bestmatch[0:num_best]]
        print("\n likely sequential peaks:")
        print(np.stack([labels_output, np.around(correlations_output, decimals=2)]).T)
        
        if gen_plot:
            data_internal_wide = np.array([ca.extract1D(peak, ca.cop_dats[i], ca.cop_unit_convs[i],sw=150, normalize=True)[1] for i in range(ca.cop_num)]).reshape(-1)
            
            cmap = ['b','r','m','y']
            
            hz,_ = ca.extract1D(peak, ca.cop_dats[0], ca.cop_unit_convs[0], sw=150, normalize=True)
            hz_long = np.array([hz+400*i for i in range(ca.cop_num)]).reshape(-1)

            fig = plt.figure(figsize=(12,4))

            for k in range(ca.cop_num):
                if k ==0:
                    plt.plot(hz+400*k, data_internal_wide[0+len(hz)*k:len(hz)*(k+1)],'k', figure = fig, label =label)
                else:
                    plt.plot(hz+400*k, data_internal_wide[0+len(hz)*k:len(hz)*(k+1)],'k', figure = fig)

            for i in range(num_best):
                for k in range(ca.cop_num):
                    if k==0:
                        plt.plot(hz+400*k, self.seq_slices_wide[index_bestmatch[i]][0+len(hz)*k:len(hz)*(k+1)]-3*i,cmap[i],label = self.tb_sequential.index[index_bestmatch[i]], figure = fig)
                    else:
                        plt.plot(hz+400*k, self.seq_slices_wide[index_bestmatch[i]][0+len(hz)*k:len(hz)*(k+1)]-3*i,cmap[i], figure = fig)

            fig.legend()
                
            return fig
        else:
            return None
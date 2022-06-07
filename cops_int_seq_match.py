from cops_analysis import cops_analyze
from cops_setup import gaussian
import pandas as pd
import numpy as np

class int_seq_match():
    def __init__(self, peak_table_dir, cops_analyzer, cops_mode='HCA'):
        
        
        self.cops_analyzer = cops_analyzer
        
        if cops_mode == 'HCA':
            self.tb = pd.read_fwf(peak_table_dir, infer_nrows=300)
            self.tb = self.tb.rename(columns={'w1':'CA','w2':'HN'})
            try: 
                self.tb = self.tb.set_index(self.tb['Assignment'])
            except: 
                raise ValueError("SPARKY table assignment column should be titled: 'Assignment'")
            self.tb['is_sequential']=np.append([False], [len(self.tb['Assignment'][i+1]) > len(self.tb['Assignment'][i]) for i in range(len(self.tb)-1)])
            self.shifts_array = self.tb[['CA', 'HN']].to_numpy(dtype=np.float32)
            
        elif cops_mode == 'HNCA':
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
        self.seq_slices = np.array([np.array([ca.extract1D(self.shifts_array[j], ca.cop_dats[i], ca.cop_unit_convs[i],sw=90, normalize=True)[1] for i in range(len(ca.copnames))]).reshape(-1) for j in range(len(self.shifts_array))])
        self.seq_slices_wide = np.array([np.array([ca.extract1D(self.shifts_array[j], ca.cop_dats[i], ca.cop_unit_convs[i],sw=150, normalize=True)[1] for i in range(len(ca.copnames))]).reshape(-1) for j in range(len(self.shifts_array))])
    
    def find_best_matches(self, peak, num_best=3, gen_plot=False):
        CA_diff = peak[0]-self.shifts_array[:,0]
        CA_likelihood = gaussian(CA_diff, 0, 0.05)
        ca = self.cops_analyzer
        data_internal = np.array([ca.extract1D(peak, ca.cop_dats[i], ca.cop_unit_convs[i],sw=90, normalize=True)[1] for i in range(len(ca.copnames))]).reshape(-1)
            
        correlations = np.corrcoef(data_internal, self.seq_slices)[1:,0]
        correlations = np.multiply(correlations, CA_likelihood)
        
        index_bestmatch = np.argsort(-correlations, axis = 0)
        
        correlations_output = correlations[index_bestmatch[0:num_best]]/np.sum(correlations[index_bestmatch[0:num_best]])
        labels_output = self.tb_sequential.index[index_bestmatch[0:num_best]]
        print("likely sequential peaks: \n")
        print(np.stack([labels_output, correlations_output]).T)
        
        if gen_plot:
            data_internal_wide = np.array([ca.extract1D(peak, ca.cop_dats[i], ca.cop_unit_convs[i],sw=150, normalize=True)[1] for i in range(len(ca.copnames))]).reshape(-1)
            
            cmap = ['b','r','m','y']
            
            hz,_ = ca.extract1D(shifts_array[0], ca.cop_dats[2], ca.cop_unit_convs[2], sw=150, normalize=True)
            hz_long = np.array([hz+400*i for i in range(len(ca.copnames))]).reshape(-1)

            fig = plt.figure(figsize=(15,5))

            for k in range(len(ca.copnames)):
                if k ==0:
                    plt.plot(hz+400*k, data_internal_wide[0+len(hz)*k:len(hz)*(k+1)],'k', figure = fig, label = "current peak")
                else:
                    plt.plot(hz+400*k, data_internal_wide[0+len(hz)*k:len(hz)*(k+1)],'k', figure = fig)

            for i in range(num_best):
                for k in range(len(ca.copnames)):
                    if k==0:
                        plt.plot(hz+400*k, data_sequential[index_bestmatch[i]][0+len(hz)*k:len(hz)*(k+1)]-3*i,cmap[i],label = self.tb_sequential.index[index_bestmatch[i]], figure = fig)
                    else:
                        plt.plot(hz+400*k, data_sequential[index_bestmatch[i]][0+len(hz)*k:len(hz)*(k+1)]-3*i,cmap[i], figure = fig)

            plt.legend(figure = fig)
                
            return fig
        else:
            return None
            
        
        
        
    
    
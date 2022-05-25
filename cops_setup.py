from cops_analysis import cops_analyze
import numpy as np
import pandas as pd


try:
    import pluqin
except:
    print("Check location of folder pluq and pluqin.py.")
 

df = pd.read_excel('./files/pyruvate_likelihood.xlsx')
df.index = df['Residue']



def print_probabilities(analyzer, CA, shifts, TMS, pyruvate_on):
    global df
    try:     
        params, error = analyzer.CalcCB(shifts, simple_output=False)
    except: 
        print("No optimum found. invalid or glycine peak selection; or check your initialization.")
    CB_predict = params[2]
    print('\n')
    print("peak shifts:", shifts)
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
        print('pyruvate fraction:', round(params[1],2))
        prediction = pyruvate_posterior(prediction, params[1],df)
    print(prediction)
    
    
def pyruvate_posterior(input_prob, pyruvate_frac, df, verbose=True):
    if verbose:
        print("prior probability from pluq:\n", input_prob)
    
    input_prob = np.array(input_prob)
    df['conditional_prob']=gaussian(100*pyruvate_frac, df['mean'].to_numpy(), df['stdev'].to_numpy())
    if verbose:
        print("likelihood of each AA given pyruvate fraction:")
        print(df.loc[list(input_prob[:,0])])
    posterior = input_prob[:,1].astype('float64')*df.loc[input_prob[:,0]]['conditional_prob']
    if np.sum(df.loc[input_prob[:,0]]['conditional_prob'])<0.05: 
        print("credence is low. amino acid type may be unlisted.")
    posterior = posterior/np.sum(posterior)
    input_prob[:,1] = posterior.round(2)
    posterior = input_prob
    posterior = posterior[np.argsort(-posterior[:,1].astype('float64'))]
    return posterior

def gaussian(x, mu, sig):
    return 1/np.sqrt(2*np.pi*sig)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    

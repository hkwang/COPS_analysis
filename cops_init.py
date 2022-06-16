#A GUI for interactively analyzing COPs spectra.

import os
from cops_analysis import cops_analyze
import numpy as np
from cops_prediction import *

try:
    import pluqin
except:
    print("Check location of folder pluq and pluqin.py.")
try: 
    import __main__
    s=__main__.main_session
except:
    pass


import tkinter.font
from tkinter import *
from tkinter.filedialog import askopenfilename, askdirectory

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)


def main():
    #run the TKinter window's main loop.
    program = COPS_GUI()
    root = Tk()
    program.main_frame(root)
    root.mainloop()
    

class InitWindow(Toplevel):
     
    def __init__(self, main, master = None):
         
        super().__init__(master = master)
        self.main = main
        self.title("Initialization Window")
        Label(self, text ="Initialization", font=self.main.font_title).grid(row=1, column=2)
        
        self.create_widgets()
        
    def create_widgets(self):
        
        ##input directories
        label1 = Label(self, text="peak list",font=self.main.font_bold)
        label1.grid(row=2,column=2)    
        self.ent1=Entry(self,font=self.main.font_text, width=40)
        self.ent1.grid(row=4,column=1, columnspan=3)
        b1=Button(self,text="Select",font=self.main.font_text,command=self.main.browsedir1)
        b1.grid(row=4,column=4)

        label2 = Label(self, text='no-COPs spectrum',font=self.main.font_bold)
        label2.grid(row=5,column=2)
        self.main.ent2=Entry(self,font=self.main.font_text, width=40)
        self.main.ent2.grid(row=6,column=1, columnspan=3)
        b2=Button(self,text="Select",font=self.main.font_text, command=self.main.browsefile)
        b2.grid(row=6,column=4)

        label4 = Label(self, text='GRADCOP 1 spectrum',font=self.main.font_bold)
        label4.grid(row=7,column=2) 
        self.ent3=Entry(self,font=self.main.font_text, width=40)
        self.ent3.grid(row=8,column=1, columnspan=3)
        b3=Button(self,text="Select",font=self.main.font_text, command=self.main.browsefile1)
        b3.grid(row=8,column=4)

        label5 = Label(self, text='GRADCOP 3 spectrum',font=self.main.font_bold)
        label5.grid(row=9,column=2) 
        self.ent4=Entry(self,font=self.main.font_text, width=40)
        self.ent4.grid(row=10,column=1, columnspan=3)
        b4=Button(self, text="Select",font=self.main.font_text, command=self.main.browsefile2)
        b4.grid(row=10,column=4)

        label6 = Label(self, text='GRADCOP 4 spectrum',font=self.main.font_bold)
        label6.grid(row=11,column=2) 
        self.ent5=Entry(self,font=self.main.font_text, width=40)
        self.ent5.grid(row=12,column=1, columnspan=3)
        b5=Button(self,text="Select",font=self.main.font_text, command=self.main.browsefile3)
        b5.grid(row=12,column=4)

        label7 = Label(self, text='GRADCOP 5 spectrum',font=self.main.font_bold)
        label7.grid(row=13,column=2) 
        self.ent6=Entry(self,font=self.main.font_text, width=40)
        self.ent6.grid(row=14,column=1, columnspan=3)
        b6=Button(self,text="Select",font=self.main.font_text, command=self.main.browsefile4)
        b6.grid(row=14,column=4)

        label8 = Label(self, text='GRADCOP 6 spectrum',font=self.main.font_bold)
        label8.grid(row=15,column=2) 
        self.ent7=Entry(self,font=self.main.font_text, width=40)
        self.ent7.grid(row=16,column=1, columnspan=3)
        b7=Button(self,text="Select",font=self.main.font_text, command=self.main.browsefile5)
        b7.grid(row=16,column=4)

        ##COPs experiment options       
        label3=Label(self, text='COPS experiment',font=self.main.font_text)
        label3.grid(row=17, column=1)
        
        options_list = ["HNCA", "HCA", 'HN(co)CA']
        self.question_menu = OptionMenu(self, self.main.cops_mode, *options_list)
        self.question_menu.grid(row=17, column=2)

        ##check boxes for options
        self.c1 = Checkbutton(self, text='pyruvate labelled',variable=self.main.pyr_on, onvalue=True, offvalue=False)
        self.c1.grid(row=19,column=2)
        
        self.c2 = Checkbutton(self, text='TMS referenced',variable=self.main.TMS, onvalue=True, offvalue=False)
        self.c2.grid(row=19,column=1)
        
        b11 = Button(self, text='save state',font=self.main.font_text, command=self.main.save_state)
        b11.grid(row=19, column=3)
        
        b12 = Button(self, text='load state',font=self.main.font_text, command=self.main.load_state)
        b12.grid(row=19,column=4)
        
        b8=Button(self,text="Initialize",font=self.main.font_text, command=self.main.init_analyzer)
        b8.grid(row=20, column=1)
        
        b9 =Button(self, text='clear directories', font=self.main.font_text,command=self.main.clear_dirs)
        b9.grid(row=18,column=3)
        
        btn1 = Button(self, text ="Exit", font=self.main.font_text, command = self.destroy)
        btn1.grid(row=19,column=5)

class PredOptionWindow(Toplevel):
     
    def __init__(self, main, master = None):
         
        super().__init__(master = master)
        self.main = main
        self.title("Options Window")
        
        self.create_widgets()
        
    def create_widgets(self):
        
        c3 = Checkbutton(self, text='prediction plot on',variable=self.main.predict_plot, onvalue=True, offvalue=False)
        c3.grid(row=2,column=1)
        
        Label(self, text='peak picking SNR:').grid(row=3, column=1)
        Entry(self, textvariable=self.main.SNR).grid(row=3, column=2)
        
        c4 = Checkbutton(self, text='calculation verbose',variable=self.main.verbose, onvalue=True, offvalue=False)
        c4.grid(row=4,column=1)
         

class COPS_GUI:
    #A GUI for interactively analyzing COPs spectra.
    
    def __init__(self):
        pass
    
    
    def main_frame(self, root):
        self.root = root
        self.root.title('GRADCOPs analysis')
        
        #Initialize TKinter window and format 
        self.font_title = tkinter.font.Font(family = "Arial", 
                                         size = 18)
        self.font_bold = tkinter.font.Font(family = "Arial", 
                                         size = 14)
        self.font_text = tkinter.font.Font(family = "Arial", 
                                         size = 12)
        
        #Initialize objects for aiding COPs analysis
        self.analyzer = None
        self.matcher = None
        
        #instances of COPs filenames
        self.copnames = []
        self.copnums = []
        self.filetypes = (("all files","*.*"),("sparky files",'*.ucsf'),("NMRpipe files","*.ft3"))
        self.tabledir = None
        
        #initialize Tkinter variables for specifying the experiment type
        self.pyr_on=BooleanVar(self.root)
        self.TMS = BooleanVar(self.root)
        self.predict_plot=BooleanVar(self.root)
        self.verbose = BooleanVar(self.root)
        
        self.cops_mode = StringVar(self.root)
        self.cops_mode.set("Select an Option")
        
        self.SNR = StringVar(self.root)
        self.SNR.set('5')
        
        self.create_widgets()

    def create_widgets(self):
        Label(self.root, text="GRADCOPs analysis", font=self.font_title).grid(row=1, column=2)
        
        ##buttons for analysis and control
        b8=Button(self.root,text="Initialize",font=self.font_text)
        b8.bind("<Button>", lambda e: InitWindow(self, self.root))
        b8.grid(row=20,column=1)
        
        b9=Button(self.root, text="Prediction Options",font=self.font_text)
        b9.bind("<Button>", lambda e: PredOptionWindow(self, self.root))
        b9.grid(row=20, column=4)

        b7=Button(self.root,text="calculate",font=self.font_text, command=self.calculate)
        b7.grid(row=20,column=2)
        
        b13=Button(self.root,text="predict",font=self.font_text, command=self.predict)
        b13.grid(row=20,column=3)

        b12=Button(self.root,text="append",font=self.font_text)#, command=calculate)
        b12.grid(row=18,column=4)

        b10=Button(self.root,text="save",font=self.font_text)#, command=calculate)
        b10.grid(row=18,column=5)
        
        btn1 = Button(self.root, text ="Exit", font=self.font_text, command = self.root.destroy)
        btn1.grid(row=19,column=5)
        
    

        
    ##functions called on file selection button press
    ##main directory
    def browsedir1(self):
        self.tabledir = askopenfilename(title="Select a file", filetypes=(("SPARKY lists", "*.list"),("all files","*.*")))
        if self.tabledir != '':
            self.ent1.insert(END, self.tabledir)     

    ##nocop spectrum
    def browsefile(self):
        nocop_dir = askopenfilename(title="Select a file", filetypes=self.filetypes)
        if nocop_dir != '':
            self.copnames.append(nocop_dir)
            self.copnums.append(0)
            self.ent2.insert(END, nocop_dir)       

    ##COP 1
    def browsefile1(self):
        cop1_dir = askopenfilename(title="Select a file", filetypes=self.filetypes)
        if cop1_dir != '':
            self.copnames.append(cop1_dir)
            self.copnums.append(1)
            self.ent3.insert(END, cop1_dir) 

    ##COP 3
    def browsefile2(self):
        cop3_dir = askopenfilename(title="Select a file", filetypes=self.filetypes)
        if cop3_dir != '':
            self.copnames.append(cop3_dir)
            self.copnums.append(2)
            self.ent4.insert(END, cop3_dir) 

    ##COP 4
    def browsefile3(self):
        cop4_dir = askopenfilename(title="Select a file", filetypes=self.filetypes)
        if cop4_dir != '':
            self.copnames.append(cop4_dir)
            self.copnums.append(3)
            self.ent5.insert(END, cop4_dir) 

    ##COP 5
    def browsefile4(self):
        cop5_dir = askopenfilename(title="Select a file", filetypes=self.filetypes)
        if cop5_dir != '':
            self.copnames.append(cop5_dir)
            self.copnums.append(4)
            self.ent6.insert(END, cop5_dir) 

    ##COP 6
    def browsefile5(self):
        cop6_dir = askopenfilename(title="Select a file", filetypes=self.filetypes)
        if cop6_dir != '':
            self.copnames.append(cop6_dir)
            self.copnums.append(5)
            self.ent7.insert(END, cop6_dir) 


    #button-press function for analyzing COPs lineshape and predicting amino acid type

    def calculate(self):
        peaks = s.selected_peaks()
        
        #filter out glycines
        for peak in peaks:
            CA = 0
            for i in peak.frequency:
                if 46<i<70:
                    CA = i
            
            if self.cops_mode.get()=='HNCA':
                freqs = np.array(peak.frequency)[[1,0,2]]
            else:
                freqs = peak.frequency
            print("___________ \n"+"New calculation")
            print(print_probabilities(self.analyzer,CA, freqs, self.TMS.get(), self.pyr_on.get(), verbose=self.verbose.get()))

    #button-press function for displaying candidate internal-sequential matches
    def predict(self):
        if self.tabledir != None:
            self.matcher = int_seq_match(self.analyzer, cops_mode=self.cops_mode.get(), peak_table_dir=self.tabledir)
        peaks = s.selected_peaks()
        try:
            SNR = float(self.SNR.get())
        except:
            print("SNR is not a number.")
        for peak in peaks:
            assign = peak.assignment
            if self.cops_mode.get()=='HNCA':
                freqs = np.array(peak.frequency)[[1,0,2]]
            else:
                freqs = peak.frequency
            print("___________ \n"+"New prediction")
            result = self.matcher.find_best_matches(freqs, gen_plot=self.predict_plot.get(), label=assign, snr=SNR)
            try:
                
                self.plot(result)
            except:
                pass 
    
    #plot lineshapes
    def plot(self, plotfigure):
        canvas = FigureCanvasTkAgg(plotfigure,
                                   master = self.root)  
        canvas.draw()
        canvas.get_tk_widget().grid(row=23, column=1, columnspan=5)
        toolbar = NavigationToolbar2Tk(canvas,
                                       self.root)
        toolbar.update()
        canvas.get_tk_widget().grid(row=24, column=3)          
            
        
    #initialize both objects
    def init_analyzer(self):
        if self.copnums[0]==0:
            self.copnums.pop(0)
        try:
            self.analyzer = cops_analyze(self.copnames, mode=self.cops_mode.get(), pyruvate_on=self.pyr_on.get(), cop_num=self.copnums)
            self.matcher = int_seq_match(self.analyzer, cops_mode=self.cops_mode.get(), peak_table_dir=self.tabledir)
            print('Initialization complete.')
        except:
            print('Initialization incomplete.')
    
    #clear directories
    def clear_dirs(self):        
        #collect all cops spectra
        self.copnames = []
        self.copnums = []
        
        #self.ent1.delete(0,END)
        self.ent2.delete(0,END)
        self.ent3.delete(0,END)
        self.ent4.delete(0,END)
        self.ent5.delete(0,END)
        self.ent6.delete(0,END)
        self.ent7.delete(0,END)
      
    #save to dictionary. Currently, there are no options to change the filename. 
    def save_state(self):
        save_dict={'copnames':self.copnames, 'copnums':self.copnums,
                   'TMS':self.TMS.get(), 'pyr_on':self.pyr_on.get(), 
                   'cops_mode':self.cops_mode.get(), 'tabledir':self.tabledir}
        np.save('COPS_savestate_1.npy',save_dict)
        print('state save complete.')
    
    #load from dictionary.
    def load_state(self):
        try:
            load_dict = np.load('COPS_savestate_1.npy', allow_pickle=True).item()
            self.copnames = load_dict.get('copnames'); self.copnums=load_dict.get('copnums');
            self.TMS.set(load_dict.get('TMS'));self.pyr_on.set(load_dict.get('pyr_on'));
            self.cops_mode.set(load_dict.get('cops_mode')); self.tabledir=load_dict.get('tabledir');
            #self.ent1.insert(END, self.tabledir)
            self.ent2.insert(END, 'COPs loaded')
            self.init_analyzer()
        except:
            print('Save file missing or renamed.')


        
                           
                        
                   

main()
import os
from cops_analysis import cops_analyze
import numpy as np
from cops_setup import *

try:
    import pluqin
except:
    print("Check location of folder pluq and pluqin.py.")

import __main__
s=__main__.main_session


import tkinter.font
from tkinter import *
from tkinter.filedialog import askopenfilename, askdirectory

def main():
    program = COPS_GUI()
    program.root.mainloop()
    

class COPS_GUI:
    def __init__(self):
        self.root = Tk()
        self.root.title('GRADCOPs analysis')
        self.font_title = tkinter.font.Font(family = "Arial", 
                                         size = 18)
        self.font_bold = tkinter.font.Font(family = "Arial", 
                                         size = 14)
        self.font_text = tkinter.font.Font(family = "Arial", 
                                         size = 12)
        self.analyzer = None
        
        #collect all cops spectra
        self.copnames = []
        self.copnums = []
        self.filetypes = (("all files","*.*"),("sparky files",'*.ucsf'),("NMRpipe files","*.ft3"))
        self.dirname = os.getcwd()
        print(self.dirname)
        
        self.pyr_on=BooleanVar(self.root)
        self.TMS = BooleanVar(self.root)
        
        self.cops_mode = StringVar(self.root)
        self.cops_mode.set("Select an Option")
        
        self.create_widgets()

    def create_widgets(self):
        Label(self.root, text="GRADCOPs analysis", font=self.font_title).grid(row=1, column=2)

        '''
        label1 = Label(self.root, text="directory",font=self.font_bold)
        label1.grid(row=2,column=2)    
        self.ent1=Entry(self.root,font=self.font_text, width=40)
        self.ent1.grid(row=4,column=2)
        b1=Button(self.root,text="Select",font=self.font_text,command=self.browsedir1)
        b1.grid(row=4,column=4)
        '''

        label2 = Label(self.root, text='no-COPs spectrum',font=self.font_bold)
        label2.grid(row=5,column=2) 
        self.ent2=Entry(self.root,font=self.font_text, width=40)
        self.ent2.grid(row=6,column=2)
        b2=Button(self.root,text="Select",font=self.font_text, command=self.browsefile)
        b2.grid(row=6,column=4)

        label4 = Label(self.root, text='GRADCOP 1 spectrum',font=self.font_bold)
        label4.grid(row=7,column=2) 
        self.ent3=Entry(self.root,font=self.font_text, width=40)
        self.ent3.grid(row=8,column=2)
        b3=Button(self.root,text="Select",font=self.font_text, command=self.browsefile1)
        b3.grid(row=8,column=4)

        label5 = Label(self.root, text='GRADCOP 3 spectrum',font=self.font_bold)
        label5.grid(row=9,column=2) 
        self.ent4=Entry(self.root,font=self.font_text, width=40)
        self.ent4.grid(row=10,column=2)
        b4=Button(self.root,text="Select",font=self.font_text, command=self.browsefile2)
        b4.grid(row=10,column=4)

        label6 = Label(self.root, text='GRADCOP 4 spectrum',font=self.font_bold)
        label6.grid(row=11,column=2) 
        self.ent5=Entry(self.root,font=self.font_text, width=40)
        self.ent5.grid(row=12,column=2)
        b5=Button(self.root,text="Select",font=self.font_text, command=self.browsefile3)
        b5.grid(row=12,column=4)

        label7 = Label(self.root, text='GRADCOP 5 spectrum',font=self.font_bold)
        label7.grid(row=13,column=2) 
        self.ent6=Entry(self.root,font=self.font_text, width=40)
        self.ent6.grid(row=14,column=2)
        b6=Button(self.root,text="Select",font=self.font_text, command=self.browsefile4)
        b6.grid(row=14,column=4)

        label8 = Label(self.root, text='GRADCOP 6 spectrum',font=self.font_bold)
        label8.grid(row=15,column=2) 
        self.ent7=Entry(self.root,font=self.font_text, width=40)
        self.ent7.grid(row=16,column=2)
        b7=Button(self.root,text="Select",font=self.font_text, command=self.browsefile5)
        b7.grid(row=16,column=4)

        ##COPs experimself.ent type       
        label3=Label(self.root, text='COPS experiment',font=self.font_text)
        label3.grid(row=17, column=1)

        options_list = ["HNCA", "HCA", 'HN(co)CA']
        self.question_menu = OptionMenu(self.root, self.cops_mode, *options_list)
        self.question_menu.grid(row=17, column=2)

        ##pyruvate on or off?
        self.c1 = Checkbutton(self.root, text='pyruvate labelled',variable=self.pyr_on, onvalue=True, offvalue=False)
        self.c1.grid(row=19,column=2)
        
        self.c1 = Checkbutton(self.root, text='TMS referenced',variable=self.TMS, onvalue=True, offvalue=False)
        self.c1.grid(row=19,column=1)

        b8=Button(self.root,text="Initialize",font=self.font_text, command=self.init_analyzer)
        b8.grid(row=20,column=1)

        b7=Button(self.root,text="calculate",font=self.font_text, command=self.calculate)
        b7.grid(row=20,column=2)

        b8=Button(self.root,text="append",font=self.font_text)#, command=calculate)
        b8.grid(row=18,column=4)
        
        b9 =Button(self.root, text='clear directories', font=self.font_text,command=self.clear_dirs)
        b9.grid(row=18,column=3)

        b10=Button(self.root,text="save",font=self.font_text)#, command=calculate)
        b10.grid(row=18,column=5)
        
        btn1 = Button(self.root, text ="Exit", font=self.font_text, command = self.root.destroy)
        btn1.grid(row=19,column=5)
        
        b11 = Button(self.root, text='save state',font=self.font_text, command=self.save_state)
        b11.grid(row=19, column=3)
        
        b12 = Button(self.root, text='load state',font=self.font_text, command=self.load_state)
        b12.grid(row=19,column=4)


    ##main directory
    def browsedir1(self):
        self.dirname = askdirectory()
        self.ent1.insert(END, self.dirname)
        os.chdir(self.dirname)

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


    #compute 

    def calculate(self):
        peaks = s.selected_peaks()
        for peak in peaks:
            CA = 0
            for i in peak.frequency:
                if 45<i<70:
                    CA = i
            
            if self.cops_mode.get()=='HNCA':
                freqs = np.array(peak.frequency)[[1,0,2]]
            else:
                freqs = peak.frequency

            print_probabilities(self.analyzer,CA, freqs, self.TMS.get(), self.pyr_on.get())


    #initialize
    def init_analyzer(self):
        os.chdir(self.dirname)

        if self.copnums[0]==0:
            self.copnums.pop(0)

        try:
            self.analyzer = cops_analyze(self.copnames, mode=self.cops_mode.get(), pyruvate_on=self.pyr_on.get(), cop_num=self.copnums)
            print('Initialization complete.')
        except:
            print('Initialization incomplete.')
    
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
        
    def save_state(self):
        save_dict={'copnames':self.copnames, 'copnums':self.copnums,
                   'TMS':self.TMS.get(), 'pyr_on':self.pyr_on.get(), 
                   'cops_mode':self.cops_mode.get()}
        np.save('save_1.npy',save_dict)
        print('state save complete.')
        
    def load_state(self):
        try:
            load_dict = np.load('save_1.npy', allow_pickle=True).item()
            self.copnames = load_dict.get('copnames'); self.copnums=load_dict.get('copnums');
            self.TMS.set(load_dict.get('TMS'));self.pyr_on.set(load_dict.get('pyr_on'));
            self.cops_mode.set(load_dict.get('cops_mode')); self.dirname=load_dict.get('dir');
            #self.ent1.insert(END, self.dirname)
            self.ent2.insert(END, 'COPs loaded')
            self.init_analyzer()
        except:
            print('Select directory, or save file missing.')


        
                           
                        
                   

main()
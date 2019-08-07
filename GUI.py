import tkinter as tk
import sys
from os import listdir
from os.path import isfile, join
import threading


class GUI:

    def __init__(self, start_routine):
        self.top = tk.Tk()
        self.topframe = tk.Frame(self.top)
        self.topframe.pack()
        self.leftframe = tk.Frame(self.topframe)
        self.leftframe.pack(side=tk.LEFT)
        self.rightframe = tk.Frame(self.topframe)
        self.rightframe.pack(side=tk.RIGHT)
        self.bottomframe = tk.Frame(self.top)
        self.bottomframe.pack(side=tk.BOTTOM)
        self.output = tk.Text(self.bottomframe)
        self.output.insert(tk.INSERT, "Starting...\n")
        self.start_routine = start_routine
        self.output.pack()


    def start(self):
        var1 = tk.IntVar()
        var2 = tk.IntVar()
        var1.set(0)
        var2.set(0)
        radiobuttons = []
        annot_radiobuts = []
        onlyfiles = [f for f in listdir("./data") if isfile(join("./data", f))]
        data_files = []
        annot_files = []
        for file in onlyfiles:
            if file[-4:] == "3end":
                radiobuttons.append(tk.Radiobutton(self.leftframe, text=file[:-12], variable=var1, value=len(radiobuttons)))
                radiobuttons[-1].pack(anchor=tk.W)
                data_files.append("data/" + file)
            if file[-5:] == "annot":
                annot_radiobuts.append(tk.Radiobutton(self.leftframe, text=file[:-6], variable=var2, value=len(annot_radiobuts)))
                annot_radiobuts[-1].pack(anchor=tk.W)
                annot_files.append("data/" + file)
        def but1():
            DATA_FILE = data_files[var1.get()]
            ANNOT_FILE = annot_files[var2.get()]
            # thread = threading.Thread(target=self.start_routine, args=(self, DATA_FILE, ANNOT_FILE))
            # thread.start()
            # thread.join()
            self.start_routine(self, DATA_FILE, ANNOT_FILE)
            self.top.mainloop()
            # self.top.destroy()
        b1 = tk.Button(self.leftframe, text="Start", command=but1)
        b1.pack(side=tk.BOTTOM)
        CheckVar1 = tk.IntVar()
        CheckVar2 = tk.IntVar()
        entries = []
        def c1():
            if CheckVar2.get() == 1:
                entries.append(tk.Entry(self.rightframe, bd=5))
                entries.append(tk.Label(self.rightframe, text="Enter full or relative path \nwith the filename: "))
                entries[1].pack()
                entries[0].pack(side=tk.BOTTOM)
            else:
                entries[0].destroy()
                entries[1].destroy()
                entries.pop(0)
                entries.pop(0)
        c_pca = tk.Checkbutton(self.rightframe, text="Show PCA", variable=CheckVar1, onvalue=1, offvalue=0, height=5, width=20)
        c_to_pdf = tk.Checkbutton(self.rightframe, text="Export to PDF", variable=CheckVar2, onvalue=1, offvalue=0, height=5, width=20,
                                  command=c1)
        c_pca.pack()
        c_to_pdf.pack()
        self.top.mainloop()

    def write_to_output(self, msg = "None"):
        self.output.insert(tk.INSERT, msg)
        self.top.update()

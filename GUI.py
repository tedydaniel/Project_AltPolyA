import tkinter as tk
import sys
from os import listdir
from os.path import isfile, join

def make_gui():
    top = tk.Tk()
    topframe = tk.Frame(top)
    topframe.pack()
    leftframe = tk.Frame(topframe)
    leftframe.pack(side=tk.LEFT)
    rightframe = tk.Frame(topframe)
    rightframe.pack(side=tk.RIGHT)
    bottomframe = tk.Frame(top)
    bottomframe.pack(side=tk.BOTTOM)
    output = tk.Text(bottomframe)
    output.insert(tk.INSERT, "Hello")
    output.pack()

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
            radiobuttons.append(tk.Radiobutton(leftframe, text=file[:-12], variable=var1, value=len(radiobuttons)))
            radiobuttons[-1].pack(anchor=tk.W)
            data_files.append("data/" + file)
        if file[-5:] == "annot":
            annot_radiobuts.append(tk.Radiobutton(leftframe, text=file[:-6], variable=var2, value=len(annot_radiobuts)))
            annot_radiobuts[-1].pack(anchor=tk.W)
            annot_files.append("data/" + file)
    def but1():
        global DATA_FILE
        global ANNOT_FILE
        DATA_FILE = data_files[var1.get()]
        ANNOT_FILE = annot_files[var2.get()]
        top.destroy()
    b1 = tk.Button(leftframe, text="Start", command=but1)
    b1.pack(side=tk.BOTTOM)
    CheckVar1 = tk.IntVar()
    CheckVar2 = tk.IntVar()
    entries = []
    def c1():
        if CheckVar2.get() == 1:
            entries.append(tk.Entry(rightframe, bd=5))
            entries.append(tk.Label(rightframe, text="Enter full or relative path \nwith the filename: "))
            entries[1].pack()
            entries[0].pack(side=tk.BOTTOM)
        else:
            entries[0].destroy()
            entries[1].destroy()
            entries.pop(0)
            entries.pop(0)
    c_pca = tk.Checkbutton(rightframe, text="Show PCA", variable=CheckVar1, onvalue=1, offvalue=0, height=5, width=20)
    c_to_pdf = tk.Checkbutton(rightframe, text="Export to PDF", variable=CheckVar2, onvalue=1, offvalue=0, height=5, width=20,
                              command=c1)
    c_pca.pack()
    c_to_pdf.pack()
    top.mainloop()
#========================
# imports
#========================

import tkinter as tk
from tkinter import ttk


# Create Instance
win = tk.Tk()
win.size()
# Title
win.title("MyFirstGUI")

#makes it resizable
# win.resizable(True,False)

#defining label text size
l_fontsize=18
l_font="Courier"

#adding a label
a_label = ttk.Label(win, text="a Label")
a_label.configure(font=(l_font, l_fontsize))
a_label.grid(column=0, row=1)

#Buttonclickeventfunction
def click_me():
    action.configure(text="***I have been clicked***")
    a_label.configure(foreground='red')
    a_label.configure(text='a Red Label')

action = ttk.Button(win,text='click me',command = click_me)
action.grid(column=0, row=2)


#====================================


name_label = ttk.Label(text='Enter a name')
name_label.configure(font=(l_font, l_fontsize))
name_label.grid(column=1, row=0)

#adding a textbox
name = tk.StringVar()
name_entered = ttk.Entry(win, width=12, textvariable = name)
name_entered.grid(column = 1 , row = 1)

#Buttonclickeventfunction
def click_me2():
    actiontext.configure(text="Hello"+name.get())
    actiontext.configure(state='disabled')

actiontext = ttk.Button(win,text='click me',command = click_me2)
actiontext.grid(column=1, row=2)



#=====================================
#making a Combobox
l_choosebox = ttk.Label(text='Choose a EXPNO:')
l_choosebox.config(font=(l_font, l_fontsize))
l_choosebox.grid(column=2, row=0)


number = tk.StringVar()
number_chosen = ttk.Combobox(win, width = 12, textvariable=number, state='readonly') #readonly means that it can onlyread the given values and no manual text entry
number_chosen['values'] = ('not selected',1,2,3,42,100)
number_chosen.grid(column=2, row=1)
number_chosen.current(0)

#======================================
colors = ['red', 'green', 'blue']

#making radboxes
def radCall():
    radsel=radVar.get()
    if radsel == 0: win.configure(background = colors[0])
    elif radsel == 1: win.configure(background = colors[1])
    elif radsel == 2: win.configure(background= colors[2])

radVar=tk.IntVar()
radVar.set(99)

for col in range(len(colors)):
    curRad = tk.Radiobutton(win, text=colors[col], variable= radVar,
                            value=col, command= radCall)
    curRad.grid(column=3, row=col, sticky=tk.W)

#======================================
#matplotlib button
import matplotlib.pyplot as plt

x = [1,2,3,4,5,6,7]
y = [2,5,6,1,6,2,3]
def plotgraph():
    plt.scatter(x,y)
    plt.show()


matbutton = ttk.Button(win, text="plot a graph", command=plotgraph)
matbutton.grid(column=1, row=1)

name_entered.focus()
#=========================
# Start GUI
#=========================

win.mainloop()
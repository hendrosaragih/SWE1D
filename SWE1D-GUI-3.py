#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 19:13:11 2020

@author: macbookpro
"""

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import Frame,Label,Entry,Button
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from call_SWE1D import *
from initial_Condition import *

class Window(Frame):

	def __init__(self, master = None):
		Frame.__init__(self, master)
		self.master = master
		self.init_window()

	def Run(self):
		self.eta_list = list()
		self.xleft   = 0.
		self.xright  = 1.
		# domain
		self.Nx      = 100
		self.dx      = 0.01
		self.x  = [0]*(self.Nx+1)
		self.x = np.linspace(self.xleft, self.xright, self.Nx + 1)
		self.eta_list = call_SWE1D(self.xleft, self.xright, self.Nx, self.dx)
		self.eta_len = len(self.eta_list)
		
		def animate(i):
			self.y = self.eta_list[i][:]
			self.line.set_data(self.x,self.y)  # update the data
			return self.line,

		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(1,1,1)
		self.ax = plt.axes(xlim=(0, 1), ylim=(-0.002, 0.002))
		self.line, = self.ax.plot([], [], lw=2)
		self.ax = plt.grid(True)

		self.canvas = FigureCanvasTkAgg(self.fig, master=self)
		self.canvas.get_tk_widget().grid(column=0,row=4)

		self.ani = animation.FuncAnimation(self.fig, animate, np.arange(1, self.eta_len), interval=25, blit=False, repeat=False)

	def Clear(self):
		print("OK")

#    def Plot(self):
#        x=0

	def init_window(self):

		self.master.title("SWE 1D Simulation")
		self.pack(fill='both', expand=1)     

#Create the controls, note use of grid

		self.labelXleft = Label(self,text="xleft",width=12)
		self.labelXleft.grid(row=0,column=1)
		self.labelXright = Label(self,text="xright",width=12)
		self.labelXright.grid(row=0,column=2)

		self.textXleft = Entry(self,width=12)
		self.textXleft.grid(row=1,column=1)
		self.textXright = Entry(self,width=12)
		self.textXright.grid(row=1,column=2)


#       self.buttonPlot = Button(self,text="Plot",command=self.Plot,width=12)
		self.buttonPlot = Button(self,text="Run",command=self.Run,width=12)
		self.buttonPlot.grid(row=2,column=1)
		self.buttonClear = Button(self,text="Clear",command=self.Clear,width=12)
		self.buttonClear.grid(row=2,column=2)

#       self.buttonClear.bind(lambda e:self.Plot)
		self.buttonClear.bind(lambda e:self.Clear)

		tk.Label(self,text="Shallow Water Equation 1D Simulation", font=("Arial Bold", 30)).grid(column=0, row=3)

		self.xleft   = 0.
		self.xright  = 1.
		self.Nx      = 100
		self.x = np.linspace(self.xleft, self.xright, self.Nx + 1)
		self.y = initial_Condition(self.x, 0.001, 0.5, 0.05)
		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(1,1,1)
		self.ax = plt.axes(xlim=(0, 1), ylim=(-0.002, 0.002))
		self.line, = self.ax.plot(self.x, self.y, lw=2)
		self.ax = plt.grid(True)
		self.ax = plt.title("Initial Condition")

		self.canvas = FigureCanvasTkAgg(self.fig, master=self)
		self.canvas.get_tk_widget().grid(column=0,row=4)

root = tk.Tk()
root.geometry("1000x600")
app = Window(root)
tk.mainloop()
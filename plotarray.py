from __future__ import print_function, division
from matplotlib.animation import FuncAnimation
from matplotlib.pyplot import subplots
from matplotlib import pyplot as pp

import numpy as np
from numpy import asarray,arange,zeros, int64
from h5py import File
from win32con import ARABIC_CHARSET



class plotanimation(object):
    
    '''
    plotanimation(X, Y, r=None, update=50)
    X = Either a single vector or a list of vectors, for first axis
    Y = Either a single vector or a list of vectors, for secong axis
    r = Optional, select indices in X to display. Default is all,
        but may be a range or a vector of indices
      
      
    Always run setanimation() before opening the plot  
    
    Example:
        win = plotanimation([z], [raw, Ha,Hb], update=100)
        
        win.set_format(0, line='-', marker='', color='blue')
        win.set_format(1, line='-', marker='', color='green')
        win.set_format(2, line='-', marker='', color='red')
        win.legend(['raw','Hs','Hf'])
        win.set_xlim(0, z[-1])
        win.set_ylim(-10, 1500)
        win.set_text('Echo',pos=(1.1,0.5),txtdata=None)
        win.grid()
        win.setanimation()
        pp.show()
        
    '''
    
    __x = None
    __xs = None
    __y = None
    __ys = None
    __ax = None
    __fig = None
    __line= None
    __i = 0
    __int = 0
    __range = None
    __text = None
    __textpos = (1.05, 0.9)
    __txtdata = None
    __view = None
    
    def __init__(self, x, y, r=None, update=50):
        ''' 
            X is a vector,  M x 1
            Y is a array,  M x N 
        '''
        
        self.__fig, self.__ax = subplots(nrows=1)
        self.__x = x
        if isinstance(x, list):
            self.__xs = [asarray(e).shape[0] for e in x]
        else:    
            
            self.__xs = self.__x.shape
        
        
        if isinstance(y, list):
            self.__y = [asarray(e) for e in y]
        else: 
            self.__y = y
        
        
        if isinstance(x, list):
            self.__ys = [asarray(e).shape[0] for e in y]
        else:    
            
            self.__ys = self.__y.shape
        
        self.__int = update
        self.__line  = []
        
        for j in range(len(self.__y)):
            l, = self.__ax.plot([],[])
            self.__line.append(l)
            
        if r == None:  
            r = arange(0, self.__ys[0], dtype=int64)
        elif asarray(r).shape[0] ==2:
            r = arange(r[0], r[1]).astype(int64)
        
        self.__range = r
        self.__view = zeros(4)
        self.__text = self.__ax.text(self.__view[0]*self.__textpos[0] , self.__view[3]*self.__textpos[1],' ')
        
    def set_format(self, i, line='-', marker='x',color='blue'):
        
        self.__line[i].set_linestyle(line)
        self.__line[i].set_marker(marker)
        self.__line[i].set_markerfacecolor(color)
        if marker == 'x':
            self.__line[i].set_markeredgecolor(color)
        else:
            self.__line[i].set_markeredgecolor('black')
            
    def set_xlabel(self,xlabel):
        self.__ax.set_xlabel(xlabel)

    def set_ylabel(self, ylabel):
        self.__ax.set_ylabel(ylabel)

    def set_text(self, title='Ping', pos = (1.05, 0.9), txtdata=None):
        '''  Set title and relative position  in view area'''
        self.__txtdata = txtdata
        self.__textpos = pos
        self.__textstring = title
        self.__text.set_position(self.__textpos)
        self.__text.set_position((self.__view[0]*self.__textpos[0],self.__view[3]*self.__textpos[1]))

    def set_xlim(self, start, stop):
        self.__view[0] = start
        self.__view[1] = stop
        
        self.__text.set_position((self.__view[0]*self.__textpos[0],self.__view[3]*self.__textpos[1]))
        
        self.__ax.set_xlim(start, stop)
        
    def set_ylim(self, start, stop):   
        self.__view[2] = start
        self.__view[3] = stop
        self.__text.set_position((self.__view[0]*self.__textpos[0],self.__view[3]*self.__textpos[1]))     
        
        self.__ax.set_ylim(start, stop)

    def legend(self,labels, loc = 'best'):
        self.__ax.legend(labels, loc=loc)
    
    def grid(self, switch=True):
        self.__ax.grid(switch)

    def init(self):
        if isinstance(self.__txtdata, type(None)):
            self.__text.set_text(self.__textstring + ': #' + str(self.__i))
        else:
            self.__text.set_text(self.__textstring + ': #' + str(self.__txtdata[self.__i]))
        
        for j in range(len(self.__line)):
            yndim = asarray(self.__y[j]).ndim
            if len(self.__xs) > 1:
                
                xndim = asarray(self.__x[j]).ndim
                            
                if xndim < yndim:
                    self.__line[j].set_xdata(self.__x[j])
                elif xndim == yndim:
                    self.__line[j].set_xdata(self.__x[j][self.__i])
            else:
                self.__line[j].set_xdata(self.__x)
            self.__line[j].set_ydata(self.__y[j][self.__i])
        
        return self.__line + [self.__text]
    
    def animate(self,i):
        self.__i = i
        
        if isinstance(self.__txtdata, type(None)):
            self.__text.set_text(self.__textstring + ': #' + str(self.__i))
        else:
            self.__text.set_text(self.__textstring + ': #' + str(self.__txtdata[self.__i]))
        for j in range(len(self.__line)):
            
            yndim = asarray(self.__y[j]).ndim
            if len(self.__xs) > 1:
                xndim = asarray(self.__x[j]).ndim
                if xndim < yndim:
                    self.__line[j].set_xdata(self.__x[j])
                elif xndim == yndim:
                    self.__line[j].set_xdata(self.__x[j][self.__i])
            #else:
            #    self.__line[j].set_xdata(self.__x)
            
            self.__line[j].set_ydata(self.__y[j][self.__i])
        
        return self.__line + [self.__text]
    
    def getaxis(self):
        return self.__ax
    
    def setanimation(self):
        ''' Run this function before opening any window! '''
        self.__ani = FuncAnimation(self.__fig, self.animate, self.__range, 
                                             init_func=self.init, interval=self.__int, blit=True)




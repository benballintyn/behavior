#"skorne example of reading the MCP3008 analog input channels and printing
# them all out.
# Author: Tony DiCola
# License: Public Domain
import time
import csv
import numpy as np
import easygui
import datetime
from matplotlib import pyplot as plt

# Import SPI library (for hardware SPI) and MCP3008 library.
import Adafruit_GPIO.SPI as SPI
import Adafruit_MCP3008

# Ask the user for today's date
date = easygui.multenterbox(msg = 'Enter todays date (all integers)', fields = ['yy', 'mm', 'dd'])
for i in range(len(date)):
	date[i] = int(date[i])

# Get the mouse's name
mice_names = easygui.multenterbox(msg = 'Enter the mice names', fields = [str(i+1) for i in range(8)])

# Make the output data files
g = [open('./Desktop/data/%02i%02i%02i%s_licks.txt' % (date[0], date[1], date[2], mice_names[i]), 'w') for i in range(8)]

# Make a list of active mice
active_chans = [0 for i in range(8)]

# Duration of experiment
duration = 65
baselineDuration = 5
# Hardware SPI configuration:
SPI_PORT   = 0
SPI_DEVICE = 0
mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))

print('Reading lick spout values, press Ctrl-C to quit...\nThis experiment will run for %s minutes.' % duration)
expt_start_times = np.ones((1,len(active_chans)))
expt_start_times *= time.time()
first_start = np.min(expt_start_times)
Vema1 = np.zeros((1,len(active_chans)))
Vema2 = np.zeros((1,len(activ_chans)))
alpha1 = .1
alpha2 = .001
baseMeans = np.zeros((1,len(active_chans)))
baseVars = np.zeros((1,len(active_chans)))
baseNs = np.zeros((1,len(active_chans)))
zscores = np.zeros((1,len(active_chans)))
inLick = np.zeros((1,len(active_chans)))
Nlicks = np.zeros((1,len(active_chans)))
plotTcounter = 0
fig, ax = plt.subplots()
rects = ax.bar(range(len(active_chans)),np.zeros((1,len(active_chans))))
while (time.time() - first_start < duration*60):
    plotTcounter+=1
    # Read all the ADC channel values in a list.
    values = [0]*8
    for i in range(8):
        # The read_adc function will get the value of the specified channel (0-7).
        values[i] = mcp.read_adc(i)
        print>>g[i], values[i], '\t', "%f" % time.time()

    for i in range(8):
        Vema1[i] = alpha1*values[i] + (1-alpha1)*Vema1[i]
        Vema2[i] = alpha2*values[i] + (1-alpha2)*Vema2[i]

    if (time.time() - first_start < baselineDuration*60):
        for i in range(8):
            baseNs[i]+=1
            oldMean = baseMeans[i]
            baseMeans[i]+=(Vema1[i] - oldMean)/baseNs[i]
            baseVars[i]+=(Vema1[i] - baseMeans[i])*(Vema1[i] - oldMean)

    if (time.time() - first_start > baselineDuration*60):
        for i in range(8):
            zscores[i] = (Vema1[i] - baseMeans[i])/np.sqrt(baseVars[i]/baseNs[i])
            if (!inLick[i] & (zscores[i] > 4) & Vema1[i] > Vema2[i]):
                inLick[i] = 1
                Nlicks[i] += 1
            else if (inLick[i] and ((zscores[i] < 4) or (Vema1[i] < Vema2[i]))):
                inLick[i] = 0
    
    if (plotTcounter == 2000):
        plotTcounter = 0
        [rect.set_height(h) for rect,h in zip(rects,Nlicks)]
        fig.canvas.draw()

    #lickdata = np.vstack((lickdata, values)
    
    # Pause for 1/2 a millisecond second.
    time.sleep(0.0005)
    
for i in range(len(g)):
    g[i].close() 


#"skorne example of reading the MCP3008 analog input channels and printing
# them all out.
# Author: Tony DiCola
# License: Public Domain
import time
import csv
import numpy as np
import easygui
import datetime

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
active_mice = [0 for i in range(8)]

# Mouse start time counter list
start_counter = [0 for i in range(8)]

# Duration of experiment
duration = 30

# Hardware SPI configuration:
SPI_PORT   = 0
SPI_DEVICE = 0
mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))

print('Reading lick spout values, press Ctrl-C to quit...\nThis experiment will run for %s minutes.' % duration)

while True:

    # open up the activity file and count the number of lines in it
    f = open('activity.txt', 'r')
    lines = 0
    for line in f.readlines():
	lines = lines + 1
    f.close()
    print lines
    for i in range(len(active_mice)):
	if active_mice[i] == 0 and lines >= i+1:
		active_mice[i] = 1
		start_counter[i] = time.time()

    # Read all the ADC channel values in a list.
    values = [0]*8
    for i in range(8):
        # The read_adc function will get the value of the specified channel (0-7).
        values[i] = mcp.read_adc(i)
        if active_mice[i] > 0:
        	print>>g[i], values[i], '\t', "%f" % time.time()
        	if time.time() - start_counter[i] >= 60*duration:
			active_mice[i] = -1
    #lickdata = np.vstack((lickdata, values))
    
    # Pause for half a second.
    time.sleep(0.001)

    print active_mice
    
    if np.all(active_mice == [-1 for i in range(8)]):
    	break

for i in range(len(g)):
    g[i].close() 





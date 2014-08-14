# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 10:06:40 2013

@author: ttsai

To understand Plume Scanner data
Takes in either pickled or unpickled data files
2. Takes the closest 0 to determine the final tape beginnings and endings

Input file columns:
    0: Time
    1: CH4_dry concentration
    2: H2O concentration
    3: Valve Mask
    4: "SpectrumID"
    5:["WS_WIND_LAT"]
    6:"WS_WIND_LON"]
    7:["CAR_SPEED"])
    8: GPS Lat
    9: GPS Long
    10: CO2

Output:
    0: Filename
    1: Carspeed
    2: Std Carspeed
    3: Weighted average lateral wind speed
    4: Std lateral wind speed
    5: Peaks realignment
    6: Startends
    7: Background
    8: Peak values
    9: Hellman factor
    10: ans
    11: sec ans
    12: plumestartends
    13: lat
    14: long
    15: FWHM
    16: Overhead flag
    17: Line integral ratio value
    18: Overhead concentration ratio value
    19: Wind dir
    20: Std Wind dir

    
v22 Has weighted average
v24 Takes into account the CO2 push gas and the CFADS upgrade 
Output changed as windspeed is swapped for weightedaverage
Average carspeed, stdcar, and stdwind ONLY includes those interpolated
Replace error with plume startends
v25 Includes mean wind dir and std of wind dir. Wind dir is determined by using WIND_N and WIND_E 
parameters taken from pickled file
Outputs also now include line integral ratio and concentration ratio to determine overheadness
"""

from numpy import *
import easygui as eg
import string
import pdb
import pylab
import scipy
import scipy.integrate
import time
from SavitzkyGolay import savitzky_golay as sv
from scipy.interpolate import griddata

import little_reader
#filenamepickle = r'C:\Users\ttsai\Desktop\Australia\data\plume_1378349356_2013_9_4__19_49_unpickled.txt'
#filenamepickle = eg.fileopenbox(msg="File to unpickle")
filedir = r'C:\Users\ttsai\Desktop\PlumeScanner6\batch'
usefulfilename = \
'plume_1381968767_2013_10_16__17_12_unpickled.txt'


filenamepickle = filedir+'\\'+usefulfilename

if filenamepickle[len(filenamepickle)-3:] == 'txt':
    filenamepickle = filenamepickle.replace("_unpickled.txt", ".dat")
    
filename =  little_reader.readfunc(filenamepickle ,'n', 'y', "CH4_dry")[0]

#filename = r'C:\Users\ttsai\Desktop\Australia\data\CFADS2206-20130910-010951Z-DataLog_User_Minimal_minimal.txt'
upside = 'n'

backgroundyn ='n'
backgroundset = \
1.908

plumeset = 'n'
plumestart =27
plumeend = 75
windstartset = 'n'
winds = \
[41,71]   

fwhmcalc = 'y'

print filename

a = 0.17 #Hellman factor for vertical wind calculations
def windfunc (a):
    """
    Power Law equation which returns scaling factor
    Anemoeter on Toyota 4Runner is measured 111 in above the ground
    """
    #pixelheights = [1.01, 2.23, 3.45, 4.67] #For Australia
    #pixelheights = [0, 1.01, 2.23, 3.45, 4.67] #For Australia
    pixelheights = [0.401, 0.401+0.635,0.401+(0.635*2), 0.401+(0.635*3), 0.401+(0.635*4), 0.401+(0.635*5)] #For 6 pixels
    #pixelheights = [0.5, 1.72, 2.94, 4.16]
    #pixelheights = [0, 0.5, 1.72, 2.94, 4.16] #To include ground in calculations
    #pixelheights = [0, 0.432, 1.09, 1.75, 2.44] #Plume Scanner Jr
    #ha = 2.82
    ha = 2.98
    w = [(i/ha)**a for i in pixelheights]
    return w, pixelheights

w = windfunc(a)[0]
pixelheights = windfunc(a)[1]   
#print "Adjust wind factors are: ", w

data = loadtxt(filename)
valves = [x for x in data[:,3]] #converts the ndarray into a list
onindex = valves.index(5) #determines the index when the valves turn on completely
bufferonindex = onindex+1 #Give ten points extra due to trigger happening so close to peak
pretape = data[:bufferonindex,:] 
print "Data truncated at row %i" %onindex
tapes = data[onindex:,:] 
samplerate = average(data[1:,0]-data[:len(data)-1,0])
#pylab.plot(tapes[:,1])
#pylab.show()

def searchtapesbackground (tape, thres, upside = 'n'):
    """
    Utilizing background to find startends. This will only work if zero gas is negative
    """
    if upside == 'y':
        d = (tape-thres)*-1
    else:
        d = tape-thres
    #inisearch = [200, 400, 600]
    #inisearch = [170, 325, 489, 643, 806]  
    inisearch = [150, 300, 450, 600, 750] 
    refinedsearch = [0]
    for i in inisearch:
        refinedsearch.append(argmin(d[i:i+100])+i)
    refinedsearch.sort()
    refinedsearch.append(900)
    #print "Refined search areas are: ", refinedsearch
    secderiv = gradient(gradient(sv(tape, 15,3)))
    startend = []        
    for i in refinedsearch:
        if i >50:
            ini = i-50
        else: 
            ini = i
        fin = i+50       
        temp = secderiv[ini:fin].argsort()[:2]+ini
        #to ensure that adjacent points aren't chosen
        count = 1
        while abs(temp[0]-temp[1]) <= 10: 
            count +=1
            temp[1] = secderiv[ini:fin].argsort()[count+1]+ini
        temp[0] = temp[0]
        temp[1] = temp[1]
        startend.extend(temp) 
#        pylab.plot(tape[ini:fin])
#        pylab.twinx()
#        pylab.plot(secderiv[ini:fin])
#        pylab.show()
#        pdb.set_trace()
    return startend
    


#startend = searchtapesbackground(tapes[:,1],min(pretape[:,1]))
#pdb.set_trace()
startend = searchtapesbackground(tapes[:,10],min(pretape[:,2]), upside)
startend = startend[1:len(startend)-1]
startend.sort()
print startend
if startend[0] <10:
    startend[0] = startend[0]+20
    print startend

def searchbetter (tape, tempstartend, threshold, buff, upside = 'n'):
    """
    Find nearest "0"
    Assuming that it is negative zero gas
    """
    #pdb.set_trace()
    startend =[]
    secderiv = gradient(gradient(tape)) 
    #threshold for determining 0ness
    for index, item in enumerate(tempstartend):
        if index%2 != 0: #odd number indices
            temptest = secderiv[item-buff:item]
            if upside == 'n':
                startend.append(item-buff+[ind for ind, it in enumerate(temptest) if abs(it) <threshold][-1]) #assuming that it is monotonically decreasing        
            else:
                startend.append(item-buff+[ind for ind, it in enumerate(temptest) if abs(it) <threshold][0]) #assuming that it is monotonically increasing        
        else:
            temptest = secderiv[item:item+buff]
            if upside == 'n':            
                startend.append(item+[ind for ind, it in enumerate(temptest) if abs(it) <threshold][0]) #assuming that it is monotonically increasing
            else:
                startend.append(item+[ind for ind, it in enumerate(temptest) if abs(it) <threshold][-1]) #assuming that it is monotonically decreasing                    
    return startend

try:    
    #startendnew = searchbetter(tapes[:,1], startend, 0.01, 10) #for nice peaks
    startendnew = searchbetter(tapes[:,10], startend, 0.05, 50, upside) #Using CO2
    #startendnew = [18.0, 155.0, 169.0, 306.0, 321.0, 457.0, 473.0, 609.0, 626.0, 759.0, 776.0, 910.0]
    #startendnew = searchbetter(tapes[:,1], startend,0.05, 15, upside) 
except:
    print "Automatic startends fail"    
    startendnew = startend

import Tkinter
root = Tkinter.Tk()
root.withdraw()
import tkSimpleDialog

#pdb.set_trace()
pylab.figure("Tape startpoints")
pylab.subplot(3,1,1)#Top plot will be with CH4
pylab.ylabel("CH4 search")
pylab.xlabel("Row numbers")
pylab.scatter(startendnew, tapes[startendnew,1])
pylab.plot(tapes[:,1])
pylab.subplot(3,1,2)#CO2
pylab.ylabel("C02 search")
pylab.scatter(startendnew, tapes[startendnew,10])
pylab.plot(tapes[:,10])
pylab.ion()
pylab.show()

tempstartend = tkSimpleDialog.askstring("Fix startends", str(startend), initialvalue = str(startend))
startendnew = [float(i.strip()) for i in tempstartend.lstrip('[').rstrip(']').split(',')]
startend = startendnew

pylab.subplot(3,1,3)
#pylab.plot(gradient(gradient(sv(tapes[:,1], 15,3))))
#pylab.scatter(startend, gradient(gradient(sv(tapes[:,1], 15,3)))[startend])
pylab.ylabel("CH4 another search")
pylab.xlabel("Row numbers")
pylab.scatter(startendnew, tapes[startendnew,1])
pylab.plot(tapes[:,1])
pylab.ioff()


tapeA = tapes[startend[1]:startend[0]:-1, :]    # TapeA is played in reverse so reversed here
tapeB = tapes[startend[2]:startend[3], :]
tapeC = tapes[startend[5]:startend[4]:-1, :]    # TapeC is played in reverse so reversed here
tapeD = tapes[startend[6]:startend[7], :]
tapeE = tapes[startend[9]:startend[8]:-1, :]    # TapeE is played in reverse so reversed here
tapeF = tapes[startend[10]:startend[11], :]
#ratio = (390.0/1070.0) #Ratio of recording flow over playback flow
ratio = 1050.0/1000.0

"""
Wind selection for average; use first derivatives in pretape to determine
where the CH4 plume occurred and use those wind measurements
"""
firstderiv = diff(pretape[:,1])
if plumeset == 'n':
    tempplumestart = argmax(firstderiv)
    tempplumeend =argmin(firstderiv[5:])+5
    plumestartend = sort(searchbetter(pretape[:,1], [tempplumestart, tempplumeend], 1, 7,'y'))
    plumestart = plumestartend[0]
    plumeend = plumestartend[1]

#delay = 4*2  # approx delay is 4 s, considering about 2 pts per second
delay = 0 
if windstartset == 'n':
    windstart = plumestart + delay
    windend = plumestart + delay + (plumeend-plumestart)
else:
    windstart = winds[0]
    plumestart = windstart-delay
    windend = winds[1]
    plumeend = windend-delay
#windstart = 21
#windend=10
#windend = len(pretape)
print "max: %i, min: %i" %(windstart, windend)
#pylab.figure("Carspeed and windspeed")
#ax1 = pylab.subplot(2,1,1)#Top plot will be with CH4
#pylab.ylabel("CH4 pretape")
#pylab.xlabel("Row numbers")
#pylab.plot(data[:100,1])
#pylab.scatter([plumestart, plumeend], [data[plumestart,1], data[plumeend,1]])
#pylab.twinx()
#pylab.plot(firstderiv, 'r')
#pylab.subplot(2,1,2, sharex = ax1)
#pylab.ylabel("Wind")
#pylab.xlabel("Row numbers")
#pylab.plot(data[:100,5])
#pylab.scatter([windstart, windend], [data[windstart,5], data[windend,5]])
#pylab.xlim([0, len(data[:100])])
#pylab.show()
stdwind = std(data[windstart:windend,5])
carspeed = average(data[plumestart:plumeend,7])
stdcar = std(data[plumestart:plumeend,7])

def scrunchtime (tape, ratio):
    difftime = (abs(tape[1:,0]-tape[:len(tape)-1,0]))*ratio
    timenew = [0]
    for index, x in enumerate(difftime):
        timenew.append(timenew[index]+x)
    for index, item in enumerate(tape[:,0]):
        tape[index,0] =timenew[index]
    return tape
    
#for plotting
tapeAplot = scrunchtime(tapeA,ratio)
tapeBplot = scrunchtime(tapeB,ratio)
tapeCplot = scrunchtime(tapeC,ratio)
tapeDplot = scrunchtime(tapeD,ratio)
tapeEplot = scrunchtime(tapeE,ratio)
tapeFplot = scrunchtime(tapeF,ratio)
#ynsave = raw_input("Save reconstructed tapes? (y/n) \n")
#ynsave = 'n'
#if ynsave == "y":    
##Save without background subtracted with respect to time
#    savetxt(filename[:len(filename)-4]+"auto_tapeA.txt", transpose(vstack([scrunchtime(tapeA,ratio)*carspeed, tapeA[:,1]])) , fmt='%10f', delimiter="\t")
#    savetxt(filename[:len(filename)-4]+"auto_tapeB.txt", transpose(vstack([scrunchtime(tapeB,ratio)*carspeed, tapeB[:,1]])), fmt='%10f', delimiter="\t")
#    savetxt(filename[:len(filename)-4]+"auto_tapeC.txt", transpose(vstack([scrunchtime(tapeC,ratio)*carspeed, tapeC[:,1]])) , fmt='%10f', delimiter="\t")
#    savetxt(filename[:len(filename)-4]+"auto_tapeD.txt", transpose(vstack([scrunchtime(tapeD,ratio)*carspeed, tapeD[:,1]])) , fmt='%10f', delimiter="\t")
#    savetxt(filename[:len(filename)-4]+"auto_tapeE.txt", transpose(vstack([scrunchtime(tapeE,ratio)*carspeed, tapeE[:,1]])) , fmt='%10f', delimiter="\t")
#    savetxt(filename[:len(filename)-4]+"auto_tapeF.txt", transpose(vstack([scrunchtime(tapeD,ratio)*carspeed, tapeF[:,1]])) , fmt='%10f', delimiter="\t")

def backgroundcalc (pretape, tapearray):
    """
    2% oflowest values. It only takes in 1D array of ch4
    """
    pretapemin = min(pretape) #lowest 2% of 50 pts is barely 1 pt
    backgroundarray = [pretapemin]
    for i in tapearray:
        if (min(i) < pretapemin*(1-0.05)) or (min(i) > pretapemin*(1+0.05)): #if the background in tapes is 5% too different then its prob zero gas
            print "Poor tape cut affecting background calc"
        else:
            #cutoff = int(round(0.02*len(i)))
            cutoff = int(round(0.5*len(i))) #Had to use 5% to account for the poor tape cutting
            backgroundarray.append(average(sort(i)[5:cutoff])) #5pt buffer for poor tape cutting
    background = min(backgroundarray)
    return background
    
"""
!!!! BACKGROUND HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
if backgroundyn == 'y':
    print "Using set background"    
    background = backgroundset
else:
    print "Using automatic background"
    background = round(backgroundcalc(pretape[:50,1], [tapeAplot[:,1], tapeBplot[:,1], tapeCplot[:,1], tapeDplot[:,1]])*1000)*10**-3 #(Rounding to the third decimal)
print "Background is %f ppm" %background

"""
Peak aligning
"""
pretaperecon =scrunchtime(pretape,1)
peaklocation = pretaperecon[argmax(pretaperecon[:,1]),0]
#temppretaperecon = pretaperecon[:len(pretaperecon)-10,:]
#peaklocation = temppretaperecon[argmax(temppretaperecon[:,1]),0]
print "peak location"
print peaklocation
#pdb.set_trace()
tapedeck = [tapeAplot, tapeBplot, tapeCplot, tapeDplot, tapeEplot, tapeFplot]
peaks = []
for i in tapedeck:
    temppeak = i[argmax(i[:,1]),0]
    peaks.append(temppeak)
    
pylab.figure("peak realignment")
pylab.subplot(2,1,1)
pylab.xlabel("Distance")
pylab.ylabel("CH4 conc (ppm)")
pylab.plot(tapeAplot[:,0], tapeAplot[:,1], label = 'tapeA')
pylab.plot(tapeBplot[:,0], tapeBplot[:,1], label = 'tapeB')
pylab.plot(tapeCplot[:,0], tapeCplot[:,1], label = 'tapeC')
pylab.plot(tapeDplot[:,0], tapeDplot[:,1], label ='tapeD')
pylab.plot(tapeEplot[:,0], tapeEplot[:,1], label ='tapeE')
pylab.plot(tapeFplot[:,0], tapeFplot[:,1], label ='tapeF')
pylab.plot(pretaperecon[:,0], pretaperecon[:,1], label = 'pretape')
#pylab.legend(loc = 'upper left')
pylab.ion()
pylab.show()

temppeaks = tkSimpleDialog.askstring("Fix Peaks", str(peaks), initialvalue = str(peaks))
peaks = [float(i.strip()) for i in temppeaks.lstrip('[').rstrip(']').split(',')]

def peakrealigner(tape, peaklocation, peak):
    """
    takes in scrunchtimed tape
    """
    for index, item in enumerate(tape[:,0]):
        tape[index,0] = item +(peaklocation - peak)
    return tape
tapeAalign = peakrealigner(tapeAplot, peaklocation, peaks[0])    
tapeBalign = peakrealigner(tapeBplot, peaklocation, peaks[1])
tapeCalign = peakrealigner(tapeCplot, peaklocation, peaks[2])
tapeDalign = peakrealigner(tapeDplot, peaklocation, peaks[3])
tapeEalign = peakrealigner(tapeEplot, peaklocation, peaks[4])
tapeFalign = peakrealigner(tapeFplot, peaklocation, peaks[5])
tapedeck = [tapeAalign, tapeBalign, tapeCalign, tapeDalign, tapeEalign, tapeFalign]


def interpolator (tapedeck, pretaperecon):
    newx = pretaperecon[:,0]
    newtapedeck = []
    minarray = []
    maxarray = []
    for i in tapedeck:
        points = i[:,0]
        values = i[:,1]
        intertape = griddata(points, values, newx, method = 'linear')
        temptape = transpose(vstack([newx, intertape]))
        temptape =temptape[~isnan(temptape).any(1)] #get rid of NAN
        #print temptape[0,0], temptape[len(temptape)-1,0]        
        minarray.append(temptape[0,0])
        maxarray.append(temptape[len(temptape)-1,0])
        newtapedeck.append(temptape)
    #pdb.set_trace()
    print "Min time and max time on pretape: %f s and %f s" %(pretape[0,0], pretape[len(pretape)-1,0])
    mintime = max(minarray)
    maxtime = min(maxarray)
    newerx = []
    indexarray = []
    for index, item in enumerate(pretape[:,0]):
        if item <= maxtime and item >= mintime:
            newerx.append(item)
            indexarray.append(index)
    print "new min and max time: %f s and %f s" %(mintime, maxtime)
    newertapedeck = []
    for i in tapedeck:
        points = i[:,0]
        values = i[:,1]
        intertape = griddata(points, values, newx, method = 'linear')
        temptape = transpose(vstack([newerx, intertape[indexarray]]))
        newertapedeck.append(temptape)
    return newertapedeck, indexarray,mintime, maxtime

tempans = interpolator([tapeAalign, tapeBalign, tapeCalign, tapeDalign, tapeEalign, tapeFalign], pretaperecon)
newtapedeck = tempans[0]
indexarray = tempans[1]
mintime = tempans[2]
maxtime = tempans[3]
pretapetrunc = pretaperecon[indexarray,:]

def distanceGPS (pretape):
    """
    calculate distance based on GPS coords
    Assuming a flat surface since elevation isn't taken into account
    Col 8 is Lat
    col 9 is long
    """
    import math
    lat = pretape[:,8]*(math.pi/180)
    lon  = pretape[:,9]*(math.pi/180)
    dist = [0]
    for index in range(len(lat)-1):
        rd = 6371*10**3 #radius of earth
        #Using the spherical law of cosines
        tempdist = math.acos(math.sin(lat[index])*math.sin(lat[index+1]) + \
        math.cos(lat[index])*math.cos(lat[index+1]) *math.cos(lon[index+1]-lon[index])) * rd
        dist.append(tempdist+dist[index])
    return dist

#distancebeee
dist = distanceGPS(pretapetrunc)
#dist = [i*1.00 for i in range(len(newtapedeck[0][:,1]))] #fake movement given that car isn't moving
tapeAdist = transpose(vstack([dist, newtapedeck[0][:,1]]))
tapeBdist = transpose(vstack([dist, newtapedeck[1][:,1]]))
tapeCdist = transpose(vstack([dist, newtapedeck[2][:,1]]))
tapeDdist = transpose(vstack([dist, newtapedeck[3][:,1]]))
tapeEdist = transpose(vstack([dist, newtapedeck[4][:,1]]))
tapeFdist = transpose(vstack([dist, newtapedeck[5][:,1]]))

def backgrounder (tape, background):
    """
    subtract background without going negative
    """
    y = []
    x = [t for t in tape[:,0]]
    tracker = 0
    for item in tape[:,1]:
        if item - background <0:
            y.append(0)
            tracker +=1 
        else:
            y.append(item-background)
    print "Percentage negative: %f" %((tracker/float(len(tape)))*100)
    return transpose(vstack([x, array(y)]))
    
tapeArecon = backgrounder(tapeAdist, background)
tapeBrecon = backgrounder(tapeBdist, background)
tapeCrecon = backgrounder(tapeCdist, background)
tapeDrecon = backgrounder(tapeDdist, background)
tapeErecon = backgrounder(tapeEdist, background)
tapeFrecon = backgrounder(tapeFdist, background)

pylab.subplot(2,1,2)
pylab.xlabel("Reconstructed tapes without background (m)")
pylab.ylabel("CH4 conc (ppm) - %f ppm" %background)
pylab.plot(tapeArecon[:,0], tapeArecon[:,1], label = 'tapeA')
pylab.plot(tapeBrecon[:,0], tapeBrecon[:,1], label = 'tapeB')
pylab.plot(tapeCrecon[:,0], tapeCrecon[:,1], label = 'tapeC')
pylab.plot(tapeDrecon[:,0], tapeDrecon[:,1], label ='tapeD')
pylab.plot(tapeErecon[:,0], tapeErecon[:,1], label ='tapeE')
pylab.plot(tapeFrecon[:,0], tapeFrecon[:,1], label ='tapeF')
#pylab.legend(loc = 'upper right')
pylab.ioff()

def truncarray (taperecon, start, end):
    truncated = []
    truncateindex = []
    for i in range(len(taperecon)):
        if (taperecon[i,0] < end and taperecon[i,0] > start):
            truncated.append(taperecon[i])
            truncateindex.append(i)
    return array(truncated), truncateindex
 
#truncate = raw_input("Trucate? \n")
truncate = 'n'
if truncate == 'y':
    startdist = float(raw_input("Start distance \n"))
    enddist = float(raw_input("End distance \n"))
    pylab.figure("Truncated tapes")
    pylab.xlabel("Distance (m)")
    pylab.ylabel("CH4 conc (ppm)")
    indexarray = truncarray(tapeArecon, startdist, enddist)[1]
    tapeArecon = truncarray(tapeArecon, startdist, enddist)[0]
    tapeBrecon = truncarray(tapeBrecon, startdist, enddist)[0]
    tapeCrecon = truncarray(tapeCrecon, startdist, enddist)[0]
    tapeDrecon = truncarray(tapeDrecon, startdist, enddist)[0]
    pylab.plot(tapeArecon[:,0], tapeArecon[:,1], label = 'tapeA')
    pylab.plot(tapeBrecon[:,0], tapeBrecon[:,1], label = 'tapeB')
    pylab.plot(tapeCrecon[:,0], tapeCrecon[:,1], label = 'tapeC')
    pylab.plot(tapeDrecon[:,0], tapeDrecon[:,1], label ='tapeD')
    pylab.plot(tapeErecon[:,0], tapeErecon[:,1], label ='tapeE')
    pylab.legend(loc = 'upper left')
    pylab.show()
###############################
"""
Wind here
"""
windindex = [x+delay for x in indexarray]
windlat = data[windindex,5]
pylab.figure("Carspeed and windspeed")
ax1 = pylab.subplot(3,1,1)#Top plot will be with CH4
pylab.ylabel("CH4")
pylab.plot(data[:bufferonindex,1])
pylab.scatter([plumestart, plumeend], [data[plumestart,1], data[plumeend,1]],s=40,marker = '^',c='r')
pylab.subplot(3,1,2, sharex = ax1)
pylab.ylabel("Carspeed")
pylab.plot(data[:bufferonindex,7])
pylab.scatter(indexarray, data[indexarray,7])
pylab.scatter([plumestart, plumeend], [data[plumestart,7], data[plumeend,7]],s=40,marker = '^',c='r')
pylab.subplot(3,1,3, sharex = ax1)
pylab.ylabel("Wind")
pylab.xlabel("Row numbers")
pylab.plot(data[:bufferonindex,5])
pylab.scatter(windindex, data[windindex,5])
pylab.scatter([windstart, windend], [data[windstart,5], data[windend,5]],s=40, marker = '^', c='r')
pylab.xlim([0, len(data[:bufferonindex])])
pylab.ion()
pylab.show()
pylab.ioff()

##################################################
"""
Weighted average
"""
def weightaverage(tapeB, windlat):   
    """
    Use tapeB since co-located with instrument monitoring inlet tube so I can 
    avoid using pretape
    For simplicity, ignore the power law
    """
    #winner = max(tapeB[:,1])    
    #numerator = sum([(x/winner) for x in tapeB[:,1]]*windlat)
    numerator = sum([(x) for x in tapeB[:,1]]*windlat)
    denominator = sum([x for x in tapeB[:,1]])
    return average(numerator/denominator)
#weightaverage = weightaverage([tapeArecon, tapeBrecon, tapeCrecon, tapeDrecon], windlat, w)
#weightaverage = weightaverage(tapeBrecon, windlat)
weightaverage = weightaverage(tapeErecon, windlat)
print "\nWeighted average is %f" %weightaverage

#Use Simpson rule for integration. The format is integrate y over x axis is simps.(y,x)
#Utilize a wind correction factor, w
def piecewiseIntegrator (tape, windfactor, windlat):
    sums = []
    for i in range(len(tape)-2):
        windbit = average(windlat[i:i+2])
        sums.append(scipy.integrate.simps(tape[i:i+2,1], tape[i:i+2,0])*windfactor*windbit)
    return sum(sums)

"""
Took mean and rms calculation from David.
Used pretape pickled data on WIND N and WIND E data
"""

d2r = math.pi/180.
def yamartinoMeanAndRMS(values):
    '''Calculate the Yamartino mean and standard deviation. Input values (and output) 
    are in degrees. '''

    if values is None or len(values)<1:
        print "WARNING: TimeSeriesCalculator: yamartinoMeanAndRMS: len(values) < 1"
        return 0., 0.
    d2r = math.pi/180.
    
    meanSin = mean( [sin(d2r*v) for v in values] )
    meanCos = mean( [cos(d2r*v) for v in values] )
    meanAngle = arctan2( meanSin, meanCos )
    epsilon = math.sqrt( 1. - (meanSin*meanSin + meanCos*meanCos) )
    rms = arcsin(epsilon) * (1. + (2./math.sqrt(3.) -1.)*epsilon*epsilon*epsilon)
 
    return meanAngle/d2r, rms/d2r
windN=[]
windE=[]
angledeg  = []
import cPickle
try:
    with open(filenamepickle, 'rb') as fileo:
        while True:
            try:
                obj = cPickle.load(fileo)
                if (obj["data"]["SpectrumID"] == 25) and (obj["data"]["ValveMask"] == 0):
                    windN.append(obj["data"]["WIND_N"])
                    windE.append(obj["data"]["WIND_E"])
                    angledeg.append(arctan2(obj["data"]["WIND_N"], obj["data"]["WIND_E"])/d2r)
            except:
                break
except:
    True
tempans = yamartinoMeanAndRMS(angledeg)
meanwinddir = tempans[0]
stdwinddir = tempans[1]             
        
#tapeAint = piecewiseIntegrator(tapeArecon, w[1], windlat)
#tapeBint = piecewiseIntegrator(tapeBrecon, w[2], windlat)
#tapeCint = piecewiseIntegrator(tapeCrecon, w[3], windlat)
#tapeDint = piecewiseIntegrator(tapeDrecon, w[4], windlat)
#tapeEint = piecewiseIntegrator(tapeErecon, w[5], windlat)
#tapeFint = piecewiseIntegrator(tapeFrecon, w[6], windlat)
#zarray = transpose(vstack([pixelheights, array([0, tapeAint, tapeBint, tapeCint, tapeDint, tapeEint, tapeFint])]))
tapeAint = piecewiseIntegrator(tapeArecon, w[0], windlat)
tapeBint = piecewiseIntegrator(tapeBrecon, w[1], windlat)
tapeCint = piecewiseIntegrator(tapeCrecon, w[2], windlat)
tapeDint = piecewiseIntegrator(tapeDrecon, w[3], windlat)
tapeEint = piecewiseIntegrator(tapeErecon, w[4], windlat)
tapeFint = piecewiseIntegrator(tapeFrecon, w[5], windlat)
zarray = transpose(vstack([pixelheights, array([tapeAint, tapeBint, tapeCint, tapeDint, tapeEint, tapeFint])]))
#pylab.figure("Y line integrations wrt height")
#pylab.xlabel("Height (m)")
#pylab.ylabel("Y line integration")
#pylab.plot(zarray[:,0], zarray[:,1])
#pylab.show()
#Use trapezoidal integration
ans = abs(scipy.integrate.trapz(zarray[:,1], x = zarray[:,0])*10**-6)
print "%f liters per sec" %(ans*10**3)
#err = sqrt((stdcar/carspeed)**2+(stdwind/windspeed)**2) #Error is calculated only for wind and car and is decimal percentage
err = -1

############Calculate line integrals to determine if plume is overhead
fakeones = array([1 for i in range(len(tapeArecon))])
Aint = piecewiseIntegrator(tapeArecon, 1.0, fakeones)
Bint = piecewiseIntegrator(tapeBrecon, 1.0, fakeones)
Cint = piecewiseIntegrator(tapeCrecon, 1.0, fakeones)
Dint = piecewiseIntegrator(tapeDrecon, 1.0, fakeones)
Eint = piecewiseIntegrator(tapeErecon, 1.0, fakeones)
Fint = piecewiseIntegrator(tapeFrecon, 1.0, fakeones)
lineoverheadratio = Fint/sum([Aint,Bint,Cint,Dint,Eint,Fint])
print "Line integral overhead ratio is (>0.13 is overheady): %f" %lineoverheadratio
if lineoverheadratio > 0.17:
    lineflag = -1
else:
    lineflag = 0
peakseval = [max(tapeArecon[:, 1]), \
    max(tapeBrecon[:, 1]),\
    max(tapeCrecon[:, 1]),\
    max(tapeDrecon[:, 1]), \
    max(tapeErecon[:, 1]), \
    max(tapeFrecon[:, 1])]
winningpixel = argmax(peakseval)
winningpeak = peakseval[winningpixel]
overheadratioc =  winningpeak/peakseval[5]
print "Overhead ratio in concentration(overhead conc ratio <1.7 shady): %f " %overheadratioc
if overheadratioc < 1.7:
    concflag = -1
else:
    concflag = 0
overheadflag = lineflag + concflag
print "Overheadflag: %s" %(overheadflag)
if winningpixel == 5:
    print 'Plume overhead: Top pixel has highest concentration'
elif overheadratioc <= 1.7:
    print 'Plume likely overhead, ratio: %f' %(overheadratioc)
    
"""
According the weather.gov, at 11:53am in Vernal, UT
Temperature = 22F or -5.6C
Wind = E at 6mph
Uintah basin has typical elevation of 5,000 to 5,500 feet or 1524 to 1676.4 m
Pick the highest point but use typical alitude/pressure equation: 101325*(1 - 2.25577*10**-5*h)**5.25588 
"""
#Converting to g/s using ideal gas law
p = 101325*(1 - 2.25577*10**-5*1676.4)**5.25588 #in Pa units
T = -5.6+273.15 #in Kelvin
R = 8.3144621 #in Pa*m^3*K-1*mol-1
M = 16.04 #g/mol molar mass of methane
mass = (p*ans*M)/(T*R)
#print "%f g/s" %mass
print "%f slpm" %(ans*10**3*60)
secans =(ans*10**3*60)
print "%f SCFH" %(ans*10**3*3600/28.32)
#secans =(ans*10**3*3600/28.32)


def widthinator (tapelist, background = 0):
    """
    tapelist must be in the order of A, B, C, and D
    Determine the FWHM of the plume with the highest peak concentration
    """
    maxarray = []
    for i in tapelist:
        maxarray.append(max(i[:,1]))
    winner = argmax(maxarray)
    hmax = ((maxarray[winner]-background)/2.0)
    print hmax
    if winner == 0:
        txtlabel = "TapeA"
        tape = tapelist[0]
    elif winner == 1:
        txtlabel = "TapeB"
        tape = tapelist[1]
    elif winner == 2:
        txtlabel = "TapeC"
        tape = tapelist[2]
    elif winner == 3:
        txtlabel = "TapeD"
        tape = tapelist[3]
    # TODO: 6 pixels
    elif winner == 4:
        txtlabel = "TapeE"
        tape = tapelist[4]
    elif winner == 5:
        txtlabel = "TapeF"
        tape = tapelist[5]
    else:
        print "I don't know what is going on here."
    peakloc = argmax(tape[:,1])
    width = []
    width.append(argsort(abs(tape[:peakloc,1]-hmax-background))[0])
    width.append(argsort(abs(tape[peakloc:,1]-hmax-background))[0]+peakloc)
#    print width
#    pylab.figure(txtlabel)
#    pylab.scatter(tape[width,0], tape[width,1])
#    pylab.scatter(tape[peakloc,0],tape[peakloc,1])
#    pylab.plot(tape[:,0], tape[:,1])
#    pylab.show()
    return tape[width[1],0]-tape[width[0],0]

if fwhmcalc == 'y':
    try:    
        fwhm = widthinator ([tapeArecon, tapeBrecon, tapeCrecon, tapeDrecon, tapeErecon, tapeFrecon], 0.0)
        print "FWHM is %f m" %(fwhm)
    except:
        print "Fail"
        fwhm = -1.0
else:
    fwhm = -1.0


"""
Getting GPS info
"""
peakindex = argmax(pretape[:,1])
gpslat = pretape[peakindex,8]
gpslong = pretape[peakindex,9]

print '\n'


if truncate == 'y':
    print filename[filename.rfind('\\')+1:]+ "\t" +\
    string.ljust(str(carspeed),12)+"\t"+ \
    string.ljust(str(stdcar),12)+"\t"+ \
    string.ljust(str(weightaverage), 12)+"\t"+ \
    string.ljust(str(stdwind), 12)+"\t"+ \
    string.ljust(str(meanwinddir), 12)+"\t"+ \
    string.ljust(str(stdwinddir), 12)+"\t"+ \
    string.ljust(str(peaks),12)+"\t"+ \
    string.ljust(str(startend), 12)+"\t"+ \
    string.ljust(str(background), 12)+"\t"+ \
    '['+"%.2f" %(max(tapeArecon[:,1]+background))+", %.2f" %(max(tapeBrecon[:,1]+background))+", %.2f" %(max(tapeCrecon[:,1])+background)+", %.2f" %(max(tapeDrecon[:,1])+background)+", %.2f" %(max(tapeErecon[:,1])+background)+", %.2f" %(max(tapeFrecon[:,1])+background)+']'+"\t"+ \
    string.ljust(str(a), 12)+"\t"+ \
    string.ljust(str(ans*10**3), 12)+"\t"+\
    string.ljust(str(secans), 12)+"\t"+ \
    string.ljust('['+str(plumestart)+','+str(plumeend)+']',12)+"\t"+ \
    string.ljust(str(gpslat), 12)+"\t"+ \
    string.ljust(str(gpslong), 12)+"\t"+ \
    string.ljust(str(fwhm),12)+"\t"+ \
    string.ljust(str(overheadflag),12)+"\t"+\
    string.ljust(str(lineoverheadratio),12)+"\t"+\
    string.ljust(str(overheadratioc),12)+"\t"+\
    string.ljust(str(startdist), 12)+"\t"+ \
    string.ljust(str(enddist), 12)
    
else:
    print filename[filename.rfind('\\')+1:]+ "\t" +\
    string.ljust(str(carspeed),12)+"\t"+ \
    string.ljust(str(stdcar),12)+"\t"+ \
    string.ljust(str(weightaverage), 12)+"\t"+ \
    string.ljust(str(stdwind), 12)+"\t"+ \
    string.ljust(str(meanwinddir), 12)+"\t"+ \
    string.ljust(str(stdwinddir), 12)+"\t"+ \
    string.ljust(str(peaks),12)+"\t"+ \
    string.ljust(str(startend), 12)+"\t"+ \
    string.ljust(str(background), 12)+"\t"+ \
    '['+"%.2f" %(max(tapeArecon[:,1]+background))+", %.2f" %(max(tapeBrecon[:,1]+background))+", %.2f" %(max(tapeCrecon[:,1])+background)+", %.2f" %(max(tapeDrecon[:,1])+background)+", %.2f" %(max(tapeErecon[:,1])+background)+", %.2f" %(max(tapeFrecon[:,1])+background)+']'+"\t"+ \
    string.ljust(str(a), 12)+"\t"+ \
    string.ljust(str(ans*10**3), 12)+"\t"+\
    string.ljust(str(secans), 12)+"\t"+ \
    string.ljust('['+str(plumestart)+','+str(plumeend)+']',12)+"\t"+ \
    string.ljust(str(gpslat), 12)+"\t"+\
    string.ljust(str(gpslong), 12)+"\t"+\
    string.ljust(str(fwhm),12)+"\t"+\
    string.ljust(str(overheadflag),12)+"\t"+\
    string.ljust(str(lineoverheadratio),12)+"\t"+\
    string.ljust(str(overheadratioc),12)

    

#plot3D = raw_input("Plot 3D? \n")
plot3D ='y'
if plot3D == 'y':
    """Setting plot dimensions"""
    #userDefinedXMin = raw_input("Input YMin \n") #X is horizontal i.e. it is y axis in Tracy convention
    userDefinedXMin ='n'
    if userDefinedXMin.isalpha() == True:
        if truncate =='y':
            minDistance = startdist
        else:
            minDistance = min(tapeAdist[:,0])
    else:
        minDistance = float(userDefinedXMin)
    #userDefinedXMax = raw_input("Input YMax \n")
    userDefinedXMax= 'n'
    if userDefinedXMax.isalpha() == True:
        if truncate =='y':
            maxDistance = enddist
        else:
            maxDistance =  max(tapeAdist[:,0])
        
    else:
        maxDistance = float(userDefinedXMax)

    if maxDistance-minDistance < 100:
        aspect = 2.
    elif maxDistance-minDistance >100 and maxDistance-minDistance <200:
        aspect = 6.
    elif maxDistance-minDistance >200 and maxDistance-minDistance<250:
        aspect = 10.
    elif maxDistance-minDistance > 351:
        aspect = 20.
    else:
        aspect =10.
    
    
    YCoords = []
    """Default heights of the pixels: plume scanner configuration from Uinta trip in February 2013"""
    #z0 = 20*2.54/100.
    #dz = 48*2.54/100.
#    z0 = 1.01 #Australia car was really high
    z0 = 0.401
    dz = 25*2.54/100.
    '''For Plume Scanner Jr'''
#    z0 = 0.432 #Include the bottom point for calculations
#    dz = 27*2.54/100.
    for i in range(6):
        YCoords.append( z0 + i*dz )    
    yPad = 0.001
    maxHeight = YCoords[len(YCoords)-1] + yPad     
    minHeight = YCoords[0] - yPad                             
    nY = 100j      
    nX = 0    
    xInterpFactor = 0.8
    
    def getPoints(tape, YCoordindex):
        """getting all the plot points: selected points that did NOT have
        background subtracted"""
        points = []
        for i in range(len(tape)):
            y = YCoords[YCoordindex]
            points.append( (tape[i,0], y) )
        return points
            
    pointsA = getPoints(tapeAdist, 0)
    pointsB = getPoints(tapeBdist, 1)
    pointsC = getPoints(tapeCdist, 2)
    pointsD = getPoints(tapeDdist, 3)
    pointsD = getPoints(tapeDdist, 3)
    pointsE = getPoints(tapeEdist, 4)
    pointsF = getPoints(tapeFdist, 5)
    points = vstack([pointsA, pointsB, pointsC, pointsD, pointsE, pointsF])
    points = [n for n in points] #convert to list
    values = hstack([tapeAdist[:, 1],
                       tapeBdist[:, 1],
                       tapeCdist[:, 1],
                       tapeDdist[:, 1],
                       tapeEdist[:, 1],
                       tapeFdist[:, 1]])
    values = [v for v in values]    # convert to list
    """Interperlating data for plotting"""
    nInterp = nX
    if nX == 0 and len(YCoords) > 0:
        nInterp = (1j) * int(xInterpFactor*(len(points)/len(YCoords)) )
        nX = nInterp
        nInterp = xInterpFactor*(1j)*int(len(points)/4)
    grid_x, grid_y = mgrid[minDistance:maxDistance:nInterp, minHeight:maxHeight:nY]
    grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
    interpGridX = grid_x
    interpGridY = grid_y
    interpGridZ = grid_z1  
    grid_x = interpGridX
    grid_y = interpGridY
    grid_z = interpGridZ      
    
    """Plotting"""
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 16
    fig = plt.figure(figsize=(20,10), dpi=80)        
    cax = plt.imshow(grid_z.T, extent=(minDistance,maxDistance,minHeight-yPad,maxHeight+yPad), aspect=aspect, origin='lower')
    plt.xlim(minDistance,maxDistance)
    plt.ylim(minHeight,maxHeight)
    points = array(points) #It has to be reconverted into array!!
    values = array(values)
    plt.plot(points[:,0], points[:,1], 'k.', ms=3)
    plt.xlabel('Horizontal Position [m]', fontsize = 'large')
    plt.ylabel('Vertical Position [m]', fontsize = 'large')
    cbar = fig.colorbar(cax, orientation='horizontal')
    cbar.set_label('CH4 Concentration [ppm]', fontsize = 'large')
    plt.draw()
    plt.figtext(0.01, 0.96, filename)
    plt.figtext(0.01, 0.93, "v25sr_piecewise_comp_6px saved "+ time.asctime(time.localtime()))
    plt.figtext(0.01, 0.90, "Carspeed %f m/s with %f m/s" %(carspeed, stdcar))
    plt.figtext(0.01, 0.87, "Windspeed %s m/s with %f m/s and weighted average of %f" %(weightaverage, stdwind, weightaverage))
    plt.figtext(0.01, 0.84, "Background %f ppm" %background)
    plt.figtext(0.01, 0.81, "GPS Lat: %f; GPS Long: %f" %(gpslat, gpslong))    
    plt.figtext(0.01, 0.78, "Flux: %f L/s or %f lpm" %((ans*10**3), secans))
    #print aspect
    #print max(maxarray)
    #savePlot = raw_input("Save plot? \n")
    savePlot = 'n'
    if savePlot == 'y':
        if truncate == 'y':
            fn = filename[:len(filename)-14]+'_truncate_%i_to_%i_m.png' %(startdist, enddist)
        else:
            fn = filename[:len(filename)-14]+'_piecewise_figure.png'
        fig.savefig(fn)
    plt.show()
    
#pdb.set_trace()
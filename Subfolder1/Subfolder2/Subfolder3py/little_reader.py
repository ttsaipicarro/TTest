# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:03:40 2013

@author: ttsai

Takes in a string of the pickled info that you would like to output and retrieve from pickled file

Takes in pickled data file and unpickles it into
0. Time (Epoch) obj["data"]["time"]
1. CH4 Dry obj["data"]["CH4_dry"]
2. H2O obj["data"]["H2O"]
3. Valve mask obj["data"]["ValveMask"]
4. Spectrum ID
5. "WS_WIND_LAT"
6. WS_WIND_LON
7. Carspeed
8. GPS Lat
9. GPS Long
10. CO2
"""
from __future__ import with_statement

import wx
import cPickle
import pdb
import string
import pylab
#import numpy


def printpretty(floatnum):
    return string.ljust(str(floatnum), 12)


def readfunc(filename, peakplot='y', savefileyn='y', outputvar=""):
    time = []
    ch4 = []
    h2o = []
    co2 = []
    carspeed = []
    windlat = []
    gpslat = []
    gpslong = []
    stdwindir = []
    var = []

    # replace .dat extension in input file with _unpickled.txt for the output filename
    savefilename = filename.replace(".dat", "_unpickled.txt")

    # TODO: What if user selected a filename that doesn't have a .dat file extension?

    if savefileyn == 'y':
        f = open(savefilename, "w")
    else:
        pass

    try:
        with open(filename, 'rb') as fileo:
            while True:
                try:
                    obj = cPickle.load(fileo)
                    if obj["data"]["SpectrumID"] == 25:
                        if savefileyn == 'y':
                            f.write(printpretty(obj["data"]["time"]) + "\t" +
                                    printpretty(obj["data"]["CH4_dry"]) + "\t" +
                                    printpretty(obj["data"]["H2O"]) + "\t" +
                                    printpretty(obj["data"]["ValveMask"]) + "\t" +
                                    printpretty(obj["data"]["SpectrumID"]) + "\t")

                            #print printpretty(obj["data"]["time"]), "\t", printpretty(obj["data"]["CH4_dry"])
                            time.append(obj["data"]["time"])
                            ch4.append(obj["data"]["CH4_dry"])
                            h2o.append(obj["data"]["H2O"])

                            try:
                                f.write(printpretty(obj["data"]["WS_WIND_LAT"]) + "\t" +
                                        printpretty(obj["data"]["WS_WIND_LON"]) + "\t" +
                                        printpretty(obj["data"]["CAR_SPEED"]) + "\t" +
                                        printpretty(obj["data"]["GPS_ABS_LAT"]) + "\t" +
                                        printpretty(obj["data"]["GPS_ABS_LONG"]) + "\t")

                                carspeed.append(obj["data"]["CAR_SPEED"])
                                windlat.append(obj["data"]["WS_WIND_LAT"])
                                stdwindir.append(obj["data"]["WIND_DIR_SDEV"])
                                gpslat.append(obj["data"]["GPS_ABS_LAT"])
                                gpslong.append(obj["data"]["GPS_ABS_LONG"])
                                try:
                                    f.write(printpretty(obj["data"]["CO2"])+"\n")
                                    co2.append(obj["data"]["CO2"])
                                    erroro = "End of file reached"
                                except:
                                    f.write("\n")
                                    erroro = 'Missing CO2 data'
                            except:
                                try:
                                    f.write(printpretty(obj["data"]["CO2"])+"\n")
                                    co2.append(obj["data"]["CO2"])
                                    erroro = "Missing Wind and Car speed data"
                                except:
                                    erroro = "Missing Wind and Car speed and CO2 data"
                                try:
                                    var.append(obj["data"][outputvar])
                                except:
                                    var.append('')
                        else:
                            erroro = 'No output'
                            time.append(obj["data"]["time"])
                            ch4.append(obj["data"]["CH4_dry"])
                            h2o.append(obj["data"]["H2O"])
                            co2.append(obj["data"]["CO2"])
                            carspeed.append(obj["data"]["CAR_SPEED"])
                            windlat.append(obj["data"]["WS_WIND_LAT"])
                            stdwindir.append(obj["data"]["WIND_DIR_SDEV"])
                            try:
                                var.append(obj["data"][outputvar])
                            except:
                                var.append('')

                except:
                    print erroro
                    if savefileyn == 'y':
                        f.close()
                    else:
                        pass
                    #pdb.set_trace()
                    if peakplot == 'y':
                        # instead of plotting epoch time directly,
                        # plot seconds since the first data point

                        # if you want epoch time, comment out the
                        # next line

                        # use exact X axis limits
                        # we're letting pylab autoscale Y limits
                        # by not setting them explicitly

                        # Concentrations page
                        pylab.figure("Concentrations")
                        pylab.subplot(3, 1, 1)
                        pylab.plot(time, ch4)
                        pylab.title("Concentrations")
                        pylab.ylabel("CH4 (ppm)")

                        pylab.subplot(3, 1, 2)
                        pylab.plot(time, h2o)
                        pylab.ylabel("H2O")

                        pylab.subplot(3, 1, 3)
                        pylab.plot(time, co2)
                        pylab.xlabel("Epoch Time")
                        pylab.ylabel("CO2")

                        # Peripherals page
                        pylab.figure("Peripherals")  # Assuming that the valves are off for first 50 pts
                        pylab.subplot(2, 1, 1)

                        # I think this is not what was wanted, this only shows the first 50
                        # points and believe we want to skip those
                        #pylab.plot(carspeed[:50])  # only shows the first 50 points
                        #pylab.plot(carspeed[50:])  # shows from point 50 to the end

                        # we'll show a vertical dotted line at 50 points

                        # Car speed
                        pylab.plot(carspeed[:50])

                        pylab.title("Car speed")

                        # Wind lat page
                        pylab.subplot(2, 1, 2)
                        pylab.plot(windlat[:50])
                        pylab.title("Wind Lat")

                        # Std Wind Dev page
                        pylab.figure("Std Wind Dev")
                        pylab.plot(stdwindir[:50])
                        pylab.show()
                    break
    except:
        True
    #print savefilename
    return savefilename, var


class ReadOptionsDialog(wx.Dialog):
    def __init__(self, peakplotState=False, outputState=True, outvarTextVal=''):
        wx.Dialog.__init__(self, None, -1, "Unpickle options", size=(200, 200))

        self.peakPlotCheckbox = wx.CheckBox(self, -1, "Peak plot", (35, 35), (150, 20))
        self.outputCheckbox = wx.CheckBox(self, -1, "File output", (35, 60), (150, 20))
        self.outvarLabel = wx.StaticText(self, -1, "outvar:", (35, 90), (55, 20))
        self.outvarText = wx.TextCtrl(self, -1, "", (80, 87), (80, 20))
        self.okButton = wx.Button(self, wx.ID_OK, "OK", (80, 140), (40, 20))

        # init control settings
        self.peakPlotCheckbox.SetValue(peakplotState)
        self.outputCheckbox.SetValue(outputState)
        self.outvarText.SetValue(outvarTextVal)

        #self.Bind(wx.EVT_BUTTON, self.onClick, self.okButton)

    def onClick(self, event):
        #self.Destroy()
        pass

    def getSettings(self):
        peakPlotState = self.peakPlotCheckbox.GetValue()
        outputState = self.outputCheckbox.GetValue()
        outvarTextVal = self.outvarText.GetValue()

        return peakPlotState, outputState, outvarTextVal


if __name__ == "__main__":
    app = wx.PySimpleApp()

    d = wx.FileDialog(None, "File to unpickle",
                      style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
                      wildcard="dat files (*.dat)|*.dat|All files (*.*)|*.*")

    if d.ShowModal() == wx.ID_OK:
        filename = d.GetPath()
        d.Destroy()

        peakplot = False
        savefileyn = True
        outputvar = ''

        # Pop a dialog with checkboxes for peak plot, output to a file,
        # and value for outputvar (I think it is an optional additional column of data)
        rod = ReadOptionsDialog(peakplotState=peakplot, outputState=savefileyn, outvarTextVal=outputvar)
        result = rod.ShowModal()

        if result == wx.ID_OK:
            peakplot, savefileyn, outputvar = rod.getSettings()
        else:
            print "options dialog canceled"

        # convert options settings to values expected by readfunc
        if peakplot is True:
            peakplot = 'y'
        else:
            peakplot = 'n'

        if savefileyn is True:
            savefileyn = 'y'
        else:
            savefileyn = 'n'

        print "peakplot=", peakplot
        print "savefileyn=", savefileyn
        print "outputvar='%s'" % outputvar

        readfunc(filename, peakplot=peakplot, savefileyn=savefileyn, outputvar=outputvar)

    else:
        d.Destroy()

    #    filename = r'C:\Users\ttsai\Desktop\Uintah Basin Field Campaign\EPA\for paper\batch\plume_1368624776_2013_5_15__6_32.dat'
    #    test = readfunc(filename,'n','n','WIND_DIR_SDEV')
    #    conc =readfunc(filename,'n','n','CH4_dry')
    #    speed = readfunc(filename,'n','n','CAR_SPEED')
    #    pylab.subplot(3,1,1)
    #    pylab.plot(conc[1][:50])
    #    pylab.subplot(3,1,2)
    #    pylab.plot(speed[1][:50])
    #    pylab.subplot(3,1,3)
    #    pylab.plot(test[1][:50])
    #    pylab.show()

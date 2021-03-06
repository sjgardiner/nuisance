from ROOT import *
import os
import sys


gColorList = [kRed, kGreen, kBlue, kYellow, kOrange, kBlack]

def DrawDataMC(keyname, filelist):

    # Extract Data
    data = None
    for readfile in filelist:
        print keyname
        data = readfile[0].Get(keyname)
        if not data: continue
        break

    if not data:
        print "Data not found for : ", keyname
        sys.exit(-1)

    # Main Data Formatting
    data.SetTitle(keyname)
    data.SetLineColor(kBlack)
    data.SetLineWidth(2)

    # Extract MC
    singlemclist   = []
    singledatalist = []
    for i, mcfile in enumerate(allfiles):

        print mcfile[0]
        # Extract individual MC
        mckey = keyname.replace("_data","_MC")
        singlemc = mcfile[0].Get(mckey)
        if singlemc: 

            singlemc = singlemc.Clone(mcfile[1]+"_MC")
            singlemc.SetLineColor( gColorList[i] )
            singlemc.SetLineWidth(2)
            singlemc.SetTitle( mcfile[1] + " (" + str(singlemc.GetTitle().strip()) + ") " )
            
            singlemclist.append(singlemc.Clone())
            del singlemc
            
        # Extra individual data (optional)
        singledata = mcfile[0].Get(keyname)
        if singledata:

            singledata = singledata.Clone(mcfile[1] + "_DATA")
            singledata.SetLineColor( kBlack )
            singledata.SetLineWidth(2)
            singledata.SetTitle( "^-> Saved Data" )
            
            singledatalist.append(singledata.Clone())
            del singledata


    # Assign Ranges
    miny = 99999.9
    maxy = 0.0
    for i in range(data.GetNbinsX()):
        miny = min([data.GetBinContent(i+1) - data.GetBinError(i+1),miny])
        maxy = max([data.GetBinContent(i+1) + data.GetBinError(i+1),maxy])
        for singlemc in singlemclist:
            miny = min([singlemc.GetMinimum(),miny])
            maxy = max([singlemc.GetMaximum(),maxy])            
        for singledata in singledatalist:
            miny = min([singledata.GetMinimum(),miny])
            maxy = max([singledata.GetMaximum(),maxy])
    widthy = maxy - miny

    # Assign Ranges to data
    if "1D" in keyname:    data.GetYaxis().SetRangeUser( miny - 0.1*widthy, maxy + 0.3*widthy)
    elif "2D" in keyname:  data.GetZaxis().SetRangeUser( miny - 0.1*widthy, maxy + 0.3*widthy)

    # Draw Plots 1D
    if "1D" in keyname:   
        data.Draw("E1")
        for mc in singlemclist:
            mc.Draw("SAME HIST")

    # Draw Plots 2D
    elif "2D" in keyname: 
        data.Draw("E1")
        for mc in singlemclist:
            mc.Draw("SAME LEGO")

    # Build Legend
    leg = gPad.BuildLegend(0.45,0.65,0.8,0.85)
    leg.SetFillStyle(0)
    leg.SetFillColorAlpha(0,0.0)
    leg.SetBorderSize(0)

    gStyle.SetOptTitle(1)
    gPad.SetGridx(0)
    gPad.SetGridy(0)
    gPad.Update()

    singlemclist.append(data)
    return singlemclist


if __name__=="__main__":
    c1 = TCanvas("c1","c1",210*4,297*4)
    c1.cd()
    

    # Make filelist
    allfiles = []

    for i in xrange(2, len(sys.argv)):
        print "Reading ", i, sys.argv[i]
        
        # split by comma
        splitname = sys.argv[i].split(",")

        # Get First
        if (os.path.isfile(splitname[0])):
          
            # Get File
            newfile = (TFile(splitname[0],"READ"))
            if not newfile: 
                print "File is not a ROOT file : ", splitname[0]
                sys.exit()
                
            # Set Name
            name = splitname[0].replace(".root","")
            if len(splitname) > 1:
                name = splitname[1]
            
            allfiles.append([newfile, name])


    print allfiles

    # Parse Unique Keys
    uniquekeys = []
    for readfile in allfiles:
        for readkey in readfile[0].GetListOfKeys():
            if not (readkey.GetName().endswith("_data")): continue
            if readkey.GetName() not in uniquekeys: uniquekeys.append(readkey.GetName())

    print uniquekeys
            
    # Setup First Page
    leg = TLegend(0.1,0.1,0.9,0.9)
    for i, readfile in enumerate(allfiles):
        hist = TH1D(readfile[1],readfile[1],1,0,1)
        hist.SetLineColor(gColorList[i % len(gColorList)])
        hist.SetLineWidth(2)

        leg.AddEntry(hist, readfile[1], "l")

    leg.Draw()
    gPad.Update()
    

    outputfile = sys.argv[1]
    c1.Print(outputfile + "(")


    # Loop through unique keys
    for readkey in uniquekeys:
        
        # Draw
        datamclist = DrawDataMC(readkey, allfiles)

        # Save
        c1.Print(outputfile)

    # Now save the legend again to close...
    leg.Draw()
    gPad.Update()
    gPad.Print(outputfile + ")")

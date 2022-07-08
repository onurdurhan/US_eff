import ROOT,os, subprocess
import time
import ctypes
from array import array
import numpy as np
# for fixing a root bug,  will be solved in the forthcoming 6.26 release.
ROOT.gInterpreter.Declare("""
#include "MuFilterHit.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

void fixRoot(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value) {
   std::map<int,float> m = aHit.GetAllSignals();
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(MuFilterHit& aHit, std::vector<TString>& key,std::vector<float>& value) {
   std::map<TString, float> m = aHit.SumOfSignals();
   std::map<TString, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}

void fixRoot(std::vector<genfit::TrackPoint*>& points, std::vector<int>& d,std::vector<int>& k) {
      for(std::size_t i = 0; i < points.size(); ++i) {
        genfit::AbsMeasurement*  m = points[i]->getRawMeasurement();
        d.push_back( m->getDetId() );
        k.push_back( int(m->getHitId()/1000) );
    }
}
""")


Tkey  = ROOT.std.vector('TString')()
Ikey   = ROOT.std.vector('int')()
Value = ROOT.std.vector('float')()
import rootUtils as ut
import shipunit as u
h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=True)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")
options = parser.parse_args()

runNr   = str(options.runNumber).zfill(6)
path = options.path
partitions = 0
#### FROM HERE
if path.find('eos')>0:
    path = os.environ['EOSSHIP']+options.path
    dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path,shell=True) )
    data = "sndsw_raw_"+runNr+".root"
    if  dirlist.find(data)<0:
        dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path+"run_"+runNr,shell=True) )
        while 1>0:
            data = "raw-"+ str(partitions).zfill(4)
            if dirlist.find(data)>0:
                partitions+=1
            else: break

else:
    dirlist  = os.listdir(options.path+"run_"+runNr)
    for x in dirlist:
        data = "raw-"+ str(partitions).zfill(4)
        if x.find(data)>0:
             partitions+=1

#### TO HERE
import SndlhcGeo
if (options.geoFile).find('../')<0: geo = SndlhcGeo.GeoInterface(path+options.geoFile)
else:                                         geo = SndlhcGeo.GeoInterface(options.geoFile[3:])
MuFilter = geo.modules['MuFilter']
Scifi       = geo.modules['Scifi']
nav = ROOT.gGeoManager.GetCurrentNavigator()

A,B,locA,locB,globA,globB    = ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3()
latex = ROOT.TLatex()

nEvents = options.nEvents
# MuFilter mapping of planes and bars 
systemAndPlanes  = {1:2,2:5,3:7}
systemAndBars     = {1:7,2:10,3:60}
def systemAndOrientation(s,plane):
   if s==1 or s==2: return "horizontal"
   if plane%2==1 or plane == 6: return "vertical"
   return "horizontal"

systemAndChannels     = {1:[8,0],2:[6,2],3:[1,0]}
sdict                     = {1:'Veto',2:'US',3:'DS'}

freq      = 160.E6
TDC2ns = 1E9/freq
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# some helper functions


#if options.runNumber>0:
#              f=ROOT.TFile.Open(path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
#              eventTree = f.rawConv
#else:
# for MC data
#              f=ROOT.TFile.Open(options.fname)
#              eventTree = f.cbmsim


if options.runNumber>0:
    eventChain = ROOT.TChain('rawConv')
    if partitions==0:
        eventChain.Add(path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
    else:
        for p in range(partitions):
            eventChain.Add(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')

eventChain.GetEvent(0)
run      = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()
ioman.SetTreeName(eventChain.GetName())
#OT = ioman.GetSink().GetOutTree() #### this newly as well
outFile = ROOT.TMemFile('dummy','CREATE')
source = ROOT.FairFileSource(eventChain.GetCurrentFile())
if partitions>0:
    for p in range(1,partitions):
        source.AddFile(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')

run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

import SndlhcTracking
trackTask = SndlhcTracking.Tracking()
trackTask.SetName('simpleTracking')
run.AddTask(trackTask)

xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()
run.Init()
if partitions>0:  eventTree = ioman.GetInChain()
else:  eventTree = ioman.GetInChain()
#eventTree = ioman.GetInTree()

eventTree.GetEvent(0)
if eventTree.GetBranch('Digi_MuFilterHit'): eventTree.Digi_MuFilterHits = eventTree.Digi_MuFilterHit
OT = ioman.GetSink().GetOutTree() 




def Scifi_track(nPlanes = 8, nClusters = 11):
# check for low occupancy and enough hits in Scifi
    clusters = trackTask.scifiCluster()
    stations = {}
    for s in range(1,6):
       for o in range(2):
          stations[s*10+o] = []
    for cl in clusters:
         detID = cl.GetFirst()
         s  = detID//1000000
         o = (detID//100000)%10
         stations[s*10+o].append(detID)
    nclusters = 0
    check = {}
    for s in range(1,6):
       for o in range(2):
            if len(stations[s*10+o]) > 0: check[s*10+o]=1
            nclusters+=len(stations[s*10+o])
    if len(check)<nPlanes or nclusters > nClusters: return -1
# build trackCandidate
    hitlist = {}
    for k in range(len(clusters)):
           hitlist[k] = clusters[k]
    theTrack = trackTask.fitTrack(hitlist)
    eventTree.ScifiClusters = clusters
    return theTrack


def DS_track():
# check for low occupancy and enough hits in DS stations
    stations = {}
    for s in systemAndPlanes:
       for plane in range(systemAndPlanes[s]): 
          stations[s*10+plane] = {}
    k=-1
    for aHit in eventTree.Digi_MuFilterHits:
         k+=1
         if not aHit.isValid(): continue
         s = aHit.GetDetectorID()//10000
         p = (aHit.GetDetectorID()//1000)%10
         bar = aHit.GetDetectorID()%1000
         plane = s*10+p
         if s==3:
           if bar<60: plane = s*10+2*p
           else:  plane = s*10+2*p+1
         stations[plane][k] = aHit
    if not len(stations[30])*len(stations[31])*len(stations[32])*len(stations[33]) == 1: return -1
# build trackCandidate
    hitlist = {}
    for p in range(30,34):
         k = list(stations[p].keys())[0]
         hitlist[k] = stations[p][k] 
    theTrack = trackTask.fitTrack(hitlist)
    return theTrack

def residual(theTrack, detID, alg=False):
    allignment={0:0.47,1:1.06,2:0.52,3:1.5,4:0.22}
    l  = (detID%10000)//1000
    MuFilter.GetPosition(detID,A,B)
    globA,locA = array('d',[A[0],A[1],A[2]]),array('d',[A[0],A[1],A[2]])
#       if trans2local:   nav.MasterToLocal(globA,locA)
    Z = A[2]
    Y = locA[1]
    if alg== True: Y = Y+ allignment[l]
    rc = extrapolate(theTrack,Z)
    res = Y-rc[1]
    return res

def expected_detID(theTrack, alg = False):
    detIDs={}
    allignment={0:0.47,1:1.06,2:0.52,3:1.5,4:0.22}
    for l in range(5):
        for bar in range(10):
            detID=int(2E4+l*1E3+bar)
            MuFilter.GetPosition(detID,A,B)
            globA,locA = array('d',[A[0],A[1],A[2]]),array('d',[A[0],A[1],A[2]])
#           if trans2local:   nav.MasterToLocal(globA,locA)
            Z = A[2]
            Y = locA[1]
            if alg == True: Y = Y+allignment[l]
            detIDs[(Y,Z)]=detID
    expected_hits = []
    for key in detIDs:
        x,y = extrapolate(theTrack, key[1])
#        print("here is ", key[0]-y,detIDs[key],eventTree.GetReadEvent())
        if abs(key[0]-y) < 3.05 : expected_hits.append(detIDs[key])
    return expected_hits

def extrapolate(theTrack,z_mid):
    fitStatus = theTrack.getFitStatus()
    state = theTrack.getFittedState()
    pos   = state.getPos()
    mom = state.getMom()
    fitStatus = theTrack.getFitStatus()
    slope_x = mom.x()/mom.z()
    slope_y = mom.y()/mom.z()
    x=pos.x()-slope_x*(pos.z()-z_mid)
    y=pos.y()-slope_y*(pos.z()-z_mid)
    optionTrack=options.trackType
    if not optionTrack == 'DS':
        x=pos.x()+slope_x*(z_mid-pos.z())
        y=pos.y()+slope_y*(z_mid-pos.z())
    return x,y

def matchedHits(theTrack,alg=False):
    matchedhits = {}
    allignment={0:0.47,1:1.06,2:0.52,3:1.5,4:0.22}
    for hit in eventTree.Digi_MuFilterHits:
        if not hit.isValid(): continue
        detID = hit.GetDetectorID()
        system = hit.GetSystem()
        l  = (detID%10000)//1000
        bar = detID%1000
        if system == 2: station = (detID-2E4)//1E3
        if not system == 2 :continue
        MuFilter.GetPosition(detID,A,B)
        globA,locA = array('d',[A[0],A[1],A[2]]),array('d',[A[0],A[1],A[2]])
#       if trans2local:   nav.MasterToLocal(globA,locA)
        Z = A[2]
        Y = locA[1]
        if alg==True: Y = Y+allignment[l]
        rc = extrapolate(theTrack,Z)
        res = Y-rc[1]
        matchedhits[hit]=res
    return matchedhits

def std_err(eff,n):
    err = ROOT.TMath.Sqrt(eff*(1-eff)/n)
    return err

from scipy.stats import beta
def binom_int(num, den, confint=0.90):
    quant = (1-confint)/2.
    low = beta.ppf(quant, num, den-num+1)
    high = beta.ppf(1-quant,num+1,den-num)
    return low, high

def check_reach(theTrack):
    reaching  = True
    z_0 = 103.72
    x,y = extrapolate(theTrack,z_0)
#    x_det = [-74.3650,8.16]
    x_det = [-60.,8.16]
    y_det = [-0.8400,57.9950]
    if x < x_det[0] or x > x_det[1]: reaching = False
    if y < y_det[0] or y > y_det[1]: reaching = False
    return reaching


def cris(investigated_layer,all_hits, expected_hits):
    hits_per_layer = {0:[],1:[],2:[],3:[],4:[]}
    expected_bar_layer={0:None,1:None,2:None,3:None,4:None}
    cris = True
    if len(all_hits)>5:cris = False
    if len(expected_hits)<5: cris = False
    for the_hit in all_hits:
        l = (the_hit%10000)//1000
        if l == investigated_layer : continue
        hits_per_layer[l].append(the_hit)
    for expected_hit in expected_hits:
        l = (expected_hit%10000)//1000
        if l == investigated_layer : continue
        expected_bar_layer[l] = expected_hit
    for key in hits_per_layer:
        if key == investigated_layer : continue
        if len(hits_per_layer[key])!=1:cris = False
        if cris == False :break ## break it already !
        if abs(hits_per_layer[key][0]-expected_bar_layer[key])> 1: cris = False
    return cris

def  gaux(x,par):
       #Fit parameters:
       #par[0]= Box_1  (constant)  box to be convoluted by gauss
       #par[1]= gaussian normalization par
       #par[2]= mean
       #par[3]= sigma
       #par[4]= Box_2 (constant) noise
       invsq2pi = 0.3989422804014
       np = 1000.0      # number of convolution steps
       sc =  5.  #5.0      # convolution extends to +-sc Gaussian sigmas
       summe = 0.0
#       xlow = max(0,x[0] - sc * par[3])
       sc = 3. 
       xlow = x[0] - sc# * par[3]
       xupp = x[0] + sc#* par[3]
       step = (xupp-xlow) / np
       i=1.0
#       if par[0]==0 or par[3]==0: return 9999
       while i<=np/2.:
           i+=1
           xx = xlow + (i-.5) * step
           fland = par[0]
           summe += fland*ROOT.TMath.Exp(-0.5*((xx-par[2])/par[3])**2)
#           summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
           xx = xupp - (i-.5) * step
           fland = par[0]
           summe += fland*ROOT.TMath.Exp(-0.5*((xx-par[2])/par[3])**2)
#           summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
       return ( step * par[1]*summe + par[4]) ####/ par[3]) + par[4])


def US_eff(cut, apply_veto=True,apply_alg=False, Nev = nEvents):
    print( eventTree.GetEntries())
    hist = {}
    stations=[0,1,2,3,4]
    bars=[0,1,2,3,4,5,6,7,8,9]
    stations_and_bars={}
    bar_eff = {}
    for l in range(5):
        for bar in range(10):
            detID=int(2E4+l*1E3+bar)
            bar_eff[detID] = [0,0] # n_success, n_failed  
    for station in stations:
        stations_and_bars[station]=bars
    ut.bookHist(hist,"yzhits","y-z hits",100,80.,280.,100,-10.,70.)
    ut.bookHist(hist,"xzhits","x-z hits",100,80.,280.,100,-80.,10.)
    ut.bookHist(hist,"xyhits","x-y hits",100,-80.,10.,100,-10.,70.)
    ut.bookHist(hist,"occupancy","occupancy",10,0.,10.)
    ut.bookHist(hist,"hitmap_US_all","US hits per track all ",100,80.,280.)
    ut.bookHist(hist,"hitmap_US_associated","US hit per track assocated",100,80.,280.)
    ut.bookCanvas(hist,"res_stations_bars", "", 1200,2000, 10 ,5)
    ut.bookCanvas(hist,'res_stations','',700,500,3,2)
    ut.bookCanvas(hist,'avg_res_per_stat','',700,700)
    ut.bookCanvas(hist,'occupancies','',700,700,3,2)
    ut.bookCanvas(hist,'bar_eff','bar eff.',700,500,3,2)
    for station in stations_and_bars:
        for bar in stations_and_bars[station]:
            ut.bookHist(hist,"res_"+"station_"+str(station)+"-bar_"+str(bar),"res_"+"station_"+str(station)+"-bar_"+str(bar),40,-20.,20.)
            hist["res_"+"station_"+str(station)+"-bar_"+str(bar)].GetXaxis().SetTitle("res [cm]")
    for station in stations:
        ut.bookHist(hist,"res_stations-"+str(station),"plane "+str(station),60,-30.,30.)
        hist["res_stations-"+str(station)].GetXaxis().SetTitle("res [cm]")
        ut.bookHist(hist,"occupancy_"+str(station)," nbr of hits ",10,0.,10.)
        hist["occupancy_"+str(station)].GetXaxis().SetTitle("nbr of hits")
    nTr=0
    n_cris=0
    noise_per_stat = {0:[0,0],1:[0,0],2:[0,0],3:[0,0],4:[0,0]}
    background_region = {0:10.,1:10.,2:10.,3:10.,4:5.}

    for event in eventTree:
        if eventTree.GetReadEvent()>Nev:break
        optionTrack=options.trackType
##        if optionTrack=='DS': theTrack = DS_track() this is old
        if optionTrack=='DS': rc = trackTask.ExecuteTask("DS")
        if eventTree.GetReadEvent()>180000: break
#        else: theTrack = Scifi_track()
        if not OT.Reco_MuonTracks.GetEntries()==1: continue #### also this is new
        theTrack = OT.Reco_MuonTracks[0]
        if not hasattr(theTrack,"getFittedState"): continue
        nTr+=1
        print("event at",event.GetReadEvent(), "track",theTrack)
        if not theTrack.getFitStatus().isFitConverged() and optionTrack!='DS':   # for H8 where only two planes / proj were avaiable
            theTrack.Delete()
            continue    
        state = theTrack.getFittedState()
        pos   = state.getPos()
        print("pos track",pos.z())
        mom = state.getMom()
        fitStatus = theTrack.getFitStatus()
        rc = matchedHits(theTrack,apply_alg)
        sorted_keys = sorted(rc, key=lambda dict_key: abs(rc[dict_key]))
        k=1
        nHit_per_stat  = {0:0,1:0,2:0,3:0,4:0}

        reach=check_reach(theTrack)
#        if reach == False: continue
        for hit in rc:
#            hit = sorted_keys[0] ##### don't forget this !
            detID = hit.GetDetectorID()
            system = hit.GetSystem()
            l  = (detID%10000)//1000
            bar = detID%1000
            station = int((detID-2E4)//1E3)
            nHit_per_stat[station]+=1
            hist["res_stations-"+str(station)].Fill(rc[hit])
            noise_per_stat[station][0]+=1
            if abs(rc[hit]) > background_region[station]: noise_per_stat[station][1]+=1
#            hist["res_"+"station_"+str(station)+"-bar_"+str(bar)].Fill(rc[hit])
        for station in nHit_per_stat:
            hist["occupancy_"+str(station)].Fill(nHit_per_stat[station])

        all_hits = []
        for hit in eventTree.Digi_MuFilterHits:
            print("hitttttt",hit)
            if not hit.isValid() : continue
            detID = hit.GetDetectorID()
            system = hit.GetSystem()
            l  = (detID%10000)//1000
            bar = detID%1000
            if system == 2: station = (detID-2E4)//1E3
            if system == 2:
                if hit.GetnSiPMs()!=8: print("Hello")
                print(hit.GetnSides())
            MuFilter.GetPosition(detID,A,B)
            globA,locA = array('d',[A[0],A[1],A[2]]),array('d',[A[0],A[1],A[2]])
#            if trans2local:   nav.MasterToLocal(globA,locA)
            Z = A[2]
            if system == 2: all_hits.append(detID)
            if system == 2: hist["hitmap_US_all"].Fill(Z)
            DSX=False
            if hit.isVertical():
                X = locA[0]
                hist["xzhits"].Fill(Z,X)
                DSX=True
            else:
                Y = locA[1]
                hist["yzhits"].Fill(Z,Y)
            if system == 3 and DSX : hist["xyhits"].Fill(X,Y)
        expected_hits = expected_detID(theTrack, apply_alg)
        if len(expected_hits)<5 : print("a bug is here ", expected_hits)
        for expected_hit in expected_hits :
            l = (expected_hit%10000)//1000
            if expected_hit == 20009 : print(all_hits, expected_hits)
            res  = residual(theTrack, expected_hit, apply_alg)
            cris_condition = cris(l, all_hits, expected_hits)
            if cris_condition == False :continue
            up = expected_hit+1
            down=expected_hit-1
            up_up = up + 1
            down_down = down - 1
            if expected_hit in all_hits and res < background_region[l]:
                bar_eff[expected_hit][0]+=1
            if expected_hit not in all_hits:
                if up not in all_hits and down not in all_hits:
                    if up_up in all_hits and res < background_region[l]:
                        bar_eff[up_up][0]+=1
                    if down_down in all_hits and res<background_region[l]:
                        bar_eff[down_down][0]+=1
                if up in all_hits and res<background_region[l]:
                    bar_eff[up][0]+=1
                if down in all_hits and res<background_region[l]:
                    bar_eff[down][0]+=1        
            if expected_hit not in all_hits and up not in all_hits and down not in all_hits and up_up not in all_hits and down_down not in all_hits:
                bar_eff[expected_hit][1]+=1
    k=1
    l=1
    latex.SetTextColor(ROOT.kRed)
    latex.SetTextSize(0.1)
    stat_eff={}
    box = ROOT.TF1("box","[0]",-15.,15.)
    f_conv = ROOT.TF1Convolution("box","gaus",-5.,5.,True)
    f_conv.SetNofPointsFFT(1000)
    f = ROOT.TF1("f1",f_conv,-15.,15.,f_conv.GetNpar())
    f.SetParameters(0.5,1.,0,1.5)



#    f = ROOT.TF1("f1",gaux,-15.,15.,5)
#    f.SetParameters(.5,1.,0,1.5,1.)

    for station in stations_and_bars:
        np=0
        stat_eff[station]=ROOT.TGraphAsymmErrors()
        stat_eff[station].GetXaxis().SetTitle("bar ID")
        stat_eff[station].GetXaxis().SetLabelSize(0.05)
        stat_eff[station].GetXaxis().SetTitleSize(0.05)
        stat_eff[station].GetXaxis().SetNdivisions(11)
        stat_eff[station].GetYaxis().SetTitle("eff")
        stat_eff[station].GetYaxis().SetTitleSize(0.05)
        stat_eff[station].GetYaxis().SetLabelSize(0.05)
        stat_eff[station].SetMarkerStyle(21)
        for bar in stations_and_bars[station]:
            hist["res_stations_bars"].cd(k)
            key=int(2E4+station*1E3+bar)
            if bar_eff[key][0]==0 and bar_eff[key][1]==0: continue
            eff = bar_eff[key][0]/(bar_eff[key][0]+bar_eff[key][1]) 
            exl,exh = 0.0, 0.0
            lowCl, highCl = binom_int(bar_eff[key][0],bar_eff[key][0]+bar_eff[key][1])
            eyl = abs(lowCl-eff)
            eyh = abs(highCl-eff)
            stat_eff[station].SetPoint(np, bar+1, eff)
            stat_eff[station].SetPointError(np,exl,exh,eyl,eyh)
            txt = str('%.2f' % eff)
            hist["res_"+"station_"+str(station)+"-bar_"+str(bar)].Draw()
            latex.DrawLatex(-10., 10.,"eff = "+txt+" %")
            k+=1
            np+=1
        hist['bar_eff'].cd(l)
        stat_eff[station].Draw('AP')
#        stat_eff[station].Fit('pol1','F')
#        ROOT.gStyle.SetOptFit(1111)
        stat_eff[station].GetHistogram().SetMaximum(1.01)
        stat_eff[station].GetHistogram().SetMinimum(0.9)
        l+=1
    k = 1
    gr = ROOT.TGraphErrors()
    gr.GetXaxis().SetTitle("plane ID")
    gr.GetXaxis().SetNdivisions(6)
    gr.GetYaxis().SetTitle("avg. res [cm]")
    gr.SetMarkerStyle(21)
    np=0
    for station in stations:
        hist["res_stations"].cd(k)
        hist["res_stations-"+str(station)].SetLineWidth(3)
        hist["res_stations-"+str(station)].GetXaxis().SetLabelSize(0.05)
        hist["res_stations-"+str(station)].GetXaxis().SetTitleSize(0.05)
        hist["res_stations-"+str(station)].GetYaxis().SetLabelSize(0.05)
        hist["res_stations-"+str(station)].Fit("f1")
        hist["res_stations-"+str(station)].Draw()
        ROOT.gStyle.SetOptFit(1111)
        hist["occupancies"].cd(k)
        hist["occupancy_"+str(station)].SetLineWidth(3)
        hist["occupancy_"+str(station)].Draw()
        k+=1
        avg = hist["res_stations-"+str(station)].GetMean()
        ex = 0.0
        ey = hist["res_stations-"+str(station)].GetStdDev()
        gr.SetPoint(np,station+1,avg)
        gr.SetPointError(np,ex,ey)
        np+=1
    hist['avg_res_per_stat'].cd()
    gr.Draw('AP')
    f = ROOT.TFile("efficiencies_US_run_"+str(options.runNumber)+"_"+str(cut)+"_cm"+".root", "RECREATE")
    f.cd()
    hist["res_stations"].Write()
    hist["res_stations_bars"].Write()
    hist['avg_res_per_stat'].Write()
    hist['occupancies'].Write()
    hist['hitmap_US_all'].Write()
    hist['hitmap_US_associated'].Write()
    hist['bar_eff'].Write()
    f.Close()
#    ROOT.gROOT.SetBatch(ROOT.kFALSE)
    eff_list = []
    nEntries={0:0,1:0,2:0,3:0,4:0}
    for key in bar_eff:
       if bar_eff[key][0]==0 and bar_eff[key][1]==0:continue
       l= (key%10000)//1000
       nEntries[l]+=bar_eff[key][0]+bar_eff[key][1]
       eff = bar_eff[key][0]*100/(bar_eff[key][0]+bar_eff[key][1])
       print(key,eff, bar_eff[key][0],  bar_eff[key][1])
       eff_list.append(eff)
#    return bar_eff, stat_eff
    print(nEntries)
    print("noise ", noise_per_stat)
    return eff_list

def plot_eff_cuts():
    h={}
    ut.bookCanvas(h,'eff_cuts','',1200,200,2,2)
    ut.bookCanvas(h,'eff_cuts_planes','',700,500)
    f = ROOT.TFile("efficiencies_vs_cuts_"+str(options.runNumber)+".root", "RECREATE")
    
    plane_mg = {1:None,2:None,3:None,4:None}
    for plane in plane_mg:
        plane_mg[plane]=ROOT.TMultiGraph()
    
    eff_graphs = {}
    for l in range(5):
        for bar in range(10):
            detID=int(2E4+l*1E3+bar)
            l = (detID%10000)//1000
            if l == 0 : eff_graphs[detID] = None
            else : eff_graphs[detID] = ROOT.TGraphErrors()
        
    plane_sg = {}
    for l in range(1,5):
        plane_sg[l] = ROOT.TGraphErrors()
    summary_graph = ROOT.TMultiGraph()

    cuts=[3,4,5,6,10]
    for i,cut in enumerate(cuts):
        bar_eff, stat_eff_graphs = US_eff(cut, apply_veto=True, apply_alg = True)
        W, EFF, SIGMA = {}, {}, {}
        for l  in range(1,5):
            W[l] = np.array([])
            EFF[l] = np.array([])
            SIGMA[l] = np.array([])
        for key in bar_eff:
            l = (key%10000)//1000
            eff_i=bar_eff[key][0]/(bar_eff[key][0]+bar_eff[key][1])
            if eff_graphs[key] == None : continue
            ex = 0.0
            ey = std_err(eff_i,bar_eff[key][0]+bar_eff[key][1])
            eff_graphs[key].SetPoint(i,cut,eff_i)
            eff_graphs[key].SetPointError(i,ex,ey)
            if ey == 0: ey = 0.0000000000001
            w_i = 1/ey**2
            W[l] = np.append(W[l],w_i)
            EFF[l] = np.append(EFF[l], eff_i)
            SIGMA[l]=np.append(SIGMA[l],ey)

        for station in stat_eff_graphs:
            fit = stat_eff_graphs[station].GetFunction("pol1")
            print("fit result for ",station, fit.GetParameter(0),"+-", fit.GetParError(0))
#           eff_avg = fit.GetParameter(0)
            ex = 0.0
#            ey = fit.GetParError(0)
            eff_avg =sum(W[station]*EFF[station])/sum(W[station])
            s_sq = sum((EFF[station]-eff_avg)**2)/9
            ey = ROOT.TMath.Sqrt(s_sq)
            print("for amg ", eff_avg,"+-", ey)
            plane_sg[station].SetPoint(i,cut,eff_avg)
            plane_sg[station].SetPointError(i,ex,ey)


    for key in eff_graphs:
        l  = (key%10000)//1000
        if l == 0 : continue
        eff_graphs[key].SetMarkerStyle(21)
        plane_mg[l].Add(eff_graphs[key],"P")
    for l in plane_mg:
        h['eff_cuts'].cd(l)
        plane_mg[l].Draw("AP pmc plc")

    for l in plane_sg:
        plane_sg[l].SetMarkerStyle(21)
        summary_graph.Add(plane_sg[l],"P")

    h['eff_cuts_planes'].cd()
    summary_graph.Draw("AP pmc plc")

    f.cd()
    h['eff_cuts'].Write()
    h['eff_cuts_planes'].Write()
    f.Close()
rc = US_eff(cut=6.,apply_veto = False, apply_alg=True)
#rc = plot_eff_cuts()
#rc2= US_eff(cut=6.,apply_veto = False, apply_alg=False)

#for i in range(len(rc)):
#    print("alligned = "+ str("%.2f" % rc[i]), "no allignment = "+str("%.2f" % rc2[i]))


















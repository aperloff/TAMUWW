import os, sys, math
import numpy as np
from OrderedDict26 import *

def yieldsFromCombinedCard():
    "This function a datacard combining the two lepton channel and makes a table of expected yields"
    leptons = ["electron","muon"]
    variables = ["KinMEBDT"]
    jetBin = [("TwoJ0B","2 Jets"),("ThreeJ0B","3 Jets"),("FourJ0B","$\geqlant$4 Jets")]
    processDict = {
        'Diboson'                               : ['WW','WZ'],
        '\Wjets'                                : ['WJets'],
        '\Zjets'                                : ['ZJets'],
        '\\ttbar'                               : ['TTbar'],
        'Single \cPqt'                          : ['STopS_T','STopS_Tbar','STopT_T','STopT_Tbar','STopTW_T','STopTW_Tbar'],
        'Multijet'                              : ['QCD_ElFULL','QCD_MuFULL'],
        'Total Background'                      : ['WW','WZ','WJets','ZJets','TTbar','STopS_T','STopS_Tbar','STopT_T','STopT_Tbar','STopTW_T','STopTW_Tbar','QCD_ElFULL','QCD_MuFULL'],
        'ggH, \HWW \MH=125\gev'                 : ['ggH125'],
        'qqH, \HWW \MH=125\gev'                 : ['qqH125'],
        'WH\_ZH\_TTH, \HWW \MH=125\gev'         : ['WH_HToWW_M125','ZH_HToWW_M125','TTH_HToWW_M125'],
        'Total \HWW'                            : ['ggH125','qqH125','WH_HToWW_M125','ZH_HToWW_M125','TTH_HToWW_M125'],
        'WH\_ZH\_TTH, \HZZ \MH=125\gev'         : ['WH_HToZZ_M125','ZH_HToZZ_M125','TTH_HToZZ_M125'],
        'WH, \Hbb \MH=125\gev'                  : ['WH125_HToBB'],
        '\\ttH, \Hbb \MH=125\gev'               : ['TTH_HToBB_M125'],
        'Total Volunteer Signal'                : ['WH_HToZZ_M125','ZH_HToZZ_M125','TTH_HToZZ_M125','WH125_HToBB','TTH_HToBB_M125'],
    }
    sort_list=('Diboson','\Wjets','\Zjets','\\ttbar','Single \cPqt','Multijet','Total Background','ggH, \HWW \MH=125\gev','qqH, \HWW \MH=125\gev','WH\_ZH\_TTH, \HWW \MH=125\gev','Total \HWW','WH\_ZH\_TTH, \HZZ \MH=125\gev','WH, \Hbb \MH=125\gev','\\ttH, \Hbb \MH=125\gev','Total Volunteer Signal','Signal\\textsubscript{\HWW}/Bkg','Signal\\textsubscript{\HWW}/$\sqrt{\text{Bkg}}$')

    output_folder = "/uscms/home/aperloff/nobackup/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/datacards/"

    for v in variables:
        rows = []
        rowNames_expanded = []
        rates_combined = {
            'Diboson'                               : [0.0]*len(jetBin),
            '\Wjets'                                : [0.0]*len(jetBin),
            '\Zjets'                                : [0.0]*len(jetBin),
            '\\ttbar'                               : [0.0]*len(jetBin),
            'Single \cPqt'                          : [0.0]*len(jetBin),
            'Multijet'                              : [0.0]*len(jetBin),
            'Total Background'                      : [0.0]*len(jetBin),
            'ggH, \HWW \MH=125\gev'                 : [0.0]*len(jetBin),
            'qqH, \HWW \MH=125\gev'                 : [0.0]*len(jetBin),
            'WH\_ZH\_TTH, \HWW \MH=125\gev'         : [0.0]*len(jetBin),
            'Total \HWW'                            : [0.0]*len(jetBin),
            'WH\_ZH\_TTH, \HZZ \MH=125\gev'         : [0.0]*len(jetBin),
            'WH, \Hbb \MH=125\gev'                  : [0.0]*len(jetBin),
            '\\ttH, \Hbb \MH=125\gev'               : [0.0]*len(jetBin),
            'Total Volunteer Signal'                : [0.0]*len(jetBin),
            'Signal\\textsubscript{\HWW}/Bkg'        : [0.0]*len(jetBin),
            'Signal\\textsubscript{\HWW}/$\sqrt{\text{Bkg}}$' : [0.0]*len(jetBin),
        }
        ijetBin=0
        for jetBinName, jetBinTitle in jetBin:
            print "JetBin =",jetBinName
            process_names_all = []
            process_names = []
            rates = []
            with open(output_folder+"DataCard_"+jetBinName+"_CombinedLepton_"+v+".txt", 'r') as inF:
                for line in inF:
                    if ('process' in line) and ('WH_HToWW_M125' in line):
                        process_names_all = line.split()[1:]
                        #process_names = process_names_all[:len(process_names_all)/2]
                    if ('rate' in line):
                        rates = line.split()[1:]
            for group in processDict:
                beingSummed = ''
                for iproc, proc in enumerate(process_names_all):
                    if proc in processDict[group]:
                        beingSummed += proc+' '
                        rates_combined[group][ijetBin] += float(rates[iproc])
                print '\t',rates_combined[group][ijetBin], "= [", beingSummed, "] =", group
            ijetBin+=1

        for i in range(0,3):
            rates_combined['Signal\\textsubscript{\HWW}/Bkg'][i] = rates_combined['Total \HWW'][i]/rates_combined['Total Background'][i]
            rates_combined['Signal\\textsubscript{\HWW}/$\sqrt{\text{Bkg}}$'][i] = rates_combined['Total \HWW'][i]/math.sqrt(rates_combined['Total Background'][i])


        print rates_combined

        print "\\begin{table}[htbp]\n\\centering"
        print "\\begin{tabular}{lccc} \\hline"
        print "\\textbf{Process} & \\textbf{2 Jets} & \\textbf{3 Jets} & \\textbf{$\geqslant$4 Jets}\\\\ \\hline"
        array = np.ndarray((len(rates_combined), len(rates_combined.values()[0])))
        #for index, val in enumerate(rates_combined.values()): # UNSORTED
        for index, val in enumerate([x[1] for x in sorted(rates_combined.items(), key=lambda pair: sort_list.index(pair[0]))]):
            array[index] = val
            if sort_list[index] in ['Total Background','Total \HWW','Total Volunteer Signal']:
                print '\\rowcolor{mygray}'
            print sort_list[index],'&'," \\\\\n".join([" & ".join(map('{0:.2f}'.format,array[index]))]),'\\\\' #replace str with '{0:.2f}'.format for change precision
        print "\end{tabular}\n\\caption{}\n\\label{tab:yields_"+v+"}"
        print "\end{table}\n\n"

def percentYieldsFromCombinedCard():
    "This function a datacard combining the two lepton channel and makes a table of expected yields"
    leptons = ["electron","muon"]
    variables = ["KinMEBDT"]
    jetBin = [("TwoJ0B","2 Jets"),("ThreeJ0B","3 Jets"),("FourJ0B","$\geqlant$4 Jets")]
    processDict = {
        'Diboson'                               : ['WW','WZ'],
        '\Wjets'                                : ['WJets'],
        '\Zjets'                                : ['ZJets'],
        '\\ttbar'                               : ['TTbar'],
        'Single \cPqt'                          : ['STopS_T','STopS_Tbar','STopT_T','STopT_Tbar','STopTW_T','STopTW_Tbar'],
        'Multijet'                              : ['QCD_ElFULL','QCD_MuFULL'],
        'Total Background'                      : ['WW','WZ','WJets','ZJets','TTbar','STopS_T','STopS_Tbar','STopT_T','STopT_Tbar','STopTW_T','STopTW_Tbar','QCD_ElFULL','QCD_MuFULL'],
        'ggH, \HWW \MH=125\gev'                 : ['ggH125'],
        'qqH, \HWW \MH=125\gev'                 : ['qqH125'],
        'WH\_ZH\_TTH, \HWW \MH=125\gev'         : ['WH_HToWW_M125','ZH_HToWW_M125','TTH_HToWW_M125'],
        'Total \HWW'                            : ['ggH125','qqH125','WH_HToWW_M125','ZH_HToWW_M125','TTH_HToWW_M125'],
        'WH\_ZH\_TTH, \HZZ \MH=125\gev'         : ['WH_HToZZ_M125','ZH_HToZZ_M125','TTH_HToZZ_M125'],
        'WH, \Hbb \MH=125\gev'                  : ['WH125_HToBB'],
        '\\ttH, \Hbb \MH=125\gev'               : ['TTH_HToBB_M125'],
        'Total Volunteer/Total \HWW'            : ['WH_HToZZ_M125','ZH_HToZZ_M125','TTH_HToZZ_M125','WH125_HToBB','TTH_HToBB_M125'],
    }
    sort_list=('Diboson','\Wjets','\Zjets','\\ttbar','Single \cPqt','Multijet','Total Background','ggH, \HWW \MH=125\gev','qqH, \HWW \MH=125\gev','WH\_ZH\_TTH, \HWW \MH=125\gev','Total \HWW','WH\_ZH\_TTH, \HZZ \MH=125\gev','WH, \Hbb \MH=125\gev','\\ttH, \Hbb \MH=125\gev','Total Volunteer/Total \HWW')

    output_folder = "/uscms/home/aperloff/nobackup/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/datacards/"

    for v in variables:
        rows = []
        rowNames_expanded = []
        rates_combined = {
            'Diboson'                               : [0.0]*len(jetBin),
            '\Wjets'                                : [0.0]*len(jetBin),
            '\Zjets'                                : [0.0]*len(jetBin),
            '\\ttbar'                               : [0.0]*len(jetBin),
            'Single \cPqt'                          : [0.0]*len(jetBin),
            'Multijet'                              : [0.0]*len(jetBin),
            'Total Background'                      : [0.0]*len(jetBin),
            'ggH, \HWW \MH=125\gev'                 : [0.0]*len(jetBin),
            'qqH, \HWW \MH=125\gev'                 : [0.0]*len(jetBin),
            'WH\_ZH\_TTH, \HWW \MH=125\gev'         : [0.0]*len(jetBin),
            'Total \HWW'                            : [0.0]*len(jetBin),
            'WH\_ZH\_TTH, \HZZ \MH=125\gev'         : [0.0]*len(jetBin),
            'WH, \Hbb \MH=125\gev'                  : [0.0]*len(jetBin),
            '\\ttH, \Hbb \MH=125\gev'               : [0.0]*len(jetBin),
            'Total Volunteer/Total \HWW'            : [0.0]*len(jetBin),
        }
        ijetBin=0
        for jetBinName, jetBinTitle in jetBin:
            print "JetBin =",jetBinName
            process_names_all = []
            process_names = []
            rates = []
            with open(output_folder+"DataCard_"+jetBinName+"_CombinedLepton_"+v+".txt", 'r') as inF:
                for line in inF:
                    if ('process' in line) and ('WH_HToWW_M125' in line):
                        process_names_all = line.split()[1:]
                        #process_names = process_names_all[:len(process_names_all)/2]
                    if ('rate' in line):
                        rates = line.split()[1:]
            for group in processDict:
                beingSummed = ''
                for iproc, proc in enumerate(process_names_all):
                    if proc in processDict[group]:
                        beingSummed += proc+' '
                        rates_combined[group][ijetBin] += float(rates[iproc])
                print '\t',rates_combined[group][ijetBin], "= [", beingSummed, "] =", group
            ijetBin+=1


        for group in rates_combined:
            for i in range(0,3):
                if group in ['Diboson','\Wjets','\Zjets','\\ttbar','Single \cPqt','Multijet']:
                    rates_combined[group][i]/=rates_combined['Total Background'][i]
                elif group in ['ggH, \HWW \MH=125\gev','qqH, \HWW \MH=125\gev','WH\_ZH\_TTH, \HWW \MH=125\gev','WH\_ZH\_TTH, \HZZ \MH=125\gev','WH, \Hbb \MH=125\gev','\\ttH, \Hbb \MH=125\gev']:
                    rates_combined[group][i]/=rates_combined['Total \HWW'][i]
        for i in range(0,3):
            rates_combined['Total Background'][i]/=rates_combined['Total Background'][i]
            rates_combined['Total Volunteer/Total \HWW'][i]/=rates_combined['Total \HWW'][i]
            rates_combined['Total \HWW'][i]/=rates_combined['Total \HWW'][i]

        print rates_combined

        print "\\begin{table}[htbp]\n\\centering"
        print "\\begin{tabular}{lccc} \\hline"
        print "\\textbf{Process} & \\textbf{2 Jets} & \\textbf{3 Jets} & \\textbf{$\geqslant$4 Jets}\\\\ \\hline"
        array = np.ndarray((len(rates_combined), len(rates_combined.values()[0])))
        for index, val in enumerate([x[1] for x in sorted(rates_combined.items(), key=lambda pair: sort_list.index(pair[0]))]):
            array[index] = val
            if sort_list[index] in ['Total Background','Total \HWW','Total Volunteer/Total \HWW']:
                print '\\rowcolor{mygray}'
            print sort_list[index],'&'," \\\\\n".join([" & ".join(map('{0:.3f}'.format,array[index]))]),'\\\\' #replace str with '{0:.2f}'.format for change precision
        print "\end{tabular}\n\\caption{}\n\\label{tab:yields_"+v+"}"
        print "\end{table}\n\n"

#yieldsFromCombinedCard()
percentYieldsFromCombinedCard()
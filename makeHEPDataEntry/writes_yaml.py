#!/usr/bin/env python
#-*- coding:utf-8 -*-
import os
from hepdata_lib import *
import numpy as np

'''
def conver_hist_1d(hist):
    points = {}
    for key in ["x", "y", "y_error"]:
        points[key] = []

    for binx in range(1,hist.GetNbinsX()+1):
        points["x"].append(hist.GetBinCenter(binx));
        points["y"].append(hist.GetBinContent(binx));
        points["y_error"].append(hist.GetBinError(binx));

    return points;

def conver_hist_2d(hist):
    points = {}
    for key in ["x", "y", "x_edges", "y_edges", "z"]:
        points[key] = []

    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsY()+1):
            if float(hist.GetBinContent(binx,biny)) < 1e-4 : continue;
            if float(hist.GetBinContent(binx,biny)) == 10.0 : continue;
            points["x"].append(hist.GetXaxis().GetBinCenter(binx));
            points["x_edges"].append((hist.GetXaxis().GetBinLowEdge(binx),hist.GetXaxis().GetBinLowEdge(binx+1)));
            points["y"].append(hist.GetYaxis().GetBinCenter(biny));
            points["y_edges"].append((hist.GetYaxis().GetBinLowEdge(biny),hist.GetYaxis().GetBinLowEdge(biny+1)));
            points["z"].append(hist.GetBinContent(binx,biny));

    return points;


def make_table_figure5_left(outidr):
    file = r.TFile("figures/combined_fit_unmasked.root","READ")

    hqcd  = file.Get("shapes_fit_b/ch1_ch4/QCD_GJ");
    hwjet = file.Get("shapes_fit_b/ch1_ch4/VGamma_GJ");
    hvgam = file.Get("shapes_fit_b/ch1_ch4/WJets_GJ");
    hqcd.Add(hwjet);
    hqcd.Add(hvgam);
    hznn = file.Get("shapes_fit_b/ch1_ch4/Znunu");
    hbkg = file.Get("shapes_fit_b/ch1_ch4/total_background");
    hdata = file.Get("shapes_fit_b/ch1_ch4/data");

    points_data = convert_graph_1d(hdata);
    points_bkg  = convert_hist_1d(hbkg);
    points_znn  = convert_hist_1d(hznn);
    points_min  = convert_hist_1d(hqcd);
    
    print points_bkg["x"],"  ",points_bkg["x_edges"],"  ",points_bkg["y"],"   ",points_bkg["y_error"]

    observable = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    observable.values = points_data["x_edges"];
    
    

def make_table_figure10_left(outdir):
    
    file = r.TFile("figures//Figure_010-a.root","READ")
    hobs = file.Get("Observed_limit");
    hexp = file.Get("Expected_limit");

    points_obs = conver_hist_2d(hobs);
    points_exp = conver_hist_2d(hexp);
    
    mmed = Variable("Mediator mass", is_independent=True, is_binned=True, units="GeV")
    mmed.values = points_obs["x_edges"];
    mdm = Variable("Dark matter mass", is_independent=True, is_binned=True, units="GeV")
    mdm.values = points_obs["y_edges"];
    
    obs = Variable("Observed limit", is_independent=False, is_binned=False, units="")
    obs.values = points_obs["z"]

    exp = Variable("Expected limit", is_independent=False, is_binned=False, units="")
    exp.values = points_exp["z"]

    table = Table("DMSimp limits for Vector Mediator")
    table.location = "Data from Figure 10 (left)"
    table.description = "Limit on the signal strength of the DM signal in a simplified model with a vector mediator."
    table.add_variable(mmed)
    table.add_variable(mdm)
    table.add_variable(obs)
    table.add_variable(exp)
    table.add_image("./figures/PDF//Figure_010-a.pdf","./submission/")
    
    return table

'''
def main():
    # Write to this directory
    print ("open main")

    outdir = "./hepdata/"
    submission = Submission()
    #print "##### Generate figure5 left #####"
    submission.add_table(make_table_figure5_left(outdir))
    #submission.add_table(make_table_figure5_right(outdir))
    #submission.add_table(make_table_figure6_top_left(outdir))
    #submission.add_table(make_table_figure6_top_right(outdir))
    #submission.add_table(make_table_figure6_bottom_left(outdir))
    #submission.add_table(make_table_figure6_bottom_right(outdir))
    ##submission.add_table(make_table_figure7_top_left(outdir))
    #submission.add_table(make_table_figure7_top_right(outdir))
    #submission.add_table(make_table_figure7_bottom_left(outdir))
    #submission.add_table(make_table_figure7_bottom_right(outdir))
    #submission.add_table(make_table_figure8_left(outdir))
    #submission.add_table(make_table_figure8_right(outdir))
    #submission.add_table(make_table_figure9_left(outdir))
    #submission.add_table(make_table_figure9_right(outdir))
    #print "##### Generate figure10 left #####"
    #submission.add_table(make_table_figure10_left(outdir))
    #submission.add_table(make_table_figure10_right(outdir))
    #submission.add_table(make_table_figure11_left(outdir))
    #submission.add_table(make_table_figure11_right(outdir))
    #submission.add_table(make_table_figure12_left(outdir))
    #submission.add_table(make_table_figure12_right(outdir))
    #submission.add_table(make_table_figure15(outdir))
    #submission.add_table(make_table_figure16_left(outdir))
    #submission.add_table(make_table_figure16_right(outdir))
    #submission.add_table(make_table_figure17(outdir))
    #submission.add_table(make_table_figure18_right(outdir))

    #submission.read_abstract("./abstract/abstract.txt")
    #submission.create_files(outdir)

if __name__ == '__main__':
    main()

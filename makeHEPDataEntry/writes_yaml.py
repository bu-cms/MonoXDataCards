#!/usr/bin/env python
import os
from hepdata_lib import *
import numpy as np

met_bins_monojet = [250,280,310,340,370,400,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1400]
met_bins_monov   = [250,300,350,400,500,600,750,1000]

################
def convert_hist_1d(hist,width=True):
    points = {}
    for key in ["x","x_edges","y","y_error"]:
        points[key] = []

    for binx in range(1,hist.GetNbinsX()+1):
        points["x"].append(hist.GetBinCenter(binx));
        points["x_edges"].append((hist.GetXaxis().GetBinLowEdge(binx),hist.GetXaxis().GetBinLowEdge(binx+1)));
        if width == True:
            points["y"].append(hist.GetBinContent(binx)*hist.GetBinWidth(binx));
            points["y_error"].append(hist.GetBinError(binx)*hist.GetBinWidth(binx));
        else:
            points["y"].append(hist.GetBinContent(binx));
            points["y_error"].append(hist.GetBinError(binx));

    return points;

################
def convert_hist_2d(hist,isPS=False):
    points = {}
    for key in ["x", "y", "x_edges", "y_edges", "z"]:
        points[key] = []

    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsY()+1):

            if not isPS and float(hist.GetBinContent(binx,biny)) <= 1e-4 : continue;
            if isPS and float(hist.GetBinContent(binx,biny)) <= 3e-3: continue;

            points["x"].append(hist.GetXaxis().GetBinCenter(binx));
            points["x_edges"].append((hist.GetXaxis().GetBinLowEdge(binx),hist.GetXaxis().GetBinLowEdge(binx+1)));
            points["y"].append(hist.GetYaxis().GetBinCenter(biny));
            points["y_edges"].append((hist.GetYaxis().GetBinLowEdge(biny),hist.GetYaxis().GetBinLowEdge(biny+1)));
            points["z"].append(hist.GetBinContent(binx,biny));

    return points;

################
def convert_corr_2d(hist,isMonoV=False):
    points = {}
    for key in ["x", "y", "x_edges", "y_edges", "z"]:
        points[key] = []

    for binx in range(0,hist.GetNbinsX()):
        for biny in range(0,hist.GetNbinsY()):

            if isMonoV == False:
                points["x"].append((met_bins_monojet[binx+1]-met_bins_monojet[binx])/2);
                points["x_edges"].append((met_bins_monojet[binx],met_bins_monojet[binx+1]));
                points["y"].append((met_bins_monojet[biny+1]-met_bins_monojet[biny])/2);
                points["y_edges"].append((met_bins_monojet[biny],met_bins_monojet[biny+1]));
                points["z"].append(hist.GetBinContent(binx+1,biny+1));            
            else:
                points["x"].append((met_bins_monov[binx+1]-met_bins_monov[binx])/2);
                points["x_edges"].append((met_bins_monov[binx],met_bins_monov[binx+1]));
                points["y"].append((met_bins_monov[biny+1]-met_bins_monov[biny])/2);
                points["y_edges"].append((met_bins_monov[biny],met_bins_monov[biny+1]));
                points["z"].append(hist.GetBinContent(binx+1,biny+1));            

    return points;




################
def convert_graph_1d(graph,width=False):
    points = {}
    for key in ["x","y","y_error_up","y_error_dw"]:
        points[key] = []
    for binx in range(0,graph.GetN()):
        x = r.Double()
        y = r.Double() 
        graph.GetPoint(binx,x,y);
        points["x"].append(x);
        if width == False:
            points["y"].append(y);
            points["y_error_up"].append(graph.GetErrorYhigh(binx));
            points["y_error_dw"].append(-1.*graph.GetErrorYlow(binx));
        else:
            points["y"].append(y*(graph.GetErrorXhigh(binx)+graph.GetErrorXlow(binx)));
            points["y_error_up"].append(graph.GetErrorYhigh(binx)*(graph.GetErrorXhigh(binx)+graph.GetErrorXlow(binx)));
            points["y_error_dw"].append(-1.*graph.GetErrorYlow(binx)*(graph.GetErrorXhigh(binx)+graph.GetErrorXlow(binx)));
            
    return points;


##########
def make_table_figure5(outidr,isMonoV):

    file = r.TFile("figures/combined_fit_unmasked.root","READ")
    channel = "";
    if isMonoV == True:
        channel = "ch2_ch4";
    else:
        channel = "ch1_ch4";
    
    hqcd  = file.Get("shapes_fit_b/"+channel+"/QCD_GJ");
    hwjet = file.Get("shapes_fit_b/"+channel+"/VGamma_GJ");
    hvgam = file.Get("shapes_fit_b/"+channel+"/WJets_GJ");
    hqcd.Add(hwjet);
    hqcd.Add(hvgam);
    hznn = file.Get("shapes_fit_b/"+channel+"/Znunu");
    hbkg_post = file.Get("shapes_fit_b/"+channel+"/total_background");
    hbkg_pre  = file.Get("shapes_prefit/"+channel+"/total_background");
    hdata = file.Get("shapes_fit_b/"+channel+"/data");

    points_data = convert_graph_1d(hdata,True);
    points_post = convert_hist_1d(hbkg_post);
    points_pre  = convert_hist_1d(hbkg_pre);
    points_znn  = convert_hist_1d(hznn);
    points_min  = convert_hist_1d(hqcd);
    
    observable = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    observable.values = points_post["x_edges"];

    ### data
    data = Variable("Observed data", is_independent=False, is_binned=False, units="")
    data.values = points_data["y"];
    #data_unc = Uncertainty("data unc.",is_symmetric=False);
    #errors = [];
    #for i in range(0,len(points_data["y_error_dw"])):
    #    errors.append((points_data["y_error_dw"][i],points_data["y_error_up"][i]));
    #data_unc.values = errors;
    #data.uncertainties.append(data_unc);
    
    ### total bkg post-fit
    bkg_postfit = Variable("Total Background post-fit",is_independent=False, is_binned=False, units="")
    bkg_postfit.values = points_post["y"];
    bkg_postfit_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_postfit_unc.values = points_post["y_error"];
    bkg_postfit.uncertainties.append(bkg_postfit_unc);
    
    ### total bkg pre-fit
    bkg_prefit = Variable("Total Background pre-fit",is_independent=False, is_binned=False, units="")
    bkg_prefit.values = points_pre["y"];
    bkg_prefit_unc = Uncertainty("pre-fit unc.",is_symmetric=True);
    bkg_prefit_unc.values = points_pre["y_error"];
    bkg_prefit.uncertainties.append(bkg_prefit_unc);

    ### gamma+jets
    bkg_gamma = Variable("$\gamma$+jets",is_independent=False, is_binned=False, units="")
    bkg_gamma.values = points_znn["y"];
    bkg_gamma_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_gamma_unc.values = points_znn["y_error"];
    bkg_gamma.uncertainties.append(bkg_gamma_unc);

    ### other bkgs
    bkg_other = Variable("Other backgrounds",is_independent=False, is_binned=False, units="")
    bkg_other.values = points_min["y"];
    bkg_other_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_other_unc.values = points_min["y_error"];
    bkg_other.uncertainties.append(bkg_other_unc);

    if isMonoV == True:
        table = Table("Photon CR event yields mono-V category")
        table.location = "Data from Figure 5 (right)"
        table.description = """Comparison between data and MC simulation in the $\gamma$+jets control sample before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. The plot shows the mono-V category. The hadronic recoil $p_{T}$ in $\gamma$+jets events is used as a proxy for $p_{T}^{miss}$ in the signal region. The last bin includes all events with hadronic recoil $p_{T}$ larger than 750 GeV in the mono-V category."""
    else:
        table = Table("Photon CR event yields monojet category")
        table.location = "Data from Figure 5 (left)"
        table.description = """Comparison between data and MC simulation in the $\gamma$+jets control sample before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. The plot shows the monojet category. The hadronic recoil $p_{T}$ in $\gamma$+jets events is used as a proxy for $p_{T}^{miss}$ in the signal region. The last bin includes all events with hadronic recoil $p_{T}$ larger than 1250 GeV in the monojet category."""


    table.add_variable(observable)
    table.add_variable(data)
    table.add_variable(bkg_postfit)
    table.add_variable(bkg_prefit)
    table.add_variable(bkg_gamma)
    table.add_variable(bkg_other)
    if isMonoV == True:
        table.add_image("figures/PDF/Figure_005-b.pdf","./submission/")
    else:
        table.add_image("figures/PDF/Figure_005-a.pdf","./submission/")

    return table

##########
def make_table_figure6_top(outidr,isMonoV):
    file = r.TFile("figures/combined_fit_unmasked.root","READ")

    channel = "";
    if isMonoV == True:
        channel = "ch2_ch2";
    else:
        channel = "ch1_ch2";
    
    hwjet = file.Get("shapes_fit_b/"+channel+"/WJets_ZM");
    hvv   = file.Get("shapes_fit_b/"+channel+"/Dibosons");
    htop  = file.Get("shapes_fit_b/"+channel+"/Top");
    htop.Add(hwjet);
    htop.Add(hvv);
    hznn = file.Get("shapes_fit_b/"+channel+"/Znunu");
    hbkg_post = file.Get("shapes_fit_b/"+channel+"/total_background");
    hbkg_pre  = file.Get("shapes_prefit/"+channel+"/total_background");
    hdata = file.Get("shapes_fit_b/"+channel+"/data");

    points_data = convert_graph_1d(hdata,True);
    points_post = convert_hist_1d(hbkg_post);
    points_pre  = convert_hist_1d(hbkg_pre);
    points_znn  = convert_hist_1d(hznn);
    points_min  = convert_hist_1d(htop);
    
    observable = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    observable.values = points_post["x_edges"];

    ### data
    data = Variable("Observed data", is_independent=False, is_binned=False, units="")
    data.values = points_data["y"];
    #data_unc = Uncertainty("data unc.",is_symmetric=False);
    #errors = [];
    #for i in range(0,len(points_data["y_error_dw"])):
    #    errors.append((points_data["y_error_dw"][i],points_data["y_error_up"][i]));
    #data_unc.values = errors;
    #data.uncertainties.append(data_unc);
    
    ### total bkg post-fit
    bkg_postfit = Variable("Total Background post-fit",is_independent=False, is_binned=False, units="")
    bkg_postfit.values = points_post["y"];
    bkg_postfit_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_postfit_unc.values = points_post["y_error"];
    bkg_postfit.uncertainties.append(bkg_postfit_unc);
    
    ### total bkg pre-fit
    bkg_prefit = Variable("Total Background pre-fit",is_independent=False, is_binned=False, units="")
    bkg_prefit.values = points_pre["y"];
    bkg_prefit_unc = Uncertainty("pre-fit unc.",is_symmetric=True);
    bkg_prefit_unc.values = points_pre["y_error"];
    bkg_prefit.uncertainties.append(bkg_prefit_unc);

    ### zmm+jets
    bkg_zmm = Variable("$Z(mm)$+jets",is_independent=False, is_binned=False, units="")
    bkg_zmm.values = points_znn["y"];
    bkg_zmm_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_zmm_unc.values = points_znn["y_error"];
    bkg_zmm.uncertainties.append(bkg_zmm_unc);

    ### other bkgs
    bkg_other = Variable("Other backgrounds",is_independent=False, is_binned=False, units="")
    bkg_other.values = points_min["y"];
    bkg_other_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_other_unc.values = points_min["y_error"];
    bkg_other.uncertainties.append(bkg_other_unc);

    if isMonoV == True:
        table = Table("Dimuon CR event yields mono-V category")
        table.location = "Data from Figure 6 (top,right)"
        table.description = """Comparison between data and MC simulation in the dimuon control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the mono-V category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, and W+jets processes."""
    else:
        table = Table("Dimuon CR event yields monojet category")
        table.location = "Data from Figure 6 (top,left)"
        table.description = """Comparison between data and MC simulation in the dimuon control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the monojet category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, and W+jets processes."""


    table.add_variable(observable)
    table.add_variable(data)
    table.add_variable(bkg_postfit)
    table.add_variable(bkg_prefit)
    table.add_variable(bkg_zmm)
    table.add_variable(bkg_other)
    if isMonoV == True:
        table.add_image("figures/PDF/Figure_006-b.pdf","./submission/")
    else:
        table.add_image("figures/PDF/Figure_006-a.pdf","./submission/")

    return table


##########
def make_table_figure6_bottom(outidr,isMonoV):
    file = r.TFile("figures/combined_fit_unmasked.root","READ")

    channel = "";
    if isMonoV == True:
        channel = "ch2_ch5";
    else:
        channel = "ch1_ch5";
    
    hznn = file.Get("shapes_fit_b/"+channel+"/WJets_ZE");
    hvv   = file.Get("shapes_fit_b/"+channel+"/Dibosons");
    htop  = file.Get("shapes_fit_b/"+channel+"/Top");
    htop.Add(hznn);
    htop.Add(hvv);
    hznn = file.Get("shapes_fit_b/"+channel+"/Znunu");
    hbkg_post = file.Get("shapes_fit_b/"+channel+"/total_background");
    hbkg_pre  = file.Get("shapes_prefit/"+channel+"/total_background");
    hdata = file.Get("shapes_fit_b/"+channel+"/data");

    points_data = convert_graph_1d(hdata,True);
    points_post = convert_hist_1d(hbkg_post);
    points_pre  = convert_hist_1d(hbkg_pre);
    points_znn  = convert_hist_1d(hznn);
    points_min  = convert_hist_1d(htop);
    
    observable = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    observable.values = points_post["x_edges"];

    ### data
    data = Variable("Observed data", is_independent=False, is_binned=False, units="")
    data.values = points_data["y"];
    #data_unc = Uncertainty("Data unc.",is_symmetric=False);
    #errors = [];
    #for i in range(0,len(points_data["y_error_dw"])):
    #    errors.append((points_data["y_error_dw"][i],points_data["y_error_up"][i]));
    #data_unc.values = errors;
    #data.uncertainties.append(data_unc);
    
    ### total bkg post-fit
    bkg_postfit = Variable("Total Background post-fit",is_independent=False, is_binned=False, units="")
    bkg_postfit.values = points_post["y"];
    bkg_postfit_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_postfit_unc.values = points_post["y_error"];
    bkg_postfit.uncertainties.append(bkg_postfit_unc);
    
    ### total bkg pre-fit
    bkg_prefit = Variable("Total Background pre-fit",is_independent=False, is_binned=False, units="")
    bkg_prefit.values = points_pre["y"];
    bkg_prefit_unc = Uncertainty("pre-fit unc.",is_symmetric=True);
    bkg_prefit_unc.values = points_pre["y_error"];
    bkg_prefit.uncertainties.append(bkg_prefit_unc);

    ### zee+jets
    bkg_zee = Variable("$Z(ee)$+jets",is_independent=False, is_binned=False, units="")
    bkg_zee.values = points_znn["y"];
    bkg_zee_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_zee_unc.values = points_znn["y_error"];
    bkg_zee.uncertainties.append(bkg_zee_unc);

    ### other bkgs
    bkg_other = Variable("Other backgrounds",is_independent=False, is_binned=False, units="")
    bkg_other.values = points_min["y"];
    bkg_other_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_other_unc.values = points_min["y_error"];
    bkg_other.uncertainties.append(bkg_other_unc);

    if isMonoV == True:
        table = Table("Dielectron CR event yields mono-V category")
        table.location = "Data from Figure 6 (bottom,right)"
        table.description = """Comparison between data and MC simulation in the dielectron control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the mono-V category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, and W+jets processes."""
    else:
        table = Table("Dielectron CR event yields monojet category")
        table.location = "Data from Figure 6 (bottom,left)"
        table.description = """Comparison between data and MC simulation in the dielectron control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the monojet category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, and W+jets processes."""


    table.add_variable(observable)
    table.add_variable(data)
    table.add_variable(bkg_postfit)
    table.add_variable(bkg_prefit)
    table.add_variable(bkg_zee)
    table.add_variable(bkg_other)
    if isMonoV == True:
        table.add_image("figures/PDF/Figure_006-d.pdf","./submission/")
    else:
        table.add_image("figures/PDF/Figure_006-c.pdf","./submission/")

    return table


##########
def make_table_figure7_top(outidr,isMonoV):

    file = r.TFile("figures/combined_fit_unmasked.root","READ")
    channel = "";
    if isMonoV == True:
        channel = "ch2_ch3";
    else:
        channel = "ch1_ch3";
    
    hzjet = file.Get("shapes_fit_b/"+channel+"/ZJets_WM");
    hvv   = file.Get("shapes_fit_b/"+channel+"/Dibosons");
    htop  = file.Get("shapes_fit_b/"+channel+"/Top");
    hqcd  = file.Get("shapes_fit_b/"+channel+"/QCD_WM");
    htop.Add(hzjet);
    htop.Add(hvv);
    htop.Add(hqcd);
    hwjet = file.Get("shapes_fit_b/"+channel+"/WJets");
    hbkg_post = file.Get("shapes_fit_b/"+channel+"/total_background");
    hbkg_pre  = file.Get("shapes_prefit/"+channel+"/total_background");
    hdata = file.Get("shapes_fit_b/"+channel+"/data");

    points_data = convert_graph_1d(hdata,True);
    points_post = convert_hist_1d(hbkg_post);
    points_pre  = convert_hist_1d(hbkg_pre);
    points_wjet  = convert_hist_1d(hwjet);
    points_min  = convert_hist_1d(htop);
    
    observable = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    observable.values = points_post["x_edges"];

    ### data
    data = Variable("Observed data", is_independent=False, is_binned=False, units="")
    data.values = points_data["y"];
    #data_unc = Uncertainty("Data unc.",is_symmetric=False);
    #errors = [];
    #for i in range(0,len(points_data["y_error_dw"])):
    #    errors.append((points_data["y_error_dw"][i],points_data["y_error_up"][i]));
    #data_unc.values = errors;
    #data.uncertainties.append(data_unc);
    
    ### total bkg post-fit
    bkg_postfit = Variable("Total Background post-fit",is_independent=False, is_binned=False, units="")
    bkg_postfit.values = points_post["y"];
    bkg_postfit_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_postfit_unc.values = points_post["y_error"];
    bkg_postfit.uncertainties.append(bkg_postfit_unc);
    
    ### total bkg pre-fit
    bkg_prefit = Variable("Total Background pre-fit",is_independent=False, is_binned=False, units="")
    bkg_prefit.values = points_pre["y"];
    bkg_prefit_unc = Uncertainty("pre-fit unc.",is_symmetric=True);
    bkg_prefit_unc.values = points_pre["y_error"];
    bkg_prefit.uncertainties.append(bkg_prefit_unc);

    ### wjet+jets
    bkg_wjet = Variable("$W(mn)$+jets",is_independent=False, is_binned=False, units="")
    bkg_wjet.values = points_wjet["y"];
    bkg_wjet_unc = Uncertainty("postfit unc.",is_symmetric=True);
    bkg_wjet_unc.values = points_wjet["y_error"];
    bkg_wjet.uncertainties.append(bkg_wjet_unc);

    ### other bkgs
    bkg_other = Variable("Other backgrounds",is_independent=False, is_binned=False, units="")
    bkg_other.values = points_min["y"];
    bkg_other_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_other_unc.values = points_min["y_error"];
    bkg_other.uncertainties.append(bkg_other_unc);

    if isMonoV == True:
        table = Table("Single-muon CR event yields for the mono-V category")
        table.location = "Data from Figure 7 (top,right)"
        table.description = """Comparison between data and MC simulation in the single-muon control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the mono-V category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, Z+jets, and QCD multijet processes."""
    else:
        table = Table("Single-muon CR event yields for the monojet category")
        table.location = "Data from Figure 7 (top,left)"
        table.description = """Comparison between data and MC simulation in the single-muon control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the monojet category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, Z+jets, and QCD multijet processes."""


    table.add_variable(observable)
    table.add_variable(data)
    table.add_variable(bkg_postfit)
    table.add_variable(bkg_prefit)
    table.add_variable(bkg_wjet)
    table.add_variable(bkg_other)
    if isMonoV == True:
        table.add_image("figures/PDF/Figure_007-b.pdf","./submission/")
    else:
        table.add_image("figures/PDF/Figure_007-a.pdf","./submission/")

    return table


##########
def make_table_figure7_bottom(outidr,isMonoV):

    file = r.TFile("figures/combined_fit_unmasked.root","READ")

    channel = "";
    if isMonoV == True:
        channel = "ch2_ch6";
    else:
        channel = "ch1_ch6";
    
    hzjet = file.Get("shapes_fit_b/"+channel+"/ZJets_WE");
    hvv   = file.Get("shapes_fit_b/"+channel+"/Dibosons");
    htop  = file.Get("shapes_fit_b/"+channel+"/Top");
    hqcd  = file.Get("shapes_fit_b/"+channel+"/QCD_WE");
    htop.Add(hzjet);
    htop.Add(hvv);
    htop.Add(hqcd);
    hwjet = file.Get("shapes_fit_b/"+channel+"/WJets");
    hbkg_post = file.Get("shapes_fit_b/"+channel+"/total_background");
    hbkg_pre  = file.Get("shapes_prefit/"+channel+"/total_background");
    hdata = file.Get("shapes_fit_b/"+channel+"/data");

    points_data = convert_graph_1d(hdata,True);
    points_post = convert_hist_1d(hbkg_post);
    points_pre  = convert_hist_1d(hbkg_pre);
    points_wjet  = convert_hist_1d(hwjet);
    points_min  = convert_hist_1d(htop);
    
    observable = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    observable.values = points_post["x_edges"];

    ### data
    data = Variable("Observed data", is_independent=False, is_binned=False, units="")
    data.values = points_data["y"];
    #data_unc = Uncertainty("Data unc.",is_symmetric=False);
    #errors = [];
    #for i in range(0,len(points_data["y_error_dw"])):
    #    errors.append((points_data["y_error_dw"][i],points_data["y_error_up"][i]));
    #data_unc.values = errors;
    #data.uncertainties.append(data_unc);
    
    ### total bkg post-fit
    bkg_postfit = Variable("Total Background post-fit",is_independent=False, is_binned=False, units="")
    bkg_postfit.values = points_post["y"];
    bkg_postfit_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_postfit_unc.values = points_post["y_error"];
    bkg_postfit.uncertainties.append(bkg_postfit_unc);
    
    ### total bkg pre-fit
    bkg_prefit = Variable("Total Background pre-fit",is_independent=False, is_binned=False, units="")
    bkg_prefit.values = points_pre["y"];
    bkg_prefit_unc = Uncertainty("pre-fit unc.",is_symmetric=True);
    bkg_prefit_unc.values = points_pre["y_error"];
    bkg_prefit.uncertainties.append(bkg_prefit_unc);

    ### zee+jets
    bkg_zee = Variable("$W(en)$+jets",is_independent=False, is_binned=False, units="")
    bkg_zee.values = points_wjet["y"];
    bkg_zee_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_zee_unc.values = points_wjet["y_error"];
    bkg_zee.uncertainties.append(bkg_zee_unc);

    ### other bkgs
    bkg_other = Variable("Other backgrounds",is_independent=False, is_binned=False, units="")
    bkg_other.values = points_min["y"];
    bkg_other_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_other_unc.values = points_min["y_error"];
    bkg_other.uncertainties.append(bkg_other_unc);

    if isMonoV == True:
        table = Table("Single-electron CR event yields for the mono-V category")
        table.location = "Data from Figure 7 (bottom,right)"
        table.description = """Comparison between data and MC simulation in the single-electron control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the mono-V category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, Z+jets, and QCD multijet processes."""
    else:
        table = Table("Single-electron CR event yields for the monojet category")
        table.location = "Data from Figure 7 (bottom,left)"
        table.description = """Comparison between data and MC simulation in the single-electron control samples before and after performing the simultaneous fit across all the control samples and the signal region assuming the absence of any signal. Plot correspond to the monojet category. The hadronic recoil $p_{T}$ in dilepton events is used as a proxy for $p_{T}^{miss}$ in the signal region. The other backgrounds include top quark, diboson, Z+jets, and QCD multijet processes."""


    table.add_variable(observable)
    table.add_variable(data)
    table.add_variable(bkg_postfit)
    table.add_variable(bkg_prefit)
    table.add_variable(bkg_zee)
    table.add_variable(bkg_other)
    if isMonoV == True:
        table.add_image("figures/PDF/Figure_007-d.pdf","./submission/")
    else:
        table.add_image("figures/PDF/Figure_007-c.pdf","./submission/")

    return table

##########
def make_table_figure8_and_9(outidr,isMonoV,isMasked):

    if isMasked:
        file = r.TFile("figures/combined_fit_masked.root","READ")
    else:
        file = r.TFile("figures/combined_fit_unmasked.root","READ")
        
    channel = "";
    if isMonoV == True:
        channel = "ch2_ch1";
    else:
        channel = "ch1_ch1";
    
    hzjet = file.Get("shapes_fit_b/"+channel+"/ZJets");
    hvv   = file.Get("shapes_fit_b/"+channel+"/Dibosons");
    htop  = file.Get("shapes_fit_b/"+channel+"/Top");
    hqcd  = file.Get("shapes_fit_b/"+channel+"/QCD");
    hgam  = file.Get("shapes_fit_b/"+channel+"/GJets");
    hgam.Add(hzjet);
    hwjet = file.Get("shapes_fit_b/"+channel+"/WJets");
    hznn  = file.Get("shapes_fit_b/"+channel+"/Znunu");
    hbkg_post = file.Get("shapes_fit_b/"+channel+"/total_background");
    hbkg_pre  = file.Get("shapes_prefit/"+channel+"/total_background");
    hdata = file.Get("shapes_fit_b/"+channel+"/data");

    if isMonoV == True:
        file_mj_av = r.TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoJ_801_0.25_catmonov_13TeV_v1.root","READ");
        file_mw_av = r.TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoW_801_0.25_catmonov_13TeV_v1.root","READ");
        file_mz_av = r.TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoZ_801_0.25_catmonov_13TeV_v1.root","READ");
    else:
        file_mj_av = r.TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoJ_801_0.25_catmonojet_13TeV_v1.root","READ");
        file_mw_av = r.TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoW_801_0.25_catmonojet_13TeV_v1.root","READ");
        file_mz_av = r.TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoZ_801_0.25_catmonojet_13TeV_v1.root","READ");

    hmj_av = file_mj_av.FindObjectAny("signal_signal_80120000001");
    hmw_av = file_mw_av.FindObjectAny("signal_signal_80120000001");
    hmz_av = file_mz_av.FindObjectAny("signal_signal_80120000001");
    hmj_av.Add(hmw_av);
    hmj_av.Add(hmz_av);
    hmj_av.Scale(1.,"width");
    
    points_data = convert_graph_1d(hdata,True);
    points_post = convert_hist_1d(hbkg_post);
    points_pre  = convert_hist_1d(hbkg_pre);
    points_wjet = convert_hist_1d(hwjet);
    points_znn  = convert_hist_1d(hznn);
    points_top  = convert_hist_1d(htop);
    points_vv  = convert_hist_1d(hvv);
    points_gam  = convert_hist_1d(hgam);
    points_qcd  = convert_hist_1d(hqcd);
    points_av = convert_hist_1d(hmj_av);

    observable = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    observable.values = points_post["x_edges"];

    ### data
    data = Variable("Observed data", is_independent=False, is_binned=False, units="")
    data.values = points_data["y"];
    #data_unc = Uncertainty("Data unc.",is_symmetric=False);
    #errors = [];
    #for i in range(0,len(points_data["y_error_dw"])):
    #    errors.append((points_data["y_error_dw"][i],points_data["y_error_up"][i]));
    #data_unc.values = errors;
    #data.uncertainties.append(data_unc);

    ### signal monojet
    dmsignal = Variable("DM signal Axial-Vector", is_independent=False, is_binned=False, units="")
    dmsignal.values = points_av["y"];
    #dmsignal_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    #dmsignal_unc.values = points_av["y_error"];
    #dmsignal.uncertainties.append(dmsignal_unc);

    ### total bkg post-fit
    bkg_postfit = Variable("Total Background post-fit",is_independent=False, is_binned=False, units="")
    bkg_postfit.values = points_post["y"];
    bkg_postfit_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_postfit_unc.values = points_post["y_error"];
    bkg_postfit.uncertainties.append(bkg_postfit_unc);
    
    ### total bkg pre-fit
    bkg_prefit = Variable("Total Background pre-fit",is_independent=False, is_binned=False, units="")
    bkg_prefit.values = points_pre["y"];
    bkg_prefit_unc = Uncertainty("pre-fit unc.",is_symmetric=True);
    bkg_prefit_unc.values = points_pre["y_error"];
    bkg_prefit.uncertainties.append(bkg_prefit_unc);

    ### znn+jets
    bkg_znn = Variable("$Z(nn)$+jets",is_independent=False, is_binned=False, units="")
    bkg_znn.values = points_wjet["y"];
    bkg_znn_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_znn_unc.values = points_wjet["y_error"];
    bkg_znn.uncertainties.append(bkg_znn_unc);

    ### wjet+jets
    bkg_wjet = Variable("$W(ln)$+jets",is_independent=False, is_binned=False, units="")
    bkg_wjet.values = points_wjet["y"];
    bkg_wjet_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_wjet_unc.values = points_wjet["y_error"];
    bkg_wjet.uncertainties.append(bkg_wjet_unc);

    ### Diboson
    bkg_vv = Variable("Dibosons",is_independent=False, is_binned=False, units="")
    bkg_vv.values = points_vv["y"];
    bkg_vv_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_vv_unc.values = points_vv["y_error"];
    bkg_vv.uncertainties.append(bkg_vv_unc);

    ### Top
    bkg_top = Variable("Top quark",is_independent=False, is_binned=False, units="")
    bkg_top.values = points_top["y"];
    bkg_top_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_top_unc.values = points_top["y_error"];
    bkg_top.uncertainties.append(bkg_top_unc);

    ### Zll+gamma+jets
    bkg_gam = Variable("$Z(ll)$+jets+$\gamma$+jets",is_independent=False, is_binned=False, units="")
    bkg_gam.values = points_gam["y"];
    bkg_gam_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_gam_unc.values = points_gam["y_error"];
    bkg_gam.uncertainties.append(bkg_gam_unc);

    ### QCD
    bkg_qcd = Variable("QCD multijets",is_independent=False, is_binned=False, units="")
    bkg_qcd.values = points_qcd["y"];
    bkg_qcd_unc = Uncertainty("post-fit unc.",is_symmetric=True);
    bkg_qcd_unc.values = points_qcd["y_error"];
    bkg_qcd.uncertainties.append(bkg_qcd_unc);

    if isMonoV == True and isMasked == True:
        table = Table("Signal region yields for the mono-V category from CR-only fit")
        table.location = "Data from Figure 8 (right)"
        table.description = """Observed $p_{T}^{miss}$ distribution in the mono-V signal region compared with the post-fit background expectations for various SM processes. The last bin includes all events with $p_{T}^{miss} > 750$ GeV for the mono-V category. The expected background distributions are evaluated after performing a combined fit to the data in all the control samples, not including the signal region. Expected signal distribution for a 2 TeV axial-vector mediator, decaying to 1 GeV DM particles, is overlaid."""
    elif isMonoV == False and isMasked == True:
        table = Table("Signal region yields for the monojet category from CR-only fit")
        table.location = "Data from Figure 8 (left)"
        table.description = """Observed $p_{T}^{miss}$ distribution in the monojet signal region compared with the post-fit background expectations for various SM processes. The last bin includes all events with $p_{T}^{miss} > 1250$ GeV for the monojet category. The expected background distributions are evaluated after performing a combined fit to the data in all the control samples, not including the signal region. Expected signal distribution for a 2 TeV axial-vector mediator, decaying to 1 GeV DM particles, is overlaid."""

    if isMonoV == True and isMasked == False:
        table = Table("Signal region yields for the mono-V category from b-only fit")
        table.location = "Data from Figure 9 (right)"
        table.description = """Observed $p_{T}^{miss}$ distribution in the mono-V signal region compared with the post-fit background expectations for various SM processes. The last bin includes all events with $p_{T}^{miss} > 750$ GeV for the mono-V category. The expected background distributions are evaluated after performing a combined fit to the data in all the control samples, not including as well as in the signal region. The fit is performed assuming the absence of any signal. Expected signal distribution for a 2 TeV axial-vector mediator, decaying to 1 GeV DM particles, is overlaid."""
    elif isMonoV == False and isMasked == False:
        table = Table("Signal region yields for the monojet category from b-only fit")
        table.location = "Data from Figure 9 (left)"
        table.description = """Observed $p_{T}^{miss}$ distribution in the monojet signal region compared with the post-fit background expectations for various SM processes. The last bin includes all events with $p_{T}^{miss} > 1250$ GeV for the monojet category. The expected background distributions are evaluated after performing a combined fit to the data in all the control samples, as well as in the signal region.  The fit is performed assuming the absence of any signal. Expected signal distribution for a 2 TeV axial-vector mediator, decaying to 1 GeV DM particles, is overlaid."""

    table.add_variable(observable)
    table.add_variable(data)
    table.add_variable(dmsignal)
    table.add_variable(bkg_postfit)
    table.add_variable(bkg_prefit)
    table.add_variable(bkg_znn)
    table.add_variable(bkg_wjet)
    table.add_variable(bkg_vv)
    table.add_variable(bkg_top)
    table.add_variable(bkg_gam)
    table.add_variable(bkg_qcd)

    if isMonoV == True and isMasked == True:
        table.add_image("figures/PDF/Figure_008-b.pdf","./submission/")
    elif isMonoV == False and isMasked == True:
        table.add_image("figures/PDF/Figure_008-a.pdf","./submission/")

    if isMonoV == True and isMasked == False:
        table.add_image("figures/PDF/Figure_009-b.pdf","./submission/")
    elif isMonoV == False and isMasked == False:
        table.add_image("figures/PDF/Figure_009-a.pdf","./submission/")
    return table
    

################
def make_table_figure10(outdir,isAV=False,type_plot="observed"):

    if isAV == False:
        file = r.TFile("figures//Figure_010-a.root","READ")
    else:
        file = r.TFile("figures//Figure_010-b.root","READ")

    if type_plot == "observed":
        hobs = file.Get("Observed_limit");
    elif type_plot == "expected":
        hobs = file.Get("Expected_limit");

    points_obs = convert_hist_2d(hobs);
    
    mmed = Variable("Mediator mass", is_independent=True, is_binned=True, units="GeV")
    mmed.values = points_obs["x_edges"];

    mdm = Variable("Dark matter mass", is_independent=True, is_binned=True, units="GeV")
    mdm.values = points_obs["y_edges"];
 
    if type_plot == "observed":
        obs = Variable("Observed limit", is_independent=False, is_binned=False, units="")
        obs.values = points_obs["z"]
    elif type_plot == "expected":
        obs = Variable("Expected limit", is_independent=False, is_binned=False, units="")
        obs.values = points_obs["z"]
    
    if isAV == False:

        if type_plot == "observed":
            table = Table("Observed limits vs mMED-mDM for a vector mediator")
        elif type_plot == "expected":
            table = Table("Expected limits vs mMED-mDM for a vector mediator")

        table.location = "Data from Figure 10 (left)"

        if type_plot == "observed":
            table.description = "Observed exclusion limits at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a vector mediator. N.B.: the granularity in the presented limit values is reduced with respect to the published result."
        elif type_plot == "expected":
            table.description = "Expected exclusion limits at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a vector mediator. N.B.: the granularity in the presented limit values is reduced with respect to the published result."
    else:

        if type_plot == "observed":
            table = Table("Observed limits vs mMED-mDM for a axial-vector mediator")
        elif type_plot == "expected":
            table = Table("Expected limits vs mMED-mDM for a axial-vector mediator")
            
        table.location = "Data from Figure 10 (left)"

        if type_plot == "observed":
            table.description = "Observed exclusion limits at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming an axial-vector mediator. N.B.: the granularity in the presented limit values is reduced with respect to the published result."
        elif type_plot == "expected":
            table.description = "Expected exclusion limits at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming an axial-vector mediator. N.B.: the granularity in the presented limit values is reduced with respect to the published result."

    table.add_variable(mmed)
    table.add_variable(mdm)
    table.add_variable(obs)

    if isAV == False:
        table.add_image("./figures/PDF//Figure_010-a.pdf","./submission/")
    else:
        table.add_image("./figures/PDF//Figure_010-b.pdf","./submission/")

    return table

################
def make_table_figure10_cont(outdir,isAV=False,type_plot="observed"):

    if isAV == False:
        file = r.TFile("figures//Figure_010-a_contour.root","READ")
    else:
        file = r.TFile("figures//Figure_010-b_contour.root","READ")

    if type_plot == "observed" :
        hobs = file.Get("Observed_contour");
    elif type_plot == "expected" :    
        hobs = file.Get("Expected_contour");

    points_obs = convert_graph_1d(hobs);

    mmed = Variable("Mediator mass", is_independent=True, is_binned=False, units="GeV")
    mmed.values = points_obs["x"];

    mdm = Variable("Dark matter mass", is_independent=False, is_binned=False, units="GeV")
    mdm.values = points_obs["y"];

    if not isAV:
        if type_plot == "observed" :
            table = Table("Observed exclusion contour for vector mediator")
            table.location = "Data from Figure 10 (left)"
            table.description = "Observed exclusion contour at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a vector mediator."
        elif type_plot == "expected" :
            table = Table("Expected exclusion contour for vector mediator")
            table.location = "Data from Figure 10 (left)"
            table.description = "Expected exclusion contour at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a vector mediator."
        
        table.add_variable(mmed)
        table.add_variable(mdm)
        table.add_image("./figures/PDF//Figure_010-a.pdf","./submission/")

    else:

        if type_plot == "observed" :
            table = Table("Observed exclusion contour for axial-vector mediator")
            table.location = "Data from Figure 10 (right)"
            table.description = "Observed exclusion contour at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming an axial-vector mediator."
        elif type_plot == "expected" :
            table = Table("Expected exclusion contour for axial-vector mediator")
            table.location = "Data from Figure 10 (right)"
            table.description = "Expected exclusion contour at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming an axial-vector mediator."
            
        table.add_variable(mmed)
        table.add_variable(mdm)
        table.add_image("./figures/PDF//Figure_010-b.pdf","./submission/")
        
    return table

################
def make_table_figure11(outdir,isPS=False,type_plot="observed"):

    if isPS == False:
        file = r.TFile("figures//Figure_011-a.root","READ")

        hobs = file.Get("Observed_limit");
        hexp = file.Get("Expected_limit");
        hexp_1s = file.Get("Expected_1sigma_band");
        hexp_2s = file.Get("Expected_2sigma_band");

        points_exp_1s = convert_graph_1d(hexp_1s);
        points_exp_2s = convert_graph_1d(hexp_2s);

        points_obs = {};
        points_exp = {};
        for key in ["x","y"]:
            points_obs[key] = [];
            points_exp[key] = [];

        for ipoint in range(0,hexp_1s.GetN()):
            x = r.Double();
            y = r.Double();
            hexp_1s.GetPoint(ipoint,x,y);
            points_obs["x"].append(x);
            points_obs["y"].append(hobs.Eval(x));
            points_exp["x"].append(x);
            points_exp["y"].append(hexp.Eval(x));

        mmed = Variable("Mediator mass", is_independent=True, is_binned=False, units="GeV")
        mmed.values = points_obs["x"];
        
        obs = Variable("Observed limit", is_independent=False, is_binned=False, units="")
        obs.values = points_obs["y"]
        
        exp_1s = Variable("Expected limit $\pm$ 1 s.d.", is_independent=False, is_binned=False, units="")
        exp_1s.values = points_exp["y"]

        exp_2s = Variable("Expected limit $\pm$ 2 s.d.", is_independent=False, is_binned=False, units="")
        exp_2s.values = points_exp["y"]
        
        exp_unc_1s = Uncertainty("1 s.d.",is_symmetric=False);
        exp_unc_2s = Uncertainty("2 s.d.",is_symmetric=False);
        error_1s = [];
        error_2s = [];
        for i in range(0,len(points_exp_1s["y"])):
            error_1s.append((points_exp_1s["y_error_dw"][i],points_exp_1s["y_error_up"][i]));
        for i in range(0,len(points_exp_2s["y"])):
            error_2s.append((points_exp_2s["y_error_dw"][i],points_exp_2s["y_error_up"][i]));

        exp_unc_1s.values = error_1s;
        exp_unc_2s.values = error_2s;
        exp_1s.uncertainties.append(exp_unc_1s);
        exp_2s.uncertainties.append(exp_unc_2s);

        table = Table("Exclusion upper limits for DM scalar mediators")
        table.location = "Data from Figure 11 (left)"
        table.description = "Exclusion limits at 95% CL on $\mu = \sigma/\sigma_{th}$ vs $m_{med}$ for $m_{DM} = 1$ GeV assuming a scalar mediator."
        
        table.add_variable(mmed)
        table.add_variable(obs)
        table.add_variable(exp_1s)
        table.add_variable(exp_2s)
        table.add_image("./figures/PDF//Figure_011-a.pdf","./submission/")

        return table

    else:

        file = r.TFile("figures//Figure_011-b.root","READ")

        if type_plot == "observed" :
            hobs = file.Get("Observed_limit");
        elif type_plot =="expected" :
            hobs = file.Get("Expected_limit");
        
        points_obs = convert_hist_2d(hobs,True);
    
        mmed = Variable("Mediator mass", is_independent=True, is_binned=True, units="GeV")
        mmed.values = points_obs["x_edges"];
        mdm = Variable("Dark matter mass", is_independent=True, is_binned=True, units="GeV")
        mdm.values = points_obs["y_edges"];
        
        if type_plot =="observed" :
            obs = Variable("Observed limit", is_independent=False, is_binned=False, units="")
        elif type_plot =="expected" :
            obs = Variable("Expected limit", is_independent=False, is_binned=False, units="")
        obs.values = points_obs["z"]
     
        if type_plot =="observed" :
            table = Table("Observed limits vs mMED-mDM for a pseudoscalar mediator")
            table.location = "Data from Figure 11 (right)"
            table.description = "Observed exclusion limits at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a pseudoscalar mediator. N.B.: the granularity in the presented limit values is reduced with respect to the published result."
        elif type_plot =="expected" :
            table = Table("Expected limits vs mMED-mDM for a pseudoscalar mediator")
            table.location = "Data from Figure 11 (right)"
            table.description = "Expected exclusion limits at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a pseudoscalar mediator. N.B.: the granularity in the presented limit values is reduced with respect to the published result."

        table.add_variable(mmed)
        table.add_variable(mdm)
        table.add_variable(obs)
        table.add_image("./figures/PDF//Figure_011-b.pdf","./submission/")

        return table

################
def make_table_figure11_cont(outdir,type_plot="observed"):

    file = r.TFile("figures//Figure_011-b_contour.root","READ")

    if type_plot == "observed" :
        hobs = file.Get("Observed_contour");
    elif type_plot == "expected" :    
        hobs = file.Get("Expected_contour");

    points_obs = convert_graph_1d(hobs);

    mmed = Variable("Mediator mass", is_independent=True, is_binned=False, units="GeV")
    mmed.values = points_obs["x"];

    mdm = Variable("Dark matter mass", is_independent=False, is_binned=False, units="GeV")
    mdm.values = points_obs["y"];

    if type_plot == "observed" :
        table = Table("Observed exclusion contour for pseudoscalar mediator")
        table.location = "Data from Figure 11 (right)"
        table.description = "Observed exclusion contour at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a pseudoscalar mediator."
    elif type_plot == "expected" :
        table = Table("Expected exclusion contour for pseudoscalar mediator")
        table.location = "Data from Figure 11 (right)"
        table.description = "Expected exclusion contour at 95% CL on $\mu = \sigma/\sigma_{th}$ in the $m_{med}-m_{DM}$ plane assuming a pseudoscalar mediator."
        
    table.add_variable(mmed)
    table.add_variable(mdm)
    table.add_image("./figures/PDF//Figure_011-b.pdf","./submission/")

    return table

#################
def make_table_figure12(outdir,isAV=False,type_plot="observed"):

    if isAV == False:
        file = r.TFile("figures//Figure_012-a.root","READ")
    else:
        file = r.TFile("figures//Figure_012-b.root","READ")

    if type_plot == "observed":
        hobs = file.Get("Observed_contour");
    elif type_plot == "expected":
        hobs = file.Get("Expected_contour");

    points_obs = convert_graph_1d(hobs);
    
    mmed = Variable("Mediator mass", is_independent=True, is_binned=False, units="GeV")
    mmed.values = points_obs["x"];
    gq = Variable("Coupling gq", is_independent=False, is_binned=False, units="")
    gq.values = points_obs["y"];
    
    if isAV == False:
        if type_plot == "observed":
            table = Table("Observed exclusion contour vs gq-mMED for a vector mediator")
            table.location = "Data from Figure 12 (left)"
            table.description = "Observed exclusion contour at 95% CL in the $g_{q}-m_{med}$ plane assuming a vector mediator."
        elif type_plot == "expected":
            table = Table("Expected exclusion contour vs gq-mMED for a vector mediator")
            table.location = "Data from Figure 12 (left)"
            table.description = "Exclusion exclusion contour at 95% CL in the $g_{q}-m_{med}$ plane assuming a vector mediator."
    else:
        if type_plot == "observed":
            table = Table("Observed exclusion contour vs gq-mMED for an AV mediator")
            table.location = "Data from Figure 12 (right)"
            table.description = "Observed exclusion contour at 95% CL in the $g_{q}-m_{med}$ plane assuming an axial-vector mediator."
        elif type_plot == "expected":
            table = Table("Expected exclusion contour vs gq-mMED for an AV mediator")
            table.location = "Data from Figure 12 (right)"
            table.description = "Exclusion exclusion contour at 95% CL in the $g_{q}-m_{med}$ plane assuming an axial-vector mediator."

    table.add_variable(mmed)
    table.add_variable(gq)

    if isAV == False:
        table.add_image("./figures/PDF//Figure_012-a.pdf","./submission/")
    else:
        table.add_image("./figures/PDF//Figure_012-b.pdf","./submission/")

    return table


##############
def make_table_figure13(outdir,isAV=False,type_plot="observed"):

    if isAV == False:
        file = r.TFile("figures//Figure_013-a.root","READ")
    else:
        file = r.TFile("figures//Figure_013-b.root","READ")

    if type_plot == "observed":
        hobs = file.Get("observed_dd");
    elif type_plot == "expected":
        hobs = file.Get("observed_dd");

    points_obs = convert_graph_1d(hobs);

    mdm = Variable("Dark matter mass", is_independent=True, is_binned=False, units="GeV")
    mdm.values = points_obs["x"];
        
    if not isAV:
        obs = Variable("Spin-independent DM-nucleon cross section ", is_independent=False, is_binned=False, units="cm$^{2}$")
    else:
        obs = Variable("Spin-dependent DM-nucleon cross section ", is_independent=False, is_binned=False, units="cms$^{2}$")

    obs.values = points_obs["y"]

    if not isAV:
        if type_plot == "observed":
            table = Table("Observed limits on the S.I. DM-nucleon scattering xsec")
            table.location = "Data from Figure 13 (left)"
            table.description = "Observed upper limits at 90% CL on $\sigma_{DM-nucleon}$ vs $m_{med}$ for vector mediator"
        elif type_plot == "expected":
            table = Table("Expected limits on the S.I. DM-nucleon scattering xsec")
            table.location = "Data from Figure 13 (left)"
            table.description = "Expected upper limits at 90% CL on $\sigma_{DM-nucleon}$ vs $m_{med}$ for vector mediator"
        
        table.add_variable(mdm)
        table.add_variable(obs)
        table.add_image("./figures/PDF//Figure_013-a.pdf","./submission/")

    else:

        if type_plot == "observed":
            table = Table("Observed limits on the S.D. DM-nucleon scattering xsec")
            table.location = "Data from Figure 13 (right)"
            table.description = "Observed upper limits at 90% CL on $\sigma_{DM-nucleon}$ vs $m_{med}$ for axial-vector mediator"
        elif type_plot == "expected":
            table = Table("Expected limits on the S.D. DM-nucleon scattering xsec")
            table.location = "Data from Figure 13 (right)"
            table.description = "Expected upper limits at 90% CL on $\sigma_{DM-nucleon}$ vs $m_{med}$ for axial-vector mediator"
        
        table.add_variable(mdm)
        table.add_variable(obs)
        table.add_image("./figures/PDF//Figure_013-b.pdf","./submission/")
        
    return table

def make_table_figure14(outdir,type_plot="observed"):

    file = r.TFile("figures//Figure_014.root","READ")
    
    if type_plot == "observed":
        hobs = file.Get("observed_dd");
    elif type_plot == "expected":
        hobs = file.Get("expected_dd");

    points_obs = convert_graph_1d(hobs);

    mdm = Variable("Dark matter mass", is_independent=True, is_binned=False, units="GeV")
    mdm.values = points_obs["x"];
        
    obs = Variable("Velocity averaged DM annihilation cross section", is_independent=False, is_binned=False, units="cm$^{2}$")
    obs.values = points_obs["y"]

    
    if type_plot == "observed":
        table = Table("Observed limits on the DM annihilation xsec")
        table.location = "Data from Figure 14"
        table.description = "Observed upper limits at 90% CL on velocity averaged DM annihilation cross section derived from those placed for pseudoscalar mediators."
    elif type_plot == "expected":
        table = Table("Expected limits on the DM annihilation xsec")
        table.location = "Data from Figure 14"
        table.description = "Expected upper limits at 90% CL on velocity averaged DM annihilation cross section derived from those placed for pseudoscalar mediators."
        
    table.add_variable(mdm)
    table.add_variable(obs)
    table.add_image("./figures/PDF//Figure_014.pdf","./submission/")
        
    return table

def make_table_figure20_and_21(outdir,isMonoV=False):

    if isMonoV == False:
        file = r.TFile("figures//Figure_020.root","READ")
    else:
        file = r.TFile("figures//Figure_021.root","READ")

    if isMonoV == False:
        hcorr = file.Get("correlation_monojet");
    else:
        hcorr = file.Get("correlation_monov");

    points_obs = convert_corr_2d(hcorr);

    obs_x = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    obs_x.values = points_obs["x_edges"];
    obs_y = Variable("Missing transverse momentum", is_independent=True, is_binned=True, units="GeV")
    obs_y.values = points_obs["y_edges"];
    
    corr_val = Variable("Correlation", is_independent=False, is_binned=False, units="")
    corr_val.values = points_obs["z"];
 
    if isMonoV == False:
        table = Table("Correlation between bins for the monojet SR")
        table.location = "Data from Figure 20"
        table.description = "Correlations between the predicted background yields in all the $p_{T}^{miss}$ bins of the monojet signal region."
    else:
        table = Table("Correlation between bins for the mono-V SR")
        table.location = "Data from Figure 21"
        table.description = "Correlations between the predicted background yields in all the $p_{T}^{miss}$ bins of the mono-V signal region."

    table.add_variable(obs_x)
    table.add_variable(obs_y)
    table.add_variable(corr_val)

    if isMonoV == False:
        table.add_image("./figures/PDF//Figure_020.pdf","./submission/")
    else:
        table.add_image("./figures/PDF//Figure_021.pdf","./submission/")

    return table



def main():

    # Write to this directory
    outdir = "./submission/"
    os.system("rm -rf "+outdir+"/*")
    submission = Submission()
    print "##### Generate Figure 5 #####"
    submission.add_table(make_table_figure5(outdir,False))
    submission.add_table(make_table_figure5(outdir,True))
    print "##### Generate Figure 6 top row #####"
    submission.add_table(make_table_figure6_top(outdir,False))
    submission.add_table(make_table_figure6_top(outdir,True))
    print "##### Generate Figure 6 bottom row #####"
    submission.add_table(make_table_figure6_bottom(outdir,False))
    submission.add_table(make_table_figure6_bottom(outdir,True))
    print "##### Generate Figure 7 top row #####"
    submission.add_table(make_table_figure7_top(outdir,False))
    submission.add_table(make_table_figure7_top(outdir,True))
    print "##### Generate Figure 7 bottom row #####"
    submission.add_table(make_table_figure7_bottom(outdir,False))
    submission.add_table(make_table_figure7_bottom(outdir,True))
    print "##### Generate Figure 8 #####"
    submission.add_table(make_table_figure8_and_9(outdir,False,True))
    submission.add_table(make_table_figure8_and_9(outdir,True,True))
    print "##### Generate Figure 9 #####"
    submission.add_table(make_table_figure8_and_9(outdir,False,False))
    submission.add_table(make_table_figure8_and_9(outdir,True,False))
    print "##### Generate Figure 10 --> limits #####"
    submission.add_table(make_table_figure10(outdir,False,"observed"))
    submission.add_table(make_table_figure10(outdir,False,"expected"))
    submission.add_table(make_table_figure10(outdir,True,"observed"))
    submission.add_table(make_table_figure10(outdir,True,"expected"))
    print "##### Generate Figure 10 --> contours #####"
    submission.add_table(make_table_figure10_cont(outdir,False,"observed"))
    submission.add_table(make_table_figure10_cont(outdir,False,"expected"))
    submission.add_table(make_table_figure10_cont(outdir,True,"observed"))
    submission.add_table(make_table_figure10_cont(outdir,True,"expected"))
    print "##### Generate Figure 11 --> limits #####"
    submission.add_table(make_table_figure11(outdir,False))
    submission.add_table(make_table_figure11(outdir,True,"observed"))
    submission.add_table(make_table_figure11(outdir,True,"expected"))
    print "##### Generate Figure 11 --> contours #####"
    submission.add_table(make_table_figure11_cont(outdir,"observed"))
    submission.add_table(make_table_figure11_cont(outdir,"expected"))
    print "##### Generate Figure 12 #####"
    submission.add_table(make_table_figure12(outdir,False,"observed"))
    submission.add_table(make_table_figure12(outdir,False,"expected"))
    submission.add_table(make_table_figure12(outdir,True,"observed"))
    submission.add_table(make_table_figure12(outdir,True,"expected"))
    print "##### Generate Figure 13 --> contours #####"
    submission.add_table(make_table_figure13(outdir,False,"observed"))
    submission.add_table(make_table_figure13(outdir,False,"expected"))
    submission.add_table(make_table_figure13(outdir,True,"observed"))
    submission.add_table(make_table_figure13(outdir,True,"expected"))
    print "##### Generate Figure 14 #####"
    submission.add_table(make_table_figure14(outdir,"observed"))
    submission.add_table(make_table_figure14(outdir,"expected"))
    print "##### Generate Figure 20 and 21 #####"
    submission.add_table(make_table_figure20_and_21(outdir,False))
    submission.add_table(make_table_figure20_and_21(outdir,True))
                         
    submission.read_abstract("./abstract/abstract.txt")
    submission.create_files(outdir)

    os.system("sed -i 's/!!python\/str//g' "+outdir+"/submission.yaml");
    os.system("tar -czvf submission.tar.gz "+outdir+"/*");

if __name__ == '__main__':
    main()

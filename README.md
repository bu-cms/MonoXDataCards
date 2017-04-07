# MonoXDataCards

Workspace with the background model are committed into the Workspace directory:
       workspace_MJ.root --> background model for monojet category
       workspace_MV.root --> background model for monoV category


Text datacards for each DM model are committed in the Datacard directory, with the following nomenclature:
     MV stands for monoV-category, MJ for monojet category
     There is one datacard for each control region i.e. ZM = Zmumu, ZE = Zee, GJ = gamma+jets, WM = Wmunu, WE = Wenu, SR = signal region.
     datacard_COMB_monoJ.txt is the combined monojet category card 
     datacard_COMB_monoV.text is the combined mono-V category card
     datacard_COMB.txt is the combined monojet + mono-V card


To run the asymptotic limit for alternative models, one needs to use datacard_COMB.txt sobstituting the existing model with the new one. Please rememeber to remove the lines which implements the statistical bin-by-bin uncertainties on the signal template i.e. MonoJ_SR_MJ_CMS_bin1_stat.
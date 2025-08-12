import ROOT
import os
import numpy as np
import tabulate
import array
print(os.getcwd())

from .hist_tools import get_max

def PidifScale(cutoff_time=0.1):
    """
    """
    # scale factors due to cutoffs
    kTau_pi = 26 # ns
    fPidifTime = cutoff_time
    return (1 - np.exp(-fPidifTime / kTau_pi))

def MudifScale(cutoff_time=0.05):
    """
    """
    # scale factors due to cutoffs
    tau_muon = 2197. #ns
    # f_50ps = 1 - np.exp(-0.05 / tau_muon)
    f_50ps = 1 - np.exp(-cutoff_time / tau_muon)
    return f_50ps

def get_time_corrections(T0, T1, verbose=False):
    tau_pi = 26.033 #(±0.005)
    tau_mu = 2196.9811 #(±0.0022)
    lambda_pi = 1. / tau_pi
    lambda_mu = 1. / tau_mu

    # time window correction for pi-mu-e
    frac_1 = lambda_pi*lambda_mu/(lambda_pi - lambda_mu)
    frac_2 = (ROOT.exp(-lambda_mu*T0) - ROOT.exp(-lambda_mu*T1))/lambda_mu 
    frac_3 = (ROOT.exp(-lambda_pi*T0) - ROOT.exp(-lambda_pi*T1))/lambda_pi
    frac_pimue = frac_1 * (frac_2 - frac_3)

    frac_pie = ROOT.exp(-lambda_pi*T0) - ROOT.exp(-lambda_pi*T1)
    R_time =  frac_pimue / frac_pie
    if verbose:
        print('pi-mu-e time window correction:', frac_pimue)
        print('pi-e time window correction:', frac_pie)
        print ('R_time = ', R_time)

    return frac_pie, frac_pimue, R_time


def fetch_histograms(rfile, hist_info, tail_fraction=False, scaling=False):
    """
    """
    # get the relevant histograms
    _hists_2d = {}

    for _name in hist_info.keys():
        if tail_fraction:
            h = rfile.Get(hist_info[_name]['tf_name'])
        else:
            h = rfile.Get(hist_info[_name]['name'])
        if scaling:
            h.Scale(hist_info[_name]['scale'])
        h.Scale(1. / hist_info[_name]['n_gen'])
        _hists_2d[_name] = h
    return _hists_2d

def scaled_histograms(hist_dict, scaling=1):
    _hist_2d_pioneer_run1 = {}
    for _name in hist_dict.keys():
        h = hist_dict[_name].Clone()
        h.SetName(h.GetName() + '_scaled_by_{}'.format(scaling))
        h.Scale(scaling)
        _hist_2d_pioneer_run1[_name] = h
    return _hist_2d_pioneer_run1

def get_time_hists(hist_dict, hist_info, energy_binning, rebin=None):
    """
    """
    hists_he = {}
    hists_le = {}

    le_bin_edges = energy_binning[0], energy_binning[1]
    he_bin_edges = energy_binning[1], energy_binning[2]
    
    for _name in hist_dict.keys():
        h_le = hist_dict[_name].ProjectionX(
            _name +'_proj_le', le_bin_edges[0], le_bin_edges[1])
        h_he = hist_dict[_name].ProjectionX(
            _name +'_proj_he', he_bin_edges[0], he_bin_edges[1])
        h_he.SetTitle(hist_info[_name]['title'])
        h_le.SetTitle(hist_info[_name]['title'])
        h_he.GetXaxis().SetTitle('Time [ns]')
        h_le.GetXaxis().SetTitle('Time [ns]')
        h_he.GetYaxis().SetTitle('Number of Events')
        h_le.GetYaxis().SetTitle('Number of Events')
        if rebin != None:
            h_le.Rebin(rebin)
            h_he.Rebin(rebin)
        hists_le[_name] = h_le
        hists_he[_name] = h_he

    return hists_le, hists_he

def get_singlebin_hists(hist_dict, hist_info, energy_binning, time_window=(5, 55,), rebin=None):
    """
    """
    hists_he = {}
    hists_le = {}

    le_bin_edges = energy_binning[0], energy_binning[1]
    he_bin_edges = energy_binning[1], energy_binning[2]
    for _name in hist_dict.keys():
        _hist_name = hist_dict[_name].GetName()
        h_le = hist_dict[_name].ProjectionX(
            _hist_name +'_proj_le_sb', le_bin_edges[0], le_bin_edges[1])
        h_he = hist_dict[_name].ProjectionX(
            _hist_name +'_proj_he_sb', he_bin_edges[0], he_bin_edges[1])

        _low_time_bin_edge  = h_le.GetXaxis().FindBin(time_window[0])
        _high_time_bin_edge = h_le.GetXaxis().FindBin(time_window[1])
        _err_le = array.array('d', [0])
        _le_int = h_le.IntegralAndError(_low_time_bin_edge, _high_time_bin_edge, _err_le)

        _low_time_bin_edge  = h_he.GetXaxis().FindBin(time_window[0])
        _high_time_bin_edge = h_he.GetXaxis().FindBin(time_window[1])
        _err_he = array.array('d', [0])
        _he_int = h_he.IntegralAndError(_low_time_bin_edge, _high_time_bin_edge, _err_he)
        h_le_1bin = ROOT.TH1F(
            h_le.GetName() + '_1bin',
            hist_info[_name]['title'],
            1,
            time_window[0],
            time_window[1])
        h_le_1bin.SetBinContent(1, _le_int)
        h_le_1bin.SetBinError(1, _err_le[0])

        h_he_1bin = ROOT.TH1F(
            h_he.GetName() + '_1bin',
            hist_info[_name]['title'],
            1,
            time_window[0],
            time_window[1])
        h_he_1bin.SetBinContent(1, _he_int)
        h_he_1bin.SetBinError(1, _err_he[0])
        
        h_he_1bin.GetXaxis().SetTitle('Time [ns]')
        h_le_1bin.GetXaxis().SetTitle('Time [ns]')
        h_he_1bin.GetYaxis().SetTitle('Number of Events')
        h_le_1bin.GetYaxis().SetTitle('Number of Events')
        hists_le[_name] = h_le_1bin
        hists_he[_name] = h_he_1bin

    return hists_le, hists_he

def get_energy_hists(hist_dict, hist_info, time_window, rebin=None):
    """
    """
    hists = {}

    if not isinstance(time_window, (list, tuple)):
        raise TypeError
    if len(time_window) != 2:
        raise ValueError('time_window has a wrong length!')
    
    t0, tf = time_window[0], time_window[1]
    _low_edge_bin = hist_dict[list(hist_dict.keys())[0]].GetYaxis().FindBin(t0)
    _high_edge_bin = hist_dict[list(hist_dict.keys())[0]].GetYaxis().FindBin(tf)

    for _name in hist_dict.keys():
        h_le = hist_dict[_name].ProjectionY(
            _name +'_proj_time', _low_edge_bin, _high_edge_bin)
        h_le.SetTitle(hist_info[_name]['title'])
        h_le.GetXaxis().SetTitle('Energy [MeV]')
        h_le.GetYaxis().SetTitle('Number of Events')
        if rebin != None:
            h_le.Rebin(rebin)
        hists[_name] = h_le

    return hists

def hist_info_dict(hist_names, hist_tf_names):
    hist_info = {
        'pienu': {
            'title': '#pi - e#nu',
            'scale': 1.23270e-4,#pienu exp
            'n_gen': 250 * 1e6,
            'tf_eff': 0.1,
            'tf_uncert': None,
            'color': ROOT.kRed + 2,
        },
        'michel': {
            'title': '#pi [DAR] - #mu [DAR] - e',
            'scale': 1.,
            'n_gen': 1000* 1e6,
            'tf_eff': 1e-7,
            'tf_uncert': 0.1,
            'color': ROOT.kBlue,
        },
        
        'pileup': {
            'title': 'pile up',
            'scale': 1,
            'n_gen': 1000 * 1e6,
            'tf_eff': 1e-5,
            'tf_uncert': 0.1,
            'color': ROOT.kGreen + 2,
        },
        'mudif': {
            'title': '#pi [DAR] - #mu [DIF] - e',
            'scale': MudifScale(),
            'n_gen': 100 * 1e6,
            'tf_eff': 1e-2,
            'tf_uncert': 0.5,
            'color': ROOT.kViolet,
        },
        'pidif': {
            'title': '#pi [DIF] - #mu [DAR] - e',
            'scale': PidifScale(),
            'n_gen': 100 * 1e6,
            'tf_eff': 1e-4,
            'tf_uncert': 0.5,
            'color': ROOT.kOrange - 3
        },
        'rbm': {
            'title': 'rbm',
            'n_gen': 1000 * 1e6,
            'scale': 1,
            'tf_eff': 1e-3,
            'tf_uncert': 0.5,
            'color': ROOT.kCyan,
            # 'color': ROOT.kAzure - 10,
        },
    }

    _hist_info_dict = {}
    for k in hist_names.keys():
        _hist_info_dict[k] = hist_info[k]
        _hist_info_dict[k]['name'] = hist_names[k]
        _hist_info_dict[k]['tf_name'] = hist_tf_names[k]

    return _hist_info_dict

import ROOT
import argparse
import array
from fit.helper import fetch_histograms, scaled_histograms, get_time_hists, get_energy_hists, hist_info_dict, get_singlebin_hists
from fit.plotter import energy_plot, time_plot, tail_fraction_2bin_plot
from fit.tabulater import tail_fraction_report, acceptance_report, yield_table_report
from fit.hist_tools import sum_hist, reset_errors, convert_to_counting

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    args = parser.parse_args()
    
    # low energy bin prescale
    le_prescale = 1000.
    # time window TFA
    time_window = [5, 55]
    # scaling
    n_pienu_exp = 2e8

    # time histograms
    energy_binning = (0, 130, 160) # 58MeV to end of spectra
    rebin = None

    # time fit binning
    time_fit_binning = [-300] + [5 * i for i in range(21)] + [200, 300, 400, 500]


    hist_names = {
        'pienu': 'pienu_et_fid_pienu',
        'michel': 'all_et_fid_mudar_pidartar',
        'pileup': 'mix_et_fid_pileup',
        'mudif': 'mudif_et_fid_mudif',
        'pidif': 'pidif_et_fid_pidif',
        'rbm': 'all_et_fid_mudar_nptr',
        }

    hist_tf_names = {
        'pienu': 'pienu_et_dedz_pienu',
        'michel': 'all_et_dedz_mudar_pidartar',
        'pileup': 'mix_et_dedz_pileup',
        'mudif': 'mudif_et_dedz_mudif',
        'pidif': 'pidif_et_dedz_pidif',
        'rbm': 'all_et_dedz_mudar_nptr',
        }
 


    rfile = ROOT.TFile(args.input)
    hist_info = hist_info_dict(hist_names, hist_tf_names)


    
    # unbiased selection
    hist_dict_all = fetch_histograms(rfile, hist_info, tail_fraction=False)

    # print (20 * '==')
    # print ()
    # print ('UNSCALED HISTOGRAM YIELDS')
    # print ('-------------------------')
    # yield_table_report(hist_dict_all)


    # normalise the histograms based on PIONEER Run1 expectations
    _scaling = n_pienu_exp / hist_dict_all['pienu'].Integral()
    scaled_hist_dict_all = scaled_histograms(hist_dict_all, scaling=_scaling)

    print ()
    print ('SCALED HISTOGRAM YIELDS')
    print ('-----------------------')
    yield_table_report(scaled_hist_dict_all, scale=_scaling)

    
    ene_hists_all = get_energy_hists(
        scaled_hist_dict_all,
        hist_info,
        time_window)

    c_all, p_all = energy_plot(ene_hists_all, hist_info, time_window, tail_fraction=False)
    c_all.Draw()
    c_all.Update()
    c_all.SaveAs('plots/energy_spectrum_all.pdf')

    # tail fraction analysis selection
    hist_dict_tf = fetch_histograms(rfile, hist_info, tail_fraction=True)
    scaled_hist_dict_tf = scaled_histograms(hist_dict_tf, scaling=_scaling)

    ene_hists_tf = get_energy_hists(
        scaled_hist_dict_tf,
        hist_info,
        time_window)

    c_tf, p_tf = energy_plot(ene_hists_tf, hist_info, time_window, tail_fraction=True)
    c_tf.Draw()
    c_tf.Update()
    c_tf.SaveAs('plots/energy_spectrum_tailfraction.pdf')

    print ('SCALED TAIL FRACTION HISTOGRAMS YIELDS')
    print ('---------------------------------------')
    yield_table_report(ene_hists_tf, scale=_scaling)
    
    tail_fraction_report(hist_dict_all['pienu'], hist_dict_tf['pienu'], energy_binning)
    # acceptance_report(hist_dict_all['pienu'], hist_dict_all['michel'])
    acceptance_report(scaled_hist_dict_all['pienu'], scaled_hist_dict_all['michel'])

    hists_le, hists_he = get_time_hists(
        scaled_hist_dict_all, hist_info,
        energy_binning,
        rebin=rebin)

    print (50 * '=')
    print ('LE bin yields')
    yield_table_report(hists_le, scale=_scaling / le_prescale)
    print (50 * '=')

    print (50 * '=')
    print ('HE bin yields')
    yield_table_report(hists_le, scale=_scaling)
    print (50 * '=')
    
    for h in  [v for _, v in hists_le.items()]:
        h.Scale(1 / le_prescale)


    for h in  [v for _, v in hists_le.items()] + [v for _, v in hists_he.items()]:
        # convert_to_counting(h)
        reset_errors(h, remove=True, poisson=False)

    data_he = sum_hist(hists_he, 'data_he')
    # convert_to_counting(data_he)
    reset_errors(data_he)

    data_le = sum_hist(hists_le, 'data_le')
    # convert_to_counting(data_le)
    reset_errors(data_le)

    c_le, p_le = time_plot(hists_le, data_le, hist_info, label='Low Energy Bin')
    c_le.Draw()
    c_le.Update()
    c_le.SaveAs('plots/time_spectrum_le.pdf')
    
    c_he, p_he = time_plot(hists_he, data_he, hist_info, label='High Energy Bin')
    c_he.Draw()
    c_he.Update()
    c_he.SaveAs('plots/time_spectrum_he.pdf')


    # Tail Fraction histograms
    hists_tf_le, hists_tf_he = get_singlebin_hists(
        scaled_hist_dict_tf, hist_info,
        energy_binning,
        time_window=time_window,
        rebin=rebin)
    print(50 * "=")
    print('Yield reports in the LE bin of the TFA')
    yield_table_report(hists_tf_le, scale=_scaling)
    print(50 * "=")

    print(50 * "=")
    print('Yield reports in the HE bin of the TFA')
    yield_table_report(hists_tf_he, scale=_scaling)
    print(50 * "=")

    # remove pileup  from TFA Low Energy Bin
    hists_tf_le.pop('pileup')

    # remove rbm and pidif from TFA high energy bin
    hists_tf_he.pop('rbm')
    hists_tf_he.pop('pidif')

    print(50 * "=")
    print('Yield reports in the LE bin of the TFA after removal of low stat samples')
    yield_table_report(hists_tf_le, scale=_scaling)
    print(50 * "=")
                       
    print(50 * "=")
    print('Yield reports in the HE bin of the TFA after removal of low stat samples')
    yield_table_report(hists_tf_he, scale=_scaling)
    print(50 * "=")

    data_tf_he = sum_hist(hists_tf_he, 'data_tf_he')
    data_tf_le = sum_hist(hists_tf_le, 'data_tf_le')

    c_tf_2bin, _ = tail_fraction_2bin_plot(hists_tf_he, hists_tf_le, data_tf_he, data_tf_le, hist_info)
    c_tf_2bin.Draw()
    c_tf_2bin.Update()
    c_tf_2bin.SaveAs('plots/twobins_tailfraction.pdf')


    for h in  [v for _, v in hists_tf_le.items()] + [v for _, v in hists_tf_he.items()]:
        # convert_to_counting(h)
        reset_errors(h, remove=True, poisson=False)
    # convert_to_counting(data_tf_he)
    reset_errors(data_tf_he)
    # convert_to_counting(data_tf_le)
    reset_errors(data_tf_le)

    
    binning_arr = array.array('d', time_fit_binning)
    for k, v in hists_le.items():
        v = v.Rebin(len(time_fit_binning) - 1, k + '_le', binning_arr)
        hists_le[k] = v
    for k, v in hists_he.items():
        v = v.Rebin(len(time_fit_binning) - 1, k + '_he', binning_arr)
        hists_he[k] = v
        
    data_le_rb = data_le.Rebin(len(time_fit_binning) - 1, data_le.GetName() + '_rebin', binning_arr)
    data_he_rb = data_he.Rebin(len(time_fit_binning) - 1, data_he.GetName() + '_rebin', binning_arr)

    fit_input = ROOT.TFile('output/fit_input.root', 'recreate')
    fit_input.cd()
    for k, v in hists_tf_he.items():
        v.SetName(k + '_tf_he')
        fit_input.Add(v)
    fit_input.Add(data_tf_he)
    
    for k, v in hists_tf_le.items():
        v.SetName(k + '_tf_le')
        fit_input.Add(v)
    fit_input.Add(data_tf_le)

    for k, v in hists_le.items():
        v.SetName(k + '_le')
        fit_input.Add(v)
    fit_input.Add(data_le_rb)

    for k, v in hists_he.items():
        v.SetName(k + '_he')
        fit_input.Add(v)
    fit_input.Add(data_he_rb)

    fit_input.Write()
    fit_input.Print()
    fit_input.Close()
    
    
        

        

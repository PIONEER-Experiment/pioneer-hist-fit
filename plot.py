import ROOT
import argparse
from fit.helper import fetch_histograms, scaled_histograms, get_time_hists, get_energy_hists, hist_info_dict, get_singlebin_hists
from fit.plotter import shape_comparison
from fit.tabulater import tail_fraction_report
from fit.hist_tools import sum_hist, reset_errors, convert_to_counting

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    args = parser.parse_args()
    
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



    # pienu vs michel
    hists_2d  = {
        'pienu': rfile.Get('pienu_edeps_fid_pienu'),
        'michel': rfile.Get('all_edeps_fid_mudar_pidartar')
        }

    c, _ = shape_comparison(hists_2d, hist_info, do_2d=True)
    c.Draw()
    c.Update()
    c.SaveAs('plots/edeps_pienu_vs_michel.png')
    
    hists_2d  = {
        'pienu': rfile.Get('pienu_voca_fid_pienu').ProjectionX(),
        'michel': rfile.Get('all_voca_fid_mudar_pidartar').ProjectionX()
        }

    c2, _ = shape_comparison(hists_2d, hist_info, log_y=True)
    c2.Draw()
    c2.Update()
    c2.SaveAs('plots/distance_of_closest_approach_pienu_vs_michel.png')

    # pienu vs mudif
    hists_2d  = {
        'pienu': rfile.Get('pienu_e5dedz_fid_pienu').ProjectionY(),
       'mudif': rfile.Get('mudif_e5dedz_fid_mudif').ProjectionY()
    }

    c3, _ = shape_comparison(hists_2d, hist_info, log_y=True)
    c3.Draw()
    c3.Update()
    c3.SaveAs('plots/dedz_pienu_vs_mudif.png')

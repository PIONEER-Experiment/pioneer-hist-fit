import ROOT

from fit.measurement import make_measurement
from fit.tabulater import fitresult_report
from fit.plotter import nuisance_params, fitresult_cat_projection, ranking_plot, correlation_plot

def fit(ws, verbose=True, n_tries=5, res_name='results'):
    """
    """
    mc = workspace.obj('ModelConfig')
    #workspace.loadSnapshot('snapshot_paramsVals_initial')
    data = workspace.obj('obsData')
    # data = workspace.obj('asimovData')
    pdf = mc.GetPdf()
    nps = mc.GetNuisanceParameters()
    fcn = pdf.createNLL(
        data, 
        ROOT.RooFit.Constrain(mc.GetNuisanceParameters()),
        ROOT.RooFit.GlobalObservables(mc.GetGlobalObservables()))

    if verbose:
        fcn.Print()
        printLevel = ROOT.Math.MinimizerOptions.DefaultPrintLevel()
    else:
        printLevel = -2
    ROOT.RooMsgService.instance().globalKillBelow()
    if( printLevel < 0 ):
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    strat = ROOT.Math.MinimizerOptions.DefaultStrategy()

    minim = ROOT.RooMinimizer(fcn)
    minim.setEps(1e-3/0.001)
    minim.setStrategy(1)
    # minim.setPrintLevel(printLevel)
    if verbose:
        minim.setPrintLevel(1)
        minim.setProfile()
    minim.setOffsetting(True)
    minim.optimizeConst(2)

    status = minim.minimize(
        ROOT.Math.MinimizerOptions.DefaultMinimizerType(),
        ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo())
    if status % 1000 == 0:
        print ('worked with default stratey!')
        print (status)
        _cov_status = minim.hesse()
        # print (_cov_status)
        # if _cov_status != 3:
        #     print ('hesse status != 3, running minos')
        _cov_status = minim.minos()
        res = minim.save(res_name)
        return res

    #print ('Try strategy 2')
    minim.setStrategy(2)
    status = minim.minimize(
            ROOT.Math.MinimizerOptions.DefaultMinimizerType(),
            ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo())
    if status % 1000 == 0:
        print ('worked with strategy 2!')
        _cov_status = minim.hesse()
        _cov_status = minim.minos()
        # if _cov_status != 3:
        #     print ('hesse status != 3, running minos')
        #     _cov_status = minim.minos()
        res = minim.save(res_name)
        return res

    #print ('trying with improve')
    status = minim.minimize('Minuit2', 'migradimproved')
    if status % 1000 == 0:
        print ('worked with migradimproved!')
        _cov_status = minim.hesse()
        # if _cov_status != 3:
        #     print ('hesse status != 3, running minos')
        #     _cov_status = minim.minos()
        res = minim.save(res_name)
        return res

    print ('return with final attempt (migradimproved)')
    _cov_status = minim.hesse()
    # if _cov_status != 3:
    #     print ('hesse status != 3, running minos')
    # _cov_status = minim.minos()
    res = minim.save(res_name)
    return res


def formulas(workspace):
    R_time =  0.23482407672465772
    R_unknown =  0.9990791438407314
    R_prescale = 1. / 1000.
    tf_corr_unknown = 1.0236239178895448
    
    pie_HE_prefit = workspace.obj('pienu_he_high_energy_shapes').createIntegral(
    workspace.obj('obs_x_high_energy')).getVal()

    michel_LE_prefit = workspace.obj('michel_le_low_energy_shapes').createIntegral(
        workspace.obj('obs_x_low_energy')).getVal()

    pie_tf_he_prefit = workspace.obj('pienu_tf_he_high_energy_tail_fraction_shapes').createIntegral(
        workspace.obj('obs_x_high_energy')).getVal()

    pie_tf_le_prefit = workspace.obj('pienu_tf_le_low_energy_tail_fraction_shapes').createIntegral(
        workspace.obj('obs_x_low_energy')).getVal()

    # print(pie_HE_prefit)
    # print(michel_LE_prefit)
    # print(pie_tf_he_prefit)
    # print(pie_tf_le_prefit)
    tf = pie_tf_le_prefit/(pie_tf_he_prefit+pie_tf_le_prefit)
    # print(tf)
    # print(pie_HE_prefit/michel_LE_prefit*R_time*R_unknown*1e4)
    # print(pie_HE_prefit/michel_LE_prefit*R_time*R_unknown*(1+tf)*1e4)
    
    r_pi_formula_notf = "(pie_HE * {pie_he_prefit}) / (pidar_mudar_LE * {pimue_le_prefit}) * {R_time} * {R_unknown} * {R_prescale}".format(
        pie_he_prefit=pie_HE_prefit,
        pimue_le_prefit=michel_LE_prefit,
        R_time=R_time,
        R_prescale=R_prescale,
        R_unknown=R_unknown)
    r_pi_notf = ROOT.RooFormulaVar(
        "R_pi_notf", 
        r_pi_formula_notf, 
        ROOT.RooArgList(workspace.obj("pie_HE"), workspace.obj("pidar_mudar_LE")))

    tf_formula = "pie_tf_LE * {pi_tf_le_prefit} / (pie_tf_HE * {pi_tf_he_prefit} + pie_tf_LE * {pi_tf_le_prefit})".format(
        pi_tf_le_prefit=pie_tf_le_prefit,
        pi_tf_he_prefit=pie_tf_he_prefit)
    tf = ROOT.RooFormulaVar(
        'tail_fraction',
        tf_formula,
        ROOT.RooArgList(workspace.obj('pie_tf_LE'), workspace.obj('pie_tf_HE')))

    tf_formula_corr = "pie_tf_LE * {pi_tf_le_prefit} / (pie_tf_HE * {pi_tf_he_prefit} + pie_tf_LE * {pi_tf_le_prefit}) * {tf_corr_unknown}".format(
        pi_tf_le_prefit=pie_tf_le_prefit,
        pi_tf_he_prefit=pie_tf_he_prefit,
        tf_corr_unknown=tf_corr_unknown)

    r_pi_formula = '{}*(1 + {})'.format(r_pi_formula_notf, tf_formula_corr)
    arg_list = ROOT.RooArgList()
    arg_list.add(workspace.obj("pie_HE"))
    arg_list.add(workspace.obj("pidar_mudar_LE"))
    arg_list.add(workspace.obj("pie_tf_LE"))
    arg_list.add(workspace.obj("pie_tf_HE"))
    r_pi = ROOT.RooFormulaVar(
        "R_pi", 
        r_pi_formula,
        arg_list)

    return r_pi, r_pi_notf, tf

def compute_np_impact(workspace, mc, r_pi):

    _res = {}
    for var in list(mc.GetNuisanceParameters()) + list(mc.GetParametersOfInterest()):
        print()
        print (20* '-')
        print (var.GetName(), var.getVal(), var.getErrorLo(), var.getErrorHi())
        workspace.loadSnapshot('snapshot_paramsVals_final')
        var.setVal(var.getVal())
        var.setConstant(True)
        res = fit(workspace)
        _nom = r_pi.getVal() 
        workspace.loadSnapshot('snapshot_paramsVals_final')
        var.setVal(var.getVal() + var.getErrorHi())
        var.setConstant(True)
        res = fit(workspace)
        _up = r_pi.getVal()
        workspace.loadSnapshot('snapshot_paramsVals_final')
        var.setVal(var.getVal() - var.getErrorLo())
        var.setConstant(True)
        res = fit(workspace)
        _do = r_pi.getVal()
        
        workspace.loadSnapshot('snapshot_paramsVals_final')
        _res[var.GetName()] = {
            'rpi_down_var': (_nom - _do)/_nom * 1e2,
            'rpi_up_var': (_up - _nom) / _nom *1e2,
            'nuis_nom': var.getVal(),
            'nuis_down': var.getErrorLo(),
            'nuis_up': var.getErrorHi(),
        }
    return _res

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--make-meas', default=False, action='store_true')
    parser.add_argument('--run-ranking', default=False, action='store_true')
    parser.add_argument('--compute-rpi', default=False, action='store_true')
    parser.add_argument('--fit-projections', default=False, action='store_true')
    args = parser.parse_args()

    if args.make_meas:
        make_measurement()
        # make_measurement(
        #     do_he=False,
        #     do_le=False,
        #     do_tf_he=False,
        #     do_tf_le=True)
    
    rfile = ROOT.TFile('output/measurement.root')
    workspace = rfile.Get('combined')
    mc  = workspace.obj('ModelConfig')
    workspace.saveSnapshot('snapshot_paramsVals_initial', mc.GetNuisanceParameters())
    params = [
        # 'alpha_mudif_LE',
        # 'alpha_pie_LE',
        # 'alpha_beam_muons_LE',
        # 'pidif_LE',
        # 'alpha_pileup_LE',
        # 'alpha_pileup_HE',
        # 'pidar_mudar_HE',
        'Lumi',
    ]
    for p in params:
        workspace.obj(p).setConstant(True)
    
    res = fit(workspace, verbose=False)
    res.Print()
    workspace.saveSnapshot('snapshot_paramsVals_final', mc.GetNuisanceParameters())


    c_np, plotables = nuisance_params(workspace, mc)
    c_np.Draw()
    c_np.Update()
    c_np.SaveAs('plots/pulls.pdf')

    if args.compute_rpi:
        r_pi, r_pi_notf, tf = formulas(workspace)
        fitresult_report(r_pi, r_pi_notf, tf, res)

    c_corr, _ = correlation_plot(res)
    c_corr.Draw()
    c_corr.Update()
    c_corr.SaveAs('plots/correlations.pdf')

    if args.fit_projections:
        fitresult_cat_projection(workspace, mc)

    if args.run_ranking:
        if args.compute_rpi:
            _res = compute_np_impact(workspace,  mc, r_pi)
            c1, _ = ranking_plot( _res)
        else:
            _res = compute_np_impact(workspace,  mc, mc.GetParametersOfInterest()[0])
            c1, _ = ranking_plot( _res, second_axis_title="#DeltaN_{#pi-e}/N_{#pi-e} [%]", second_axis_color=ROOT.kRed, second_axis_scale=160)

        c1.Draw()
        c1.SaveAs('plots/ranking.pdf')


import ROOT
import array
from .hist_tools import get_max
ROOT.gROOT.SetStyle('ATLAS')
ROOT.gROOT.SetBatch(True)

def energy_plot(ene_hists, hist_info, time_window, tail_fraction=False):
    c = ROOT.TCanvas()
    c.SetTopMargin(0.1)
    c.SetLogy()
    _max = get_max(ene_hists)
    for i_k, k in enumerate(ene_hists.keys()):
        ene_hists[k].GetYaxis().SetRangeUser(0.1, 100 * _max)
        ene_hists[k].SetLineColor(hist_info[k]['color'])
        if i_k == 0:
            ene_hists[k].Draw('HIST')
        else:
            ene_hists[k].Draw('SAMEHIST')

    leg = ROOT.TLegend(c.GetLeftMargin(), 1 - c.GetTopMargin(), 1 - c.GetRightMargin(), 1)
    leg.SetTextSize(0.3*leg.GetTextSize())
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    for k in ['pienu', 'mudif', 'pileup', 'michel', 'pidif', 'rbm']:
        if k in ene_hists.keys():
            leg.AddEntry(ene_hists[k], hist_info[k]['title'], 'l')
    leg.SetNColumns(3)
    leg.Draw()
            
    label = ROOT.TLatex()
    label.SetTextSize(0.9*label.GetTextSize())
    label.SetNDC(True)
    label.DrawLatex(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.06, '#bf{PIONEER} Phase I Simulation')
    if tail_fraction:
        label.DrawText(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.1, 'Tail Fraction Analysis Region')
    else:
        label.DrawText(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.1, 'Inclusive Region')
    label.DrawText(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.14, 'Time Window: [{}, {}] ns'.format(
        time_window[0], time_window[1]))
    
    _drawables = [leg, label]
    # stack = ROOT.THStack()
    # stack.Add(ene_hists['pileup'])
    # stack.Add(ene_hists['mudif'])
    # stack.Add(ene_hists['pienu'])
    # stack.Add(ene_hists['beam_muons'])
    # stack.Add(ene_hists['pidif'])
    # stack.Add(ene_hists['michel'])
    # stack.Draw('HIST')
    # stack.GetYaxis().SetTitle(ene_hists['pienu'].GetYaxis().GetTitle())
    # stack.GetXaxis().SetTitle(ene_hists['pienu'].GetXaxis().GetTitle())
    return c, _drawables
            


def time_plot(hists, data_hist, hist_info, label='High Energy Bin'):

    # high energy bin
    c = ROOT.TCanvas()
    c.SetTopMargin(0.15)
    c.cd()
    c.SetLogy()
    for k, h in hists.items():
        h.SetLineColor(hist_info[k]['color'])

    samples = ['pienu', 'mudif', 'pileup', 'pidif', 'rbm', 'michel']
    if 'High Energy Bin' in label:
        samples = ['rbm', 'pidif', 'mudif', 'michel', 'pileup', 'pienu']
    stack = ROOT.THStack()
    _max = get_max(hists)
    for s in samples:
        if s in hists.keys():
            stack.Add(hists[s])


    data_hist.SetMarkerColor(1)
    data_hist.SetMarkerStyle(20)
    data_hist.SetMarkerSize(0.5)
    stack.SetMinimum(1)
    stack.SetMaximum(100 * _max)
    stack.Draw("HIST")
    stack.GetXaxis().SetTitle(hists[list(hists.keys())[0]].GetXaxis().GetTitle())
    stack.GetYaxis().SetTitle(hists[list(hists.keys())[0]].GetYaxis().GetTitle())
    data_hist.Draw('samePE')

    leg = ROOT.TLegend(c.GetLeftMargin(), 1 - c.GetTopMargin(), 1 - c.GetRightMargin(), 1)
    leg.SetTextSize(0.3 * leg.GetTextSize())
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    for k in samples:
        if k in hists.keys():
            leg.AddEntry(hists[k], hist_info[k]['title'], 'l')
    leg.AddEntry(data_hist, 'Pseudo Data', 'p')
    leg.SetNColumns(4)
    leg.Draw()

    _ttext = ROOT.TLatex()
    _ttext.SetTextSize(0.9*_ttext.GetTextSize())
    _ttext.SetNDC(True)
    _ttext.DrawLatex(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.06, '#bf{PIONEER} Phase I Simulation')
    _ttext.DrawText(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.11, label)

    _drawables = [stack, leg, _ttext]
    return c, _drawables


def tail_fraction_2bin_plot(hists_he_tf, hists_le_tf, data_he_tf, data_le_tf, hist_info):
    c2 = ROOT.TCanvas()
    c2.Divide(2, 1)
    pad1 = c2.GetPad(1)
    pad1.SetTopMargin(0.1)
    pad1.cd()
    stack_2 = ROOT.THStack()


    for k, v in hists_he_tf.items():
        v.SetLineColor(hist_info[k]['color'])

    for k, v in hists_le_tf.items():
        v.SetLineColor(hist_info[k]['color'])

    samples = ['mudif', 'pileup', 'pidif', 'beam_muons', 'michel', 'pienu']
    for s in samples:
        if s in hists_le_tf.keys():
            stack_2.Add(hists_le_tf[s])


    data_le_tf.SetMarkerColor(1)
    data_le_tf.SetMarkerStyle(20)
    data_le_tf.SetMarkerSize(1)
    # data_le_tf.SetLineWidth(0)
    
    stack_2.SetMinimum(1e-1)
    stack_2.SetMaximum(data_le_tf.GetBinContent(data_le_tf.GetMaximumBin())*1.2)
    stack_2.Draw("HIST")
    stack_2.GetXaxis().SetTitle(hists_le_tf[list(hists_le_tf.keys())[0]].GetXaxis().GetTitle())
    stack_2.GetYaxis().SetTitle(hists_le_tf[list(hists_le_tf.keys())[0]].GetYaxis().GetTitle())
    data_le_tf.Draw('samePE')
    pad1.Update()
    _ttext = ROOT.TText()
    _ttext.SetNDC(True)
    _ttext.DrawText(
        0.3,  
        #    0.5,
        0.95,
        'Tail Fraction: Low Energy Bin')


    pad2 = c2.GetPad(2)
    pad2.SetTopMargin(0.1)
    pad2.cd()
    pad2.SetLogy()
    stack_3 = ROOT.THStack()
    samples = ['pidif', 'beam_muons', 'michel', 'pileup', 'mudif', 'pienu']
    for s in samples:
        if s in hists_he_tf.keys():
            stack_3.Add(hists_he_tf[s])
            

    data_he_tf.SetMarkerColor(1)
    data_he_tf.SetMarkerStyle(20)
    data_he_tf.SetMarkerSize(1)

    stack_3.SetMinimum(1e-1)
    stack_3.SetMaximum(data_he_tf.GetBinContent(data_he_tf.GetMaximumBin())*10)
    stack_3.Draw("HIST")
    stack_3.GetXaxis().SetTitle(hists_he_tf[list(hists_he_tf.keys())[0]].GetXaxis().GetTitle())
    stack_3.GetYaxis().SetTitle(hists_he_tf[list(hists_he_tf.keys())[0]].GetYaxis().GetTitle())
    data_he_tf.Draw('samePE')
    _ttext.DrawText(
        0.3, 
        #    0.5,
        0.95,
        'Tail Fraction: High Energy Bin')

    _drawables = [stack_2, stack_3, pad1, pad2, _ttext]
    return c2, _drawables

def shape_comparison(hists,  hist_info, log_y=False, do_2d=False):
    if do_2d:
        c_dummy = ROOT.TCanvas()
        _graphs = {}
        for k, v in hists.items():
            v.Draw('CONT LIST')
            c_dummy.Update()
            contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
            _graph_list = []
            for i in range(contours.GetSize()):
                _list = contours.At(i)
                for j in range(_list.GetSize()):
                    gr = _list.At(j).Clone()
                    gr.SetLineColor(ROOT.kRed)
                    _graph_list += [gr]
            _graphs[k] = _graph_list

    c = ROOT.TCanvas()
    if log_y:
        c.SetLogy()
    c.SetTopMargin(0.1)

    _max = get_max(hists)
    for k, v in hists.items():
        v.SetLineColor(hist_info[k]['color'])
        v.SetMarkerColor(hist_info[k]['color'])
        if isinstance(v, ROOT.TH1):
            v.GetYaxis().SetRangeUser(0, 1.1*_max)
            if log_y:
                v.GetYaxis().SetRangeUser(1e-6, 10*_max)
        v.Scale(1. / v.Integral())
        
    if isinstance(v, ROOT.TH1):
        hists['pienu'].Draw('HIST')
        for k , v in hists.items():
            if k == 'pienu':
                continue
            v.Draw('sameHIST')

    if do_2d:
        hist_template = hists['pienu'].Clone()
        hist_template.Reset()
        hist_template.Draw()
        hist_template.Print()
        for k, v in _graphs.items():
            for _gr in v:
                _gr.SetLineColor(hist_info[k]['color'])
                _gr.Draw('sameL')

    leg = ROOT.TLegend(c.GetLeftMargin(), 1 - c.GetTopMargin(), 1 - c.GetRightMargin(), 1)
    leg.SetTextSize(0.3*leg.GetTextSize())
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(2)
    for k in hists.keys():
        leg.AddEntry(hists[k], hist_info[k]['title'], 'lp')
    leg.Draw()
            
    label = ROOT.TLatex()
    label.SetTextSize(0.9*label.GetTextSize())
    label.SetNDC(True)
    label.DrawLatex(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.06, '#bf{PIONEER} Phase I Simulation')
    label.DrawText(c.GetLeftMargin() + 0.02, 1 - c.GetTopMargin() - 0.1, 'Inclusive Region')

    c.Update()
    _drawables = [leg, label]
    if do_2d:
        _drawables += [_graphs, hist_template]
    return c, _drawables


def nuisance_params(workspace, mc):
    """
    """
    np_vals = {}
    for p in list(mc.GetNuisanceParameters()) + list(mc.GetParametersOfInterest()):
        if p.GetName() == 'Lumi':
            continue
            # LumiRelError = 1e-6
            # pull    = (p.getVal() - workspace.var("nominalLumi").getVal()) / (workspace.var("nominalLumi").getVal())
            # errorHi = p.getErrorHi() / (workspace.var("nominalLumi").getVal() * LumiRelError)
            # errorLo = p.getErrorLo() / (workspace.var("nominalLumi").getVal() * LumiRelError)
            # np_vals[p.GetName()] = (errorLo, pull, errorHi)
        else:
            np_vals[p.GetName()] = {
                'lo': p.getErrorLo(),
                'nom': p.getVal(),
                'hi': p.getErrorHi(),
                'err': p.getError()
                }
        print (p.GetName(), np_vals[p.GetName()], p.getError())

            
    haxis= ROOT.TH1F('axis', 'axis', len(np_vals.keys()), 0, len(np_vals.keys()))
    np_graph = ROOT.TGraphAsymmErrors()
    for i_k, k in enumerate(sorted(np_vals.keys())):
        np_graph.SetPoint(i_k, i_k+0.5, np_vals[k]['nom'])
        np_graph.SetPointError(i_k, 0, 0, abs(np_vals[k]['err']), np_vals[k]['err']) 
        haxis.GetXaxis().SetBinLabel(i_k + 1, k)
    haxis.GetYaxis().SetRangeUser(-3, 3)
    haxis.GetYaxis().SetTitle('Pull [#sigma] or Norm Factor')
    haxis.LabelsOption("v")
    box_2sig = ROOT.TBox(0, -2, len(np_vals.keys()), 2)
    box_1sig = ROOT.TBox(0, -1, len(np_vals.keys()), 1)
    box_2sig.SetFillColor(ROOT.kYellow)
    box_1sig.SetFillColor(ROOT.kGreen)
    c = ROOT.TCanvas()
    c.SetBottomMargin(0.5)
    haxis.Draw()
    box_2sig.Draw('same')
    box_1sig.Draw('same')
    haxis.Draw('same')
    np_graph.Draw('samePE')
    c.RedrawAxis()
    _plotables = [haxis, np_graph, box_2sig, box_1sig]
    return c, _plotables


def fitresult_cat_projection(workspace, mc):
    pdf = mc.GetPdf()
    cat = pdf.indexCat()
    data = workspace.obj('obsData')

    iterator = cat.typeIterator()
    tt = iterator.Next()
    canvases = []
    while (tt):
        tt.Print()
        cat_pdf = pdf.getPdf(tt.GetName())
        obs = cat_pdf.getObservables(mc.GetObservables())
        cat_data = data.reduce('{}=={}::{}'.format(cat.GetName(), cat.GetName(), tt.GetName()))
        #      localData->reduce(Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName()));
        #  RooArgSet*  obstmp = pdftmp->getObservables(*mc->GetObservables());
        # RooRealVar* obs    = ((RooRealVar*)obstmp->first());
        
        c = ROOT.TCanvas()
        rooplot = obs[0].frame()
        rooplot.SetYTitle("Events");
        #cat_pdf.plotOn(
        #    rooplot, 
        #    ROOT.RooFit.FillColor(ROOT.kOrange), 
        #    ROOT.RooFit.LineWidth(2), 
        #    ROOT.RooFit.LineColor(ROOT.kBlue), 
        #    ROOT.RooFit.VisualizeError(res, 1),
        #    ROOT.RooFit.Normalization(postFitIntegral, ROOT.RooAbsReal.NumEvent), 
        #    ROOT.RooFit.Name("FitError"))
        workspace.loadSnapshot('snapshot_paramsVals_initial')
        preFitIntegral = cat_pdf.expectedEvents(obs)
        cat_pdf.plotOn(
        rooplot, 
            ROOT.RooFit.LineWidth(2), 
            ROOT.RooFit.Normalization(preFitIntegral, ROOT.RooAbsReal.NumEvent), 
            ROOT.RooFit.Name("Prefit"))
        workspace.loadSnapshot('snapshot_paramsVals_final')
        postFitIntegral = cat_pdf.expectedEvents(obs)
        cat_pdf.plotOn(
            rooplot, 
            ROOT.RooFit.LineWidth(2), 
            ROOT.RooFit.LineColor(ROOT.kRed),
            ROOT.RooFit.Normalization(postFitIntegral, ROOT.RooAbsReal.NumEvent), 
            ROOT.RooFit.Name("Postfit"))
        
        modelName1  = tt.GetName() + '_model'
        pdfmodel = (cat_pdf.getComponents()).find(modelName1)
        funcList1 = pdfmodel.funcList()
        for ic, comp in enumerate(funcList1):
            Ntemp = (comp.createIntegral(*obs)).getVal()
            #print(tt.GetName(), comp.GetName(), Ntemp)
            cat_pdf.plotOn(
                rooplot, 
                ROOT.RooFit.LineWidth(2), 
                ROOT.RooFit.LineColor(ic + 1),
                ROOT.RooFit.Components(comp),
                ROOT.RooFit.Normalization(Ntemp, ROOT.RooAbsReal.NumEvent), 
                ROOT.RooFit.Name(comp.GetName()))
            
        #cat_data.Print()
        cat_data.plotOn(
            rooplot,
            ROOT.RooFit.MarkerSize(1),
            #ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),
            #ROOT.RooFit.DataError(0),
            ROOT.RooFit.Name("Data"))
        
        rooplot.Draw()
        #c.SetLogy()
        c.Draw()
        c.SaveAs('plots/plot_{}.pdf'.format(tt.GetName()))
        canvases += [c]
        tt = iterator.Next()

def ranking_plot(
        res,
        second_axis_title="#DeltaR_{e/#mu}/R_{e/#mu} [%]",
        second_axis_color=ROOT.kBlue,
        second_axis_scale=80):
    """
    """
    _ranked_nps = sorted(res.keys(), key=lambda r: max(abs(res[r]['rpi_down_var']), abs(res[r]['rpi_up_var'])))
    _ranked_nps = list(filter(lambda np: 'alpha' in np, _ranked_nps))
    for _np in _ranked_nps:
        print(_np)

    c = ROOT.TCanvas()
    c.SetLeftMargin(0.4)
    c.SetTopMargin(0.2)
    _offset = 5
    h_template = ROOT.TH2F('h_temp','', 1, -2, 2, len(_ranked_nps) + _offset, 0, len(_ranked_nps) + _offset)
    h_template.GetYaxis().SetRangeUser(0, len(_ranked_nps) + _offset)
    gr_poi_impact = ROOT.TGraphAsymmErrors()
    gr_np_pull_and_constraint = ROOT.TGraphAsymmErrors()
    _scale = second_axis_scale
    # _scale = 160#N_pie
    # _scale = 4000 #N_pimue
    for i_np, _np in enumerate(_ranked_nps):
        h_template.GetYaxis().SetBinLabel(i_np + _offset, _np)
        gr_poi_impact.SetPoint(i_np, 0, i_np + _offset - 0.5)
        gr_poi_impact.SetPointError(
            i_np,
            abs(res[_np]['rpi_down_var']) * _scale,
            abs(res[_np]['rpi_up_var']) * _scale
            ,  0.45, 0.45,)
        if 'alpha' in _np:
            gr_np_pull_and_constraint.SetPoint(i_np, res[_np]['nuis_nom'], i_np + _offset - 0.5)
            gr_np_pull_and_constraint.SetPointError(
                i_np, 
                abs(res[_np]['nuis_down']),
                abs(res[_np]['nuis_up']), 0, 0)
        else:
            gr_np_pull_and_constraint.SetPoint(i_np, res[_np]['nuis_nom']-1, i_np + _offset - 0.5)
            gr_np_pull_and_constraint.SetPointError(
                i_np, 
                res[_np]['nuis_down'],
                res[_np]['nuis_up'], 0, 0)

    gr_poi_impact.SetFillColor(ROOT.kBlue)
    gr_poi_impact.SetLineColor(ROOT.kBlue)
    gr_poi_impact.SetFillStyle(3001)


    h_template.Draw()
    h_template.GetXaxis().SetTitle('Pull = (#hat{#theta} - #theta_{0})/#Delta#theta')
    gr_poi_impact.Draw('sameE2')
    gr_np_pull_and_constraint.Draw('samePE')
    ROOT.gStyle.SetPadTickX(0)

    c.Draw()
    # Create a second axis
    xmin = array.array('d', [0])
    xmax = array.array('d', [0])
    ymin = array.array('d', [0])
    ymax = array.array('d', [0])
    c.GetRangeAxis(xmin, ymin, xmax, ymax)
    #print(xmin, ymin, xmax, ymax)
    # Create a second axis
    second_axis = ROOT.TGaxis(xmin[0], ymax[0], 
                              xmax[0], ymax[0], 
                              -2/ _scale, 2 / _scale, 
                              405, "-")
    second_axis.SetTitle(second_axis_title)
    # second_axis.SetTitle("#DeltaN_{#pi-e}/N_{#pi-e} [%]")
    # second_axis.SetTitle("#DeltaN_{#pi-#mu-e}/N_{#pi-mu-e} [%]")
    second_axis.SetLineColor(second_axis_color)
    second_axis.SetLabelColor(second_axis_color)
    second_axis.SetTitleColor(second_axis_color)
    second_axis.Draw()
    
    line = ROOT.TLine()
    line.SetLineStyle(ROOT.kDashed)
    line.DrawLine(-1, 0, -1, len(_ranked_nps) + _offset)
    line.DrawLine(1, 0, 1, len(_ranked_nps) + _offset)
    leg = ROOT.TLegend(c.GetLeftMargin() , c.GetBottomMargin()+ 0.01, 1, c.GetBottomMargin() + 0.1)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.AddEntry(gr_np_pull_and_constraint, "Pull", 'lp')
    # leg.AddEntry(gr_poi_impact, '#pm 1#sigma impact on R_{e/#mu}', 'f')
    leg.AddEntry(gr_poi_impact, '#pm 1#sigma impact', 'f')
    leg.Draw()

    plotables = [leg, line, second_axis, gr_np_pull_and_constraint, gr_poi_impact, h_template]
    return c, plotables


def correlation_plot(res):

    c = ROOT.TCanvas()
    c.SetLeftMargin(0.3)
    c.SetBottomMargin(0.4)
    c.SetRightMargin(0.15)
    h = res.correlationHist()
    h.LabelsOption("V")
    ROOT.gStyle.SetPaintTextFormat('1.2f')
    h.Draw("colzTEXT01")

    return c, h

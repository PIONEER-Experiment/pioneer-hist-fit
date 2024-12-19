import ROOT


def make_measurement(do_he=True, do_le=True, do_tf_he=True, do_tf_le=True):
    InputFile = 'output/fit_input.root'
    meas = ROOT.RooStats.HistFactory.Measurement("meas", "meas")
    meas.SetOutputFilePrefix( "./output/" )
    
    meas.SetLumi( 1.0 )
    _lumi_rel_err = 1e-7
    meas.SetLumiRelErr( _lumi_rel_err )
    meas.SetExportOnly( True )
    

    if do_he:
        # Create a channel
        chan_he = ROOT.RooStats.HistFactory.Channel( "high_energy" )
        chan_he.SetData( "data_he_rebin", InputFile )
        chan_he.SetStatErrorConfig(0.0001, 'Poisson')
    
        pienu_he = ROOT.RooStats.HistFactory.Sample("pienu_he", "pienu_he", InputFile)
        pienu_he.SetNormalizeByTheory(False)
        pienu_he.AddNormFactor( "pie_HE", 1, 0., 4)
        chan_he.AddSample( pienu_he )
        
        michel_he = ROOT.RooStats.HistFactory.Sample("michel_he", "michel_he", InputFile)
        michel_he.SetNormalizeByTheory(False)
        michel_he.AddNormFactor( "pidar_mudar_HE", 1, 0., 4)
        chan_he.AddSample( michel_he )
        
        mudif_he = ROOT.RooStats.HistFactory.Sample("mudif_he", "mudif_he", InputFile)
        mudif_he.SetNormalizeByTheory(True)
        # mudif_he.AddNormFactor( "mudif_HE",  1, 0.95, 1.05 )
        mudif_he.AddOverallSys( "mudif_HE",  0.9, 1.1)
        chan_he.AddSample( mudif_he )
        
        pileup_he = ROOT.RooStats.HistFactory.Sample("pileup_he", "pileup_he", InputFile)
        pileup_he.SetNormalizeByTheory(True)
        # pileup_he.AddNormFactor('pileup_HE', 1, 0.9, 1.1)
        pileup_he.AddOverallSys('pileup_HE', 0.999, 1.001)
        chan_he.AddSample( pileup_he )

        rbm_he = ROOT.RooStats.HistFactory.Sample("rbm_he", "rbm_he", InputFile)
        rbm_he.SetNormalizeByTheory(True)
        rbm_he.AddOverallSys('rbm_HE', 0.999, 1.001)
        # rbm_he.AddNormFactor('rbm_HE', 1, 0.9, 2)
        chan_he.AddSample( rbm_he)

    
        meas.AddChannel( chan_he )

    if do_le:
        # Create a channel
        chan_le = ROOT.RooStats.HistFactory.Channel( "low_energy" )
        chan_le.SetData( "data_le_rebin", InputFile )
        chan_le.SetStatErrorConfig( 0.0001, "Poisson" )
        
        pienu_le = ROOT.RooStats.HistFactory.Sample("pienu_le", "pienu_le", InputFile)
        pienu_le.SetNormalizeByTheory(True)
        # pienu_le.AddNormFactor( "pie_LE", 1, 0., 4)
        # pienu_le.AddOverallSys( "pie_LE", 0.999, 1.001)
        chan_le.AddSample( pienu_le )
        
        michel_le = ROOT.RooStats.HistFactory.Sample("michel_le", "michel_le", InputFile)
        michel_le.SetNormalizeByTheory(False)
        michel_le.AddNormFactor( "pidar_mudar_LE", 1, 0, 4)
        chan_le.AddSample( michel_le )
        
        mudif_le = ROOT.RooStats.HistFactory.Sample("mudif_le", "mudif_le", InputFile)
        mudif_le.SetNormalizeByTheory(True)
        # mudif_le.AddNormFactor( "mudif_LE", 1., 0.5, 1.5)
        # mudif_le.AddOverallSys( "mudif_LE", 0.999, 1.001)
        chan_le.AddSample( mudif_le )
        
        pidif_le = ROOT.RooStats.HistFactory.Sample("pidif_le", "pidif_le", InputFile)
        pidif_le.SetNormalizeByTheory(True)
        # pidif_le.AddNormFactor('pidif_LE', 1, 0, 4)
        pidif_le.AddOverallSys('pidif_LE', 0.99, 1.01)
        chan_le.AddSample( pidif_le)

        pileup_le = ROOT.RooStats.HistFactory.Sample("pileup_le", "pileup_le", InputFile)
        pileup_le.SetNormalizeByTheory(True)
        # pileup_le.AddNormFactor('pileup_LE', 1, 0, 4)
        # pileup_le.AddOverallSys('pileup_LE', 0.999, 1.001)
        chan_le.AddSample( pileup_le)

        rbm_le = ROOT.RooStats.HistFactory.Sample("rbm_le", "rbm_le", InputFile)
        rbm_le.SetNormalizeByTheory(True)
        # rbm_le.AddOverallSys('rbm_LE', 0.999, 1.001)
        # rbm_le.AddNormFactor('rbm_LE', 1, 0.9, 4)
        chan_le.AddSample( rbm_le)

        meas.AddChannel( chan_le )
    
    if do_tf_le:
        # Create a channel
        chan_tf_le = ROOT.RooStats.HistFactory.Channel( "low_energy_tail_fraction" )
        chan_tf_le.SetStatErrorConfig(0.0001, "Poisson" )
        chan_tf_le.SetData( "data_tf_le", InputFile )
        
        pienu_tf_le = ROOT.RooStats.HistFactory.Sample("pienu_tf_le", "pienu_tf_le", InputFile)
        pienu_tf_le.SetNormalizeByTheory(False)
        pienu_tf_le.AddNormFactor( "pie_tf_LE", 1, 0, 4)
        chan_tf_le.AddSample( pienu_tf_le )
        
        michel_tf_le = ROOT.RooStats.HistFactory.Sample("michel_tf_le", "michel_tf_le", InputFile)
        michel_tf_le.SetNormalizeByTheory(True)
        michel_tf_le.AddOverallSys('pidar_mudar_tf_uncert', 0.85, 1.15)
        chan_tf_le.AddSample( michel_tf_le )
        
        mudif_tf_le = ROOT.RooStats.HistFactory.Sample("mudif_tf_le", "mudif_tf_le", InputFile)
        mudif_tf_le.SetNormalizeByTheory(True)
        mudif_tf_le.AddOverallSys('mudif_tf_uncert', 0.85, 1.15)
        chan_tf_le.AddSample( mudif_tf_le )
        
        # pileup_tf_le = ROOT.RooStats.HistFactory.Sample("pileup_tf_le", "pileup_tf_le", InputFile)
        # pileup_tf_le.SetNormalizeByTheory(True)
        # pileup_tf_le.AddOverallSys('pileup_tf_uncert', 0.9999, 1.0001)
        # chan_tf_le.AddSample( pileup_tf_le )
        
        pidif_tf_le = ROOT.RooStats.HistFactory.Sample("pidif_tf_le", "pidif_tf_le", InputFile)
        pidif_tf_le.SetNormalizeByTheory(True)
        pidif_tf_le.AddOverallSys('pidif_tf_uncert', 0.85, 1.15)
        chan_tf_le.AddSample( pidif_tf_le )

        rbm_tf_le = ROOT.RooStats.HistFactory.Sample("rbm_tf_le", "rbm_tf_le", InputFile)
        rbm_tf_le.SetNormalizeByTheory(True)
        rbm_tf_le.AddOverallSys('rbm_tf_uncert', 0.85, 1.15)
        chan_tf_le.AddSample( rbm_tf_le)
        
        meas.AddChannel( chan_tf_le)
        
    if do_tf_he:
        # Create a channel
        chan_tf_he = ROOT.RooStats.HistFactory.Channel( "high_energy_tail_fraction" )
        chan_tf_he.SetData( "data_tf_he", InputFile )
        chan_tf_he.SetStatErrorConfig( 0.0001, "Poisson" )
    
        # Now, create some samples
        pienu_tf_he = ROOT.RooStats.HistFactory.Sample("pienu_tf_he", "pienu_tf_he", InputFile)
        pienu_tf_he.SetNormalizeByTheory(False)
        pienu_tf_he.AddNormFactor( "pie_tf_HE", 1, 0, 4)
        chan_tf_he.AddSample( pienu_tf_he )

        michel_tf_he = ROOT.RooStats.HistFactory.Sample("michel_tf_he", "michel_tf_he", InputFile)
        michel_tf_he.SetNormalizeByTheory(True)
        michel_tf_he.AddOverallSys('pidar_mudar_tf_uncert', 0.85, 1.15)
        chan_tf_he.AddSample( michel_tf_he )

        mudif_tf_he = ROOT.RooStats.HistFactory.Sample("mudif_tf_he", "mudif_tf_he", InputFile)
        mudif_tf_he.SetNormalizeByTheory(True)
        mudif_tf_he.AddOverallSys('mudif_tf_uncert', 0.85, 1.15)
        chan_tf_he.AddSample( mudif_tf_he )

        pileup_tf_he = ROOT.RooStats.HistFactory.Sample("pileup_tf_he", "pileup_tf_he", InputFile)
        pileup_tf_he.SetNormalizeByTheory(True)
        pileup_tf_he.AddOverallSys('pileup_tf_uncert', 0.85, 1.15)
        chan_tf_he.AddSample( pileup_tf_he )

        meas.AddChannel( chan_tf_he)

    # set poi
    if do_he:
        meas.SetPOI( 'pie_HE' )
    elif do_le:
        meas.SetPOI( 'pidar_mudar_LE' )
    elif do_tf_he:
        meas.SetPOI( 'pie_tf_HE')
    elif do_tf_le:
        meas.SetPOI( 'pie_tf_LE')
    else:
        raise ValueError
        
        
    
    # for channel in meas.GetChannels():
    #     for sample in channel.GetSamples():
    #         sample.ActivateStatError()
    
    # Collect the histograms from their files,
    # print some output,
    meas.CollectHistograms()
    #meas.PrintTree()
    # meas.PrintXML( "xmlFromPy", meas.GetOutputFilePrefix() )
    workspace = ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast( meas )

    rfile = ROOT.TFile('output/measurement.root', 'recreate')
    meas.writeToFile(rfile)
    workspace.Write()
    rfile.Close()
    # return meas

if __name__ == '__main__':
    make_measurement()

import tabulate
import array
import math
import ROOT

from .helper import get_time_corrections

def ratio_absolute_error(num, den, err_num, err_den):

    rel_err_num_sq = 0
    if num != 0:
        rel_err_num_sq = math.pow(err_num / num, 2)
    rel_err_den_sq = 0
    if den != 0:
        rel_err_den_sq = math.pow(err_den / den , 2)
    rel_err = math.sqrt(rel_err_num_sq + rel_err_den_sq)
    if den == 0:
        return 0.

    return num / den * rel_err

    

def tail_fraction_report(pienu_all, pienu_tf, energy_binning):

    # threshold_ene = 57.5
    # bin_edge = pienu_all.GetYaxis().FindBin(threshold_ene)
    threshold_ene = pienu_all.GetYaxis().GetBinLowEdge(energy_binning[1])
    err_all_den = array.array('d', [0])
    err_all_num = array.array('d', [0])
    
    pienu_all_den = pienu_all.IntegralAndError(
        0, pienu_all.GetNbinsX() + 1,
        energy_binning[0], energy_binning[2],
        err_all_den)

    pienu_all_num = pienu_all.IntegralAndError(
        0, pienu_all.GetNbinsX() + 1,
        energy_binning[0], energy_binning[1],
        # 0, bin_edge,
        err_all_num)

    err_tf_den = array.array('d', [0])
    err_tf_num = array.array('d', [0])

    pienu_tf_den = pienu_tf.IntegralAndError(
        0, pienu_tf.GetNbinsX() + 1,
        energy_binning[0], energy_binning[2],
        # 0, pienu_tf.GetNbinsY() + 1,
        err_tf_den)

    pienu_tf_num = pienu_tf.IntegralAndError(
        0, pienu_tf.GetNbinsX() + 1,
        energy_binning[0], energy_binning[1],
        # 0, bin_edge,
        err_tf_num)


    tf_all = pienu_all_num / pienu_all_den
    tf_all_err = ratio_absolute_error(
        pienu_all_num, pienu_all_den, err_all_num[0], err_all_den[0])

    tf_tf = pienu_tf_num / pienu_tf_den
    tf_tf_err = ratio_absolute_error(
        pienu_tf_num, pienu_tf_den, err_tf_num[0], err_tf_den[0])

    ratio = tf_tf / tf_all
    err_ratio = ratio_absolute_error(tf_tf, tf_all, tf_tf_err, tf_all_err)
    table = [
        ['True Value', 100 * tf_all, 100 * tf_all_err],
        ['Tail Fraction Analysis Region', 100 * tf_tf, 100 * tf_tf_err],
        ['Analysis Region / True Value', 100 * ratio, 100 * err_ratio]
        ]

    headers = ['Tail Fraction Calculation', 'Value (%)', 'MC stat uncertainty']

    print('TAIL FRACTION REPORT with Energy threshold = {} MeV'.format(threshold_ene))
    print (tabulate.tabulate(table, headers=headers, floatfmt='1.5f'))
    print()
    print('TF unknown correction = ', 1. /ratio)
    print ('------------------------------------------------')

def yield_table_report(hists_tf_le, scale=1):
    _samples = []
    _yields = []
    _uncert = []
    for k, v in hists_tf_le.items():
        # if k == 'beam_muons':
        #     continue
        _samples += [k]
        err = array.array('d', [0])
        if isinstance(v, (ROOT.TH1F, ROOT.TH1D)):
            _yield = v.IntegralAndError(0, v.GetNbinsX() + 1, err)
        elif isinstance(v, ROOT.TH2):
            _yield = v.IntegralAndError(0, v.GetNbinsX() + 1, 0, v.GetNbinsY() + 1, err)
        else:
            raise TypeError(type(v))
        _yields += [_yield]
        _uncert += [err[0]]

    _tot = sum(_yields)

    table = []
    for _s, _y, _u in zip(_samples, _yields, _uncert):
        if _y == 0:
            _uncert = '<{:1.2f} (95% CL)'.format(3 * scale)
        else:
            _uncert = '{:1.2f} ({:1.2f})%'.format(
                _u, 100 * (_u /_y))
        _frac = 100 * _y/_tot
        table.append([
            _s, _y, _uncert, _frac])
    print()
    print (tabulate.tabulate(table, headers=['Sample', 'Yields', 'Uncert. (Rel. uncert. in %)', 'Composition (%)'], floatfmt=['1.2f', '1.2f', '1.2f', '1.4f']))
    print()

def acceptance_report(hist_pienu, hist_michel):

    T0 = 5 #ns
    T1 = 500 #ns
    frac_pienu, frac_michel, r_time = get_time_corrections(T0, T1)

    nbins_x = hist_pienu.GetNbinsX()
    nbins_y = hist_pienu.GetNbinsY()

    err_pienu = array.array('d', [0])
    err_michel = array.array('d', [0])
    n_pienu = hist_pienu.IntegralAndError(0, nbins_x, 0, nbins_y, err_pienu)
    n_michel= hist_michel.IntegralAndError(0, nbins_x, 0, nbins_y, err_michel)

    n_pienu_corr = n_pienu / frac_pienu
    m_michel_corr = n_michel / frac_michel

    rpi_sm = 1.23524
    err_rpi_sm = 0.00015
    rpi_pienu = 1.2327
    err_rpi_pienu = 0.023
    rpi_pioneer = n_pienu / n_michel * r_time
    err_rpi_pioneer = ratio_absolute_error(n_pienu, n_michel, err_pienu[0], err_michel[0])

    ratio = rpi_pioneer * 1e4 / rpi_sm
    err_ratio = ratio_absolute_error(rpi_pioneer * 1e4, rpi_sm, err_rpi_pioneer * 1e4, err_rpi_sm)
    table = [
        ['N(pienu events)', n_pienu, err_pienu[0]],
        ['N(pimue events)', n_michel, err_michel[0]],
        tabulate.SEPARATING_LINE,
        ['pienu time corr [{}, {}] ns'.format(T0, T1), frac_pienu, None],
        ['pimue time corr [{}, {}] ns'.format(T0, T1), frac_michel, None],
        ['R_time', r_time, None],
        tabulate.SEPARATING_LINE,
        ['R(e/mu) (True PIONEER) (X 1e4)', rpi_pioneer*1e4, err_rpi_pioneer*1e4 ],
        ['R(e/mu) (SM) (X 1e4)', rpi_sm, err_rpi_sm],
        ['R(e/mu) (PIENU EXP.) (X 1e4)', rpi_pienu, err_rpi_pienu],
        ['R(e/mu) (True PIONEER) / R(e/mu) (SM)', ratio, err_ratio]
        ]

    print()
    print('ACCEPTANCE REPORT')
    print (tabulate.tabulate(table, headers=['Quantity', 'Value', 'Uncert'], floatfmt=('', '1.4f', '1.5f')))
    print()
    print ('R_time = ', r_time)
    print ('R_unknown = ', 1. / ratio)


def fitresult_report(r_pi, r_pi_notf, tf, res):
    table = [
        [
            'PIONEER',
            r_pi.getVal()*1e4, 
            r_pi.getPropagatedError(res)*1e4,
            r_pi.getPropagatedError(res) / r_pi.getVal()*100,
        ],
        [
            'PIONEER (no TF)',
            r_pi_notf.getVal()*1e4, 
            r_pi_notf.getPropagatedError(res)*1e4,
            r_pi_notf.getPropagatedError(res) / r_pi.getVal()*100,
        ],
        
        ['pienu', 1.23270, 0.00230, 0.00230/1.23270 * 100],
        ['sm',  1.23524, 0.00015, 0.00015/1.23524 * 100],
    ]
    
    print (tabulate.tabulate(
        table, headers=[
            'Source', 'R_pi x 1e4', 'Delta(R_pi) x 1e4', 'Delta(R_pi) / R_pi (%)']))
    print(tabulate.tabulate(
        [        [
            'Tail Fraction',
            tf.getVal()*1e2,
            tf.getPropagatedError(res)*1e2,
            tf.getPropagatedError(res)/tf.getVal() * 1e2
        ]],
        headers=['Quantity x 1e2', 'Uncertainty x 1e2', 'Relative Uncertainty (%)']))

    

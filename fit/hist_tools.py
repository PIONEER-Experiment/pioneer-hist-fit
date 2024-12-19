import ROOT

def get_max(hist_list):
    if isinstance(hist_list, (list, tuple)):
        return max([h.GetBinContent(h.GetMaximumBin()) for h in hist_list])
    elif isinstance(hist_list, dict):
        return max([h.GetBinContent(h.GetMaximumBin()) for _, h in hist_list.items()])
        


def sum_hist(hist_list, name=None):
    if isinstance(hist_list, dict):
        _hist_list = [v for _, v in hist_list.items()]
    else:
        _hist_list = hist_list

    h_tot = _hist_list[0].Clone()
    if name is None:
        h_tot.SetName('tot')
    else:
        h_tot.SetName(name)
    for h in _hist_list[1:]:
        h_tot.Add(h)
    return h_tot

def reset_errors(h, remove=False, poisson=True, prescale=1):
    for ibin in range(h.GetNbinsX()):
        if remove:
            h.SetBinError(ibin, 0)
        elif poisson:
            h.SetBinError(ibin, ROOT.sqrt(prescale * h.GetBinContent(ibin)))
        else:
            raise ValueError

def convert_to_counting(h):
    for ibin in range(h.GetNbinsX()+1):
        h.SetBinContent(ibin, int(h.GetBinContent(ibin)))

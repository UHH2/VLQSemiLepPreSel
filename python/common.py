import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')


import varial
import varial.history
import varial.tools
import varial.analysis
import varial.generators as gen
import glob
import os
import itertools
import cPickle


def label_axes(wrps):
    for w in wrps:
        if 'TH1' in w.type:
            w.histo.GetXaxis().SetTitle(w.histo.GetTitle())
            w.histo.GetYaxis().SetTitle('Events')
            w.histo.SetTitle('')
        yield w


signal_indicators = ['_TH_', 'TpTp_',]


def get_samplename(wrp):
    if hasattr(wrp, 'sample') and wrp.sample:
        return wrp.sample
    fname = os.path.basename(wrp.file_path)
    if fname.startswith('uhh2'):
        return fname.split('.')[-2]
    else:
        return os.path.splitext(fname)[0]


def get_legend(wrp, sig_ind):
    smpl = get_samplename(wrp)
    if 'Run20' in smpl:
        return 'Data'
    elif any(s in smpl for s in sig_ind):
        mass = smpl.split('_')[-1]
        if mass.startswith('M'):
            mass = mass[1:]
        hnd = 'rh' if '_RH_' in smpl else 'lh'
        return 'T_{%s}(%d)#rightarrowtH' % (hnd, int(mass))
    else:
        return smpl


def get_sys_info(w):
    def get_info(tok):
        if tok in w.file_path:
            return next(s
                        for s in w.file_path.split('/') 
                        if s.endswith(tok))
        else:
            return ''
    return get_info('__minus') or get_info('__plus')


def add_wrp_info(wrps, sig_ind=None):
    sig_ind = sig_ind or signal_indicators
    return varial.generators.gen_add_wrp_info(
        wrps,
        sample=get_samplename,
        legend=lambda w: get_legend(w, sig_ind),
        is_signal=lambda w: any(s in w.sample for s in sig_ind),
        is_data=lambda w: 'Run20' in w.sample,
        variable=lambda w: w.in_file_path.split('/')[-1],
        sys_info=get_sys_info,
    )


def merge_decay_channels(wrps, postfixes=('_Tlep', '_NonTlep'), suffix='', print_warning=True, yield_orig=False):
    """histos must be sorted!!"""

    @varial.history.track_history
    def merge_decay_channel(w):
        return w

    def do_merging(buf):
        res = varial.operations.merge(buf)
        res.sample = next(res.sample[:-len(p)]
                          for p in postfixes
                          if res.sample.endswith(p))+suffix
        try:
            res.legend = next(res.legend[:-len(p)]
                              for p in postfixes
                              if res.legend.endswith(p))+suffix
        except StopIteration:
            pass
        res.file_path = ''
        del buf[:]
        return merge_decay_channel(res)

    buf = []
    for w in wrps:
        if any(w.sample.endswith(p) for p in postfixes):
            buf.append(w)
            if len(buf) == len(postfixes):
                yield do_merging(buf)
            if yield_orig:
                yield w
        else:
            if buf:
                if print_warning:
                    print 'WARNING In merge_decay_channels: buffer not empty.\n' \
                          'postfixes:\n' + str(postfixes) + '\n' \
                          'Flushing remaining items:\n' + '\n'.join(
                        '%s, %s' % (w.sample, w.in_file_path) for w in buf
                    )
                yield do_merging(buf)
            yield w
    if buf:
        yield do_merging(buf)

def is_signal(name):
    return any(t in name for t in signal_indicators)


def yield_n_objs(wrps, n=50):
    i = 0
    for w in wrps:
        i += 1
        if i < n:
            yield w
        else:
            break


######################################################### limit calculation ###
try:
    from varial.extensions.limits import ThetaLimits
except ImportError:
    from varial import monitor
    monitor.message(
        'UHH2.VLQSemiLepPreSel.common',
        'WARNING Theta is not working'
    )
    ThetaLimits = object


class TriangleLimitPlots(varial.tools.Tool):
    def __init__(self,
        name=None,
        limit_rel_path=''
    ):
        super(TriangleLimitPlots, self).__init__(name)
        self.limit_rel_path = limit_rel_path

    # def make_tri_hist(self, wrps, mass_ind):
    #     tri_hist = ROOT.TH2F("triangular_limits", ";br to th;br to tz", 11, -0.05, 1.05, 11, -0.05, 1.05)
    #     for w in wrps:
    #         br_th = float(w.brs['th'])
    #         br_tz = float(w.brs['tz'])
    #         tri_hist.Fill(br_th, br_tz, w.res_exp_y[mass_ind])
    #     return varial.wrappers.HistoWrapper(tri_hist,
    #             legend='M-'+str(w.masses[mass_ind]),
    #             mass=w.masses[mass_ind])


    def run(self):
        if isinstance(self.limit_rel_path, str):
            self.limit_rel_path = [self.limit_rel_path]
        theta_tools = []
        for lim_path in self.limit_rel_path:
            if lim_path.startswith('..'):
                theta_tools += glob.glob(os.path.join(self.cwd, lim_path))
            else:
                theta_tools += glob.glob(lim_path)
        self.wrps = list(self.lookup_result(k) for k in theta_tools)
        if any(not a for a in self.wrps):
            self.wrps = gen.dir_content(self.limit_rel_path, 'result.info', 'result')
        # theta_tools = glob.glob(os.path.join(self.cwd+'..', self.limit_rel_path))
        # print theta_tools
        # self.wrps = list(self.lookup_result(k) for k in theta_tools)
        # print self.wrps
        # self.filename = os.path.join(varial.analysis.cwd, self.name + ".root")

class TriangleCSLimitPlots(TriangleLimitPlots):
    def __init__(self,
        name=None,
        limit_rel_path=''
    ):
        super(TriangleCSLimitPlots, self).__init__(name, limit_rel_path)

    def make_tri_hist(self, wrps, mass_ind):
        tri_hist = ROOT.TH2F("triangular_limits", ";br to th;br to tz", 21, -0.025, 1.025, 21, -0.025, 1.025)
        for w in wrps:
            br_th = float(w.brs['h'])
            br_tz = float(w.brs['z'])
            tri_hist.Fill(br_th, br_tz, w.res_exp_y[mass_ind])
        return varial.wrappers.HistoWrapper(tri_hist,
                legend='M-'+str(w.masses[mass_ind]),
                mass=w.masses[mass_ind])


    def run(self):
        super(TriangleMassLimitPlots, self).run()
        # filename = os.path.join(varial.analysis.cwd, self.name + ".root")
        # f = ROOT.TFile.Open(filename, "RECREATE")
        # f.cd()
        list_hists=[]
        for i, m in enumerate(wrps[0].masses):
            list_hists.append(self.make_tri_hist(wrps, i))
        # tri_hist.Write()
        self.result = list_hists
        # f.Close()

class TriangleMassLimitPlots(TriangleLimitPlots):
    def __init__(self,
        name=None,
        limit_rel_path='',
        leg_x='BR(T #rightarrow tH)',
        leg_y='BR(T #rightarrow tZ)',
    ):
        super(TriangleMassLimitPlots, self).__init__(name, limit_rel_path)
        self.leg_x = leg_x
        self.leg_y = leg_y

    def make_tri_hist(self, wrps):
        tri_hist_exp = ROOT.TH2F('triangular_limits_exp', ';'+self.leg_x+';'+self.leg_y+';95% expected T quark mass limit (GeV)', 6, -0.1, 1.1, 6, -0.1, 1.1)
        tri_hist_obs = ROOT.TH2F('triangular_limits_obs', ';'+self.leg_x+';'+self.leg_y+';95% observed T quark mass limit (GeV)', 6, -0.1, 1.1, 6, -0.1, 1.1)
        tri_hist_exp.SetContour(15)
        tri_hist_obs.SetContour(15)
        tri_hist_exp.SetMinimum(700.)
        tri_hist_obs.SetMinimum(700.)
        if not all(isinstance(w, varial.wrappers.WrapperWrapper) for w in wrps):
            wrps = itertools.ifilter(lambda w: not w.legend == 'Theory', wrps)
            wrps = sorted(wrps, key=lambda w: w.save_name)
            wrps = gen.group(wrps, lambda w: w.save_name)
        for w in wrps:
            wrp = w[0]
            br_th = float(wrp.brs['h'])
            br_bw = float(wrp.brs['w'])
            if wrp.exp_mass_excl:
                lim_exp = float(wrp.exp_mass_excl)
            else:
                lim_exp = 1e-19
            if wrp.obs_mass_excl:
                lim_obs = float(wrp.obs_mass_excl)
            else:
                lim_obs = 1e-19
            tri_hist_exp.Fill(br_th, br_bw, lim_exp)
            tri_hist_obs.Fill(br_th, br_bw, lim_obs)
        return varial.wrappers.HistoWrapper(tri_hist_exp, save_name='lim_exp'), varial.wrappers.HistoWrapper(tri_hist_obs, save_name='lim_obs')


    def run(self):
        super(TriangleMassLimitPlots, self).run()
        list_hists = list(self.make_tri_hist(self.wrps))
        # tri_hist.Write()
        self.result = list_hists
        # f.Close()


class CrossSectionTables(varial.tools.Tool):
    def __init__(self,
        name=None,
        cs_rel_path='',
        limit_rel_path='',
        mass_rel_path='',
        br_combos=None,
        mass_points=None
    ):
        super(CrossSectionTables, self).__init__(name)
        self.limit_rel_path = limit_rel_path
        self.mass_rel_path = mass_rel_path
        self.cs_rel_path = cs_rel_path
        self.br_combos = br_combos
        self.mass_points = mass_points


    def run(self):
        rel_path = self.cs_rel_path
        if self.cs_rel_path.startswith('..'):
            rel_path = os.path.join(self.cwd, self.cs_rel_path)
        limits = os.listdir(rel_path)
        limits = list(itertools.ifilter(lambda w: w in list(i for i, _ in self.br_combos), limits))

        lines_exp = []
        lines_exp.append(r"\begin{tabular}{|l "
            + len(limits)*"| r "
            + r"|}\hline")
        lines_exp.append("Mass/BRs & " + " & ".join(br for _, br in self.br_combos)
            + r"\\ \hline")
        # lines_obs = list(lines_exp)

        for m in self.mass_points:
            line_exp = varial.analysis.get_pretty_name(m) + ' '
            # line_obs = varial.analysis.get_pretty_name(m) + ' '
            for lim in limits:
                line_exp += "& "
                # line_obs += "& "
                theta_path = os.path.join(rel_path, lim, self.mass_rel_path, m, self.limit_rel_path)
                theta_res = self.lookup_result(theta_path)
                if not theta_res:
                    theta_res = gen.dir_content(theta_path, 'result.info', 'result')
                assert len(theta_res) == 1, 'Too many Theta results!'
                exp_lim = cPickle.loads(theta_res[0].res_exp)
                obs_lim = cPickle.loads(theta_res[0].res_obs)
                line_exp += "%.2g (%.2g)" % (obs_lim.y[0], exp_lim.y[0])
                # line_obs += "%.2g" % obs_lim.y[0]
            line_exp += r" \\"
            # line_obs += r" \\"
            lines_exp.append(line_exp)
            # lines_obs.append(line_obs)

        lines_exp.append(r"\hline")
        # lines_obs.append(r"\hline")
        lines_exp.append(r"\end{tabular}")
        # lines_obs.append(r"\end{tabular}")
        lines_exp = '\n'.join(lines_exp)
        # lines_obs = '\n'.join(lines_obs)

        with open(self.cwd+'cs_limits.tex', 'w') as f:
            f.write(lines_exp)

        # with open(self.cwd+'obs_cs_limits.tex', 'w') as f:
        #     f.write(lines_obs)

        # theta_tools = glob.glob(os.path.join(self.cwd+'..', self.limit_rel_path))
        # print theta_tools
        # self.wrps = list(self.lookup_result(k) for k in theta_tools)
        # print self.wrps
        # self.filename = os.path.join(varial.analysis.cwd, self.name + ".root")


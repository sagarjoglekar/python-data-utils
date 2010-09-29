from pylab import *
from stats import *

import math
import random
import bisect

import collections
import itertools

from numpy.core import mean,std
from numpy.lib import median 

import read_objects as ro
import long_tail as lt

linestyles = lambda : itertools.cycle(["-","--",":","-."])
linemarkers = lambda: itertools.cycle(['+', ',', '.','1','2', '3', '4'])
no_lstyles = no_lmarkers = itertools.cycle([""])

def write_line(filnam,datapts,xlbl="X",ylbl="Y",comment=None):
    import csv
    #convert from [xvec,yvec] vec to [(xval,yval)...]
    if len(datapts) == 2: datapts = zip(*datapts) 
    with open(filnam,'wb') as f:
        w = csv.writer(f)
        if comment: w.writerow(["#" + comment,"N=%s points" % len(datapts)])
        w.writerow(["#" + xlbl,ylbl])
        for x,y in datapts:
            w.writerow([x,y])

def read_line(filnam):
    import csv
    with open(filnam) as f:
        r = csv.reader(f)
        x,y = list(),list()
        for _x,_y in r:
            if _x.startswith("#"): continue
            x.append(_x)
            y.append(_y)
        return x,y

def fit_exponential(samples):
    from rpy import r
    samples = [double(n) for n in samples]#because rpy does not like longs!
    r.library('MASS')
    f = r.fitdistr(samples,'exponential')
    rat = f['estimate']['rate']
    qp = r.qexp(r.ppoints(samples),rate=rat)
    
    return qp, rat

def qq_plot(samples,qp):
    samples = sorted(samples)
    lim = max(samples)
    setup_fig("predicted","actual")
    plot([0,lim],[0,lim],c='k')
    
    x_y = zip(qp,samples)
    x_y = filter(lambda xy: xy[0] < lim and xy[1] < lim, x_y)
    step_siz = len(x_y)/50000 #take at most 50k points!
    if step_siz > 1:
        x_y = x_y[::step_siz]
    x,y = zip(*x_y)

    x,y = zip(*filter(lambda xy: xy[0] < lim and xy[1] < lim, x_y))
    
    plot(x,y,'k,')
    xlim(0,lim); ylim(0,lim)

def qq_plot_labels(samples,qp,lgnd,lstyle, hold=False):
    samples = sorted(samples)
    lim = max(samples)
    if not hold:
        setup_fig("predicted","actual")
        plot([0,lim],[0,lim],c='k')
    
    x_y = zip(qp,samples)
    x_y = filter(lambda xy: xy[0] < lim and xy[1] < lim, x_y)
    step_siz = len(x_y)/50000
    if step_siz > 1:
        x_y = x_y[::step_siz]
    x,y = zip(*x_y)
    plot(x,y,lstyle, label=lgnd)
    xlim(0,lim); ylim(0,lim)

def fit_poisson(samples):
    from rpy import r
    r.library('MASS')
    f = r.fitdistr(samples,'poisson')
    l = f['estimate']['lambda'] #predicted mean
    qp = r.qpois(r.ppoints(samples),l)
    return qp,l

def fit_nbinom(samples):
    from rpy import r
    r.library('MASS')
    f = r.fitdistr(samples,'negative binomial')
    s,m = f['estimate']['size'],f['estimate']['mu']
    qp = r.qnbinom(r.ppoints(samples),size=s,mu=m)
    return qp,s,m

def fit_gamma(samples):
    from rpy import r
    samples = [double(n) for n in samples if n > 0]#because rpy does not like longs!
    r.library('MASS')
    f = r.fitdistr(samples,'gamma')
    shap,rat = f['estimate']['shape'],f['estimate']['rate']
    qp = r.qgamma(r.ppoints(samples),shape=shap,rate=rat)
    return qp,shape,rat

def fit_weibull(samples):
    from rpy import r
    #samples = [double(n) for n in samples if n > 0]#because rpy does not like longs!
    r.library('MASS')
    f = r.fitdistr(samples,'weibull')
    sc,sh = f['estimate']['scale'],f['estimate']['shape']
    qp = r.qweibull(r.ppoints(samples),scale=sc,shape=sh)
    return qp,sc,sh


#### helper funcs (DO NOT CORRESPOND TO FIGS in the paper)
def set_timecdf_axes(ylbl="P[T<t]",lgnd=False):
    axis(xmax=86400*7*10)
    xticks(arange(0,6000000,step=86400*7),['wk %d' % d for d in range(10)])
    xlabel("time (weeks)")
    ylabel(ylbl)
    if lgnd:
        legend(loc='lower right')
    show()

## adjusting figuresize (it is 252.0pt for IEEEtran)
def figsize(col_width_pt=252.0):
    """Get fig_width_pt from LaTeX using \showthe\columnwidth"""
    fig_width_pt=col_width_pt
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height =fig_width*golden_mean       # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size
params = {'backend': 'pdf',
#          'text.usetex': True,
          'font.family':'font.sans-serif',
          'axes.labelsize': 9,
          'text.fontsize': 9,
          'legend.fontsize': 9,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': figsize()}
rcParams.update(params)

def setup_fig(xlbl,ylbl):
    figure()
    #axes([0.15,0.25,0.95-0.15,0.95-0.25]) # leave space for labels
    axes([0.15,0.18,0.95-0.15,0.95-0.25]) # leave space for labels
    xlabel(xlbl)
    ylabel(ylbl)

def plt_cdf(qty,lbl,linestyle="-"):
    x,y = ecdf(qty)
    plot(x,y,ls=linestyle, label=lbl,c='k')

def plot_black_white(results,lbls,ls=linestyles(),lm=linemarkers(),lw=1,
                     plot_errorbars=False):
    for x,y,yerr in results:
        if plot_errorbars:
            #x = 0 is not liked by errorbar because of semilogx()
            #so if a fairly obscure exception occurs below,
            #remove points for x=0!
            errorbar(x,y,yerr,marker=lm.next(),label=lbls.next())
        else:
            plot(x,y,linestyle=ls.next(),marker=lm.next(),
                 linewidth=lw,c='k',label=lbls.next())


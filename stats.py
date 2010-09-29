import numpy as np
import random
import operator,itertools
import collections
import bisect

#constants
one_hr = 3600; half_hr = one_hr/2;
half_day = 12*one_hr; one_day = 2*half_day;
one_week = 7*one_day; two_weeks = 2*one_week
    

#random util funcs
def chainlists(listoflists): return list(itertools.chain(*listoflists))

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    try:
        b.next()
    except StopIteration:
        pass
    return itertools.izip(a, b)

sub_successive_elts = lambda liszt: map(operator.sub,liszt[1:],liszt[:-1])

#stats and sampling funcs
def rezero(liszt,newzero): 
    return [x-newzero for x in liszt]

def sample_upto_n(seq,n): return random.sample(seq,min(n,len(seq)))
    
def ecdf(x,num_missing=0):
    """ 
    returns empirical cdf of x, suitable for pylab plot
    use a different num_missing if you know some of the paths didnt complete et
    """
    x2 = np.repeat(np.sort(x),2)
    n = len(x)
    max_fraction = 1-float(num_missing)/(n+num_missing)
    y2 = np.hstack([0.0,np.repeat(np.arange(1,n)/float(n+num_missing),2),max_fraction])
    return x2,y2

def ccdf(x):
    x,y=ecdf(x)
    return x,1-y


def quantiles(x,m_sugg=None):
    """
    separates x into m_sugg quantiles.
    If m_sugg is not suggested, len(x) quantiles are returned
    returns x,quant, where quant[i] \approx fraction of data < x[i]
    """
    def interpolate(m_q,vals,tru_q):
        def xpolate(x):
            pos = bisect.bisect(tru_q,x)
            x1,x2 = tru_q[pos-1],tru_q[pos]
            y1,y2 = vals[pos-1],vals[pos]
            m = float(y2-y1)/(x2-x1)
            c = y1-m*x1
            return m*x + c

        return [xpolate(x) for x in m_q],m_q
    
    x = sorted(x)
    n = len(x)
    tru_q = [(i+0.5)/n for i in range(n)] #==i-0.5/n but i in [0,n-1] now...

    if m_sugg:
        assert m_sugg <= n, "can only interpolate downwards"
        m_q = [(i+0.5)/m_sugg for i in range(m_sugg)]
        x,tru_q = interpolate(m_q,x,tru_q)
        
    return x,tru_q

def repeated_measure(repeats, expt_func, *args):
    results = [expt_func(*args) for i in xrange(repeats)]
    return np.mean(results, 0),2*np.std(results,0)

def repeated_measure_with_alignment(repeats, expt_func, *args):
    return align_results([expt_func(*args) for i in xrange(repeats)])
    
def align_results(results):
    #test if results are of the form (x,y).
    #If so, use x value to align. Else use the order of occurance to align (enumerate the y values to get xvals)
    if not isinstance(results[0][0],tuple):
        results = [[(i,val) for i,val in enumerate(row)] for row in results]
    results = [dict(row) for row in results]
    collected_results = collections.defaultdict(list)
    for row in results:
        for xval in row:
            collected_results[xval].append(row[xval])
    retvals = []
    for xval in collected_results:
        yvals = collected_results[xval]
        retvals.append((xval, np.mean(yvals), 2*np.std(yvals)))
    retvals = retvals[::20]
    return zip(*sorted(retvals,key=lambda x:x[0]))# (x,y,yerr) sorted by x vals

def make_monotone(results,ascending=True):
    """
    only allow increasing values. or decreasing values.
    To be used with align_results or repeated measure_with_alignment
    because sometimes, when aligned and ordered by xvals, things like dratio
    can decrease with increasing number of steps (statistically, neighbour
    to the right of a point  could have lesser dratio because
    it was a different trial). We want to remove outliers like that.

    Typical call: make_monotone(align_results(...),ascending=True)
    """
    comparator=ascending and operator.le or operator.ge
    yvals = operator.itemgetter(2)
    res = zip(*results); res.reverse()
    monotone = [res[0]]
    for r in res:
        if comparator(yvals(r), yvals(monotone[-1])):
            monotone.append(r)
    monotone.reverse()
    return zip(*monotone)

def count_distribution(items,universe_size,normalize=True):
    """
    give a distribution of the form: a random item occurs k times with probability p(k). If normalize=True, then we return random item occurs fraction f of the trace of items with probability p(f).
    items - a sequence of items with repeats.
    universe_size - the total number of items that could occur (some may occur zero times in the list items, so we need this param separately.
    """
    occurs = collections.defaultdict(int)
    for item in items:
        occurs[item] += 1
    norm_factor = normalize and float(len(items)) or 1
    
    occur_cnts = collections.defaultdict(int)
    for item in occurs:
        occur_cnts[occurs[item]/norm_factor] += 1

    universe_size=float(universe_size)
    return zip(*sorted([(oc,oc_cnt/universe_size)
                        for oc,oc_cnt in occur_cnts.items()]))


def stable_rank(liszt,sort_key_func=None):
    ranks = dict()
    prev = None
    prev_rank = 0 #this initial val will not be used
    for cnt,item in enumerate(sorted(liszt,key=sort_key_func)):
        if item != prev:
            prev_rank = cnt
            prev = item
        ranks[item] = prev_rank
    return ranks

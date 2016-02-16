from ROOT import RooStats, Double, RooArgSet, RooFit, RooDataHist, TH1F, gPad
from array import array

# profiled likelihood limit
def plcLimit(obs_, poi_, model, ws, data, CL = 0.95, verbose = False):
    # obs : observable variable or RooArgSet of observables
    # poi : parameter of interest or RooArgSet of parameters
    # model : RooAbsPdf of model to consider including any constraints
    # data : RooAbsData of the data
    # CL : confidence level for interval
    # returns a dictionary with the upper and lower limits for the first/only
    # parameter in poi_ as well as the interval object and status flag
    
    obs = RooArgSet(obs_)
    obs.setName('observables')
    poi = RooArgSet(poi_)
    poi.setName('poi')
    poi.setAttribAll('Constant', False)
    nuis = model.getParameters(obs)
    nuis.remove(poi)
    nuis.remove(nuis.selectByAttrib('Constant', True))
    nuis.setName('nuisance')

    if verbose:
        print 'observables'
        obs.Print('v')
        print 'parameters of interest'
        poi.Print('v')
        print 'nuisance parameters'
        nuis.Print('v')
 
    mc = RooStats.ModelConfig('mc')
    mc.SetWorkspace(ws)
    mc.SetPdf(model)
    mc.SetObservables(obs)
    mc.SetParametersOfInterest(poi)
    mc.SetNuisanceParameters(nuis)

    plc = RooStats.ProfileLikelihoodCalculator(data, mc)
    plc.SetConfidenceLevel(CL)

    interval = plc.GetInterval()

    upperLimit = Double(999.)
    lowerLimit = Double(0.)
    Limits = {}

    paramIter = poi.createIterator()
    param = paramIter.Next()
    while param:
        ok = interval.FindLimits(param, lowerLimit, upperLimit)
        Limits[param.GetName()] = {'ok' : ok, 'upper' : float(upperLimit),
                                   'lower' : float(lowerLimit)}
        param = paramIter.Next()

    if verbose:
        print '%.0f%% CL limits' % (interval.ConfidenceLevel() * 100)
        print Limits

    Limits['interval'] = interval
    return Limits

def expectedPlcLimit(obs_, poi_, model, ws, ntoys = 30, CL = 0.95,
                     binData = False):
    # obs : observable variable or RooArgSet of observables
    # poi : parameter of interest or RooArgSet of parameters
    # model : RooAbsPdf of model to consider including any constraints
    #         the parameters should have the values corresponding to the
    #         background-only hypothesis which will be used to  estimate the
    #         expected limit.
    # ntoys : number of toy datsets to generate to get expected limit
    # CL : confidence level for interval
    # returns a dictionary containing the expected limits and their 1 sigma
    # errors for the first/only parameter in poi_ and a list of the results
    # from the individual toys.

    from math import sqrt
    obs = RooArgSet(obs_)
    obs.setName('observables')
    mPars = model.getParameters(obs)
    genPars = mPars.snapshot()

    print "parameters for generating toy datasets"
    genPars.Print("v")

    limits = []

    upperLimits = TH1F("upperLimits_%s" % poi_.GetName(),
                       "", 100, poi_.getMin(), poi_.getMax())
    lowerLimits = TH1F("lowerLimits_%s" % poi_.GetName(),
                       "", 100, poi_.getMin(), poi_.getMax())
    probs = array('d', [0.022, 0.16, 0.5, 0.84, 0.978])
    upperQs = array('d', [0.]*len(probs))
    lowerQs = array('d', [0.]*len(probs))
    
    for i in range(0,ntoys):
        print 'generate limit of toy %i of %i' % (i+1, ntoys)
        mPars.assignFast(genPars)

        toyData = model.generate(obs, RooFit.Extended())
        if binData:
            toyData = RooDataHist('data_obs_%i' % i, 'data_obs_%i' % i,
                                  obs, toyData)
        toyData.SetName('data_obs_%i' % i)
        toyData.Print()

        limits.append(plcLimit(obs_, poi_, model, ws, toyData, CL))

        #print limits[-1]
        if limits[-1][poi_.GetName()]['ok'] and \
            ((poi_.getMax()-limits[-1][poi_.GetName()]['upper']) > 0.001*poi_.getMax()):
            upperLimits.Fill(limits[-1][poi_.GetName()]['upper'])
        if limits[-1][poi_.GetName()]['ok'] and \
            ((limits[-1][poi_.GetName()]['lower']-poi_.getMin()) > 0.001*abs(poi_.getMin())):
            lowerLimits.Fill(limits[-1][poi_.GetName()]['lower'])

        toyData.IsA().Destructor(toyData)

    upperLimits.GetQuantiles(len(probs), upperQs, probs)
    # upperLimits.Print()
    print 'expected upper limit quantiles using %i toys: ['%(upperLimits.GetEntries()),
    for q in upperQs:
        print '%0.4f' % q,
    print ']'
    lowerLimits.GetQuantiles(len(probs), lowerQs, probs)
    # lowerLimits.Print()
    print 'expected lower limit quantiles using %i toys: [' % (lowerLimits.GetEntries()),
    for q in lowerQs:
        print '%0.4f' % q,
    print ']'
    expLimits = {'upper' : upperQs[2],
                 'upperErr' : sqrt((upperQs[2]-upperQs[1])*(upperQs[3]-upperQs[2])),
                 'lower' : lowerQs[2],
                 'lowerErr' : sqrt((lowerQs[2]-lowerQs[1])*(lowerQs[3]-lowerQs[2])),
                 'ntoys': upperLimits.GetEntries(),
                 'upperQuantiles': upperQs,
                 'lowerQuantiles': lowerQs,
                 'quantiles': probs
                 }
    return (expLimits, limits)

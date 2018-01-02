from htmd.ui import *

queue = LocalGPUQueue()
queue.datadir = './data'
ad = AdaptiveMD()
ad.app = queue
# min number of concurrent jobs
ad.nmin = 1
# max number of concurrent jobs
ad.nmax = 3
# number of epochs to run
ad.nepochs = 50
# set up projection to use for markov modelling
protsel = 'protein and name CA'
ad.projection = MetricSelfDistance(protsel)
ad.updateperiod = 1200 # execute every 20 minutes
ad.run()

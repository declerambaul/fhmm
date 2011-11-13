from time import time


from hmm import *
from db import *
from utils import *


alphabet = [ 0, 1, 2, 3]
realalphabet = ['A','C','G','T']

""" REALISTIC PROBABLITY WOULD BE 0.001 """
pEnter = 0.001
#pEnter = 0.01

""" POSITION WEIGHT MATRIX
hmm/pwm.txt : simple, length 10 pwm
"""
#pwmFile = "hmm/general15.txt"
pwmFile = "hmm/pwm.txt"
#pwmFile = "hmm/pwm1.txt"


backgroundFile = "hmm/background3020.txt"
#backgroundFile = "hmm/background25.txt"

stl = time()

motif = Motif(alphabet,pwmFile,backgroundFile,pEnter)


""" DATABASE  
hmm.1e6 : one million caraters /  hmm/pwm.txt / 30/20/20/30 background (hmm/background3020.txt)  / pEnter = 0.001
nohmm.1e6 : one million caraters /  no hmm / 30/20/20/30 background (hmm/background3020.txt)  / pEnter = 0.00
"""

#db = DB(alphabet,'hmm.1e6')
#db = DB(alphabet,'nohmm.1e6')

tm = []
h = []
db = DB(alphabet,'hmm.1e6')        
index = Index(db,wordlength=8)

for i in arange(0.5,1.4,0.1):
    stl = time()
    fastHMM = FastHMM(motif,index,cutoffWords=i)
    fastHMM.run()    
    tm.append(time()-stl)
    h.append(len(fastHMM.hits))

print h
print tm

""" PLOTTING  """
clf()
plot(arange(0.5,1.4,0.1),tm)
xlabel("Cutoff threshold for index word")
ylabel("Running Time (s)")
title("Time per index word cutoff threshold")

#show()
savefig("graphs/TimeCutoffWord")

clf()
plot(arange(0.5,1.4,0.1),h)
xlabel("Cutoff threshold for index word")
ylabel("Number of matches")
title("Number of matches per index word cutoff threshold")

#show()
savefig("graphs/NoMachtesCutoffWord")



"""
print "There are %d hits" % len(fastHMM.hits)
for hit in fastHMM.hits:     
    print '\t%s : Pos=%d , ScorePWM=%.4f ,  ScoreBG=%.4f  , Ratio=%.4f' % (printWord(db.db[hit[0]:(hit[0]+motif.length)]), hit[0],hit[1],hit[2],hit[1]/hit[2]),
    print
"""
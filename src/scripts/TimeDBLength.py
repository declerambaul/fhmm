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
for i in range(2,7):
    db = DB(alphabet,'hmm.1e'+str(i))        
    index = Index(db,wordlength=8)
    stl = time()
    fastHMM = FastHMM(motif,index)
    fastHMM.run()    
    tm.append(time()-stl)

print tm
clf()
plot(range(2,7),tm)
xlabel("Length of Database (10 ** x")
ylabel("Running Time (s)")
title("Time per Database")
savefig("DecisionList")
#show()
savefig("graphs/TimeDBLength")
""" INDEX """ 



""" SCORE HMM  """
#stl = time()
#fastHMM = FastHMM(motif,index)
#fastHMM.run()
#print 'Running Time for FastHMM : %5fs' % (time()-stl)
#print
#
##print motif
##hits : [indexWords, indexWindow, scorePWM[i],scoreBG[i]]
#print "There are %d hits" % len(fastHMM.hits)
#for hit in fastHMM.hits:     
#    print '\t%s : Pos=%d , ScorePWM=%.4f ,  ScoreBG=%.4f  , Ratio=%.4f' % (printWord(db.db[hit[0]:(hit[0]+motif.length)]), hit[0],hit[1],hit[2],hit[1]/hit[2]),
#    print

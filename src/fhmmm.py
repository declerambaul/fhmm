from time import time


from hmm import *
from db import *
from utils import *

""" BIG SPEED UP IF WORKING WITH [1,2,3,4] """
alphabet = [ 0, 1, 2, 3]
realalphabet = ['A','C','G','T']


""" TRANISTION PROBABILITY TO ENTER MOTIF (REALISTIC VALUE WOULD BE 0.001)"""
pEnter = 0.001
#pEnter = 0.01


""" POSITION WEIGHT MATRIX
hmm/pwm.txt : simple, length 10 pwm
"""
pwmFile = "hmm/pwm.txt"
#pwmFile = "hmm/general15.txt"

""" BACKGROUND  """
backgroundFile = "hmm/background3020.txt"
#backgroundFile = "hmm/background25.txt"

stl = time()

""" LOAD THE MOTIF """
motif = Motif(alphabet,pwmFile,backgroundFile,pEnter)


""" DATABASE  
hmm.1e6 : one million caraters /  hmm/pwm.txt / 30/20/20/30 background (hmm/background3020.txt)  / pEnter = 0.001
nohmm.1e6 : one million caraters /  no hmm / 30/20/20/30 background (hmm/background3020.txt)  / pEnter = 0.001
"""

""" LOAD DATABASE 
db = DB(alphabet,'hmm.1e6')
"""
""" ...OR NEW DATABASE """
db = DB(alphabet,'testDatabase')
db.createDB(motif,length=1e4)


""" DISPLAY STATISTICS 
db.stats()
"""

""" INDEX """ 
index = Index(db,wordlength=8)

print 'Loading Time: %5fs' % (time()-stl)
print


""" DISPLAY MAPPINGS 
for i,seq in enumerate(index.allwords):
    if index.mappings[i] != []:
        print '%d: %s : %s' % (i,seq,index.mappings[i])
"""


""" Fast HMM  """
stl = time()
fastHMM = FastHMM(motif,index)

fastHMM.run()
print 'Running Time for FastHMM : %5fs' % (time()-stl)
print

""" Print Fast HMM Matches  """ 
#print motif
#hits : [indexDB, scorePWM[i],scoreBG[i]]
print "There are %d hits" % len(fastHMM.hits)
for hit in fastHMM.hits:     
    print '\t%s : Pos=%6d , ScorePWM=%.4f ,  ScoreBG=%.4f  , Ratio=%.4f' % (printWord(db.db[hit[0]:(hit[0]+motif.length)]), hit[0],hit[1],hit[2],hit[1]/hit[2]),
    print




""" FORWARD BACKWARD 
stl = time()
print 

forback = ForwardBackward(motif,db)
forback.run()

print 'Running Time for Forward Backward : %5fs' % (time()-stl)
print
#i = scorePWM.index(min(scorePWM))

#print 'Matches found with forward backward at index:'
mfb = []
for i,p in enumerate(forback.pSobs):
    
    if forback.pSobs[i].argmax(0)[0]==1:
        if forback.pSobs[i+1].argmax(0)[0]==2:
            mfb.append(i)        
            #print '\t%s' % i

 """

 
""" Display all argmax of the likelyhood for all states (watch out if db is large)
print 'Index: ',
for i,p in enumerate(forback.pSobs):
        print  '%3s' % i,
print 
print 'State: ',
for i,p in enumerate(forback.pSobs):
        print  '%3s' % p.argmax(0)[0],
"""



""" COMPARE MATCHES FOR FAST HMM AND FORWARD BACKWARD 
mfhmm = [] 
for hit in fastHMM.hits:
    mfhmm.append(hit[0])
mfhmm.sort()


print "There are %d matches for the Fast HMM algorithm" % len(mfhmm)
print "There are %d matches for the Forward Backward algorithm" % len(mfb)
if len(mfhmm) == len(mfb):
    print "Index of matches found"
    for i in range(0,len(mfb)):
        print '\t%2s: FHMM = %6d\tF/B = %6d' %(i,mfhmm[i],mfb[i])
else:
    print mfhmm
    print mfb
    
"""
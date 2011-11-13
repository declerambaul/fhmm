#
#    Fast HMM Mapping
#    Fabian , dec08
#
#    hmm.py : A class representing the motif HMM that we are looking for in the db
#
#

from pylab import *
import numpy 
import scipy 

from hmm import *
from db import *
from utils import *


#alphabet = [ 1, 2, 3, 4]
#alphabet = ['A','C','G','T']

""" 
Class describing a motif in form of an hmm
The hmm has a Background state and an couple of states represented by a Position Weight Matrix.
The transition probability to change from the Background to the pwm should be small, the transition probabilities between the states in the pwm are all equal to one. 
"""
class Motif:
    def __init__(self,alphabet, pwmFile,backgroundFile,pEnter):
        
        self.alphabet = alphabet
        self.pPWM = pEnter
        self.pBG = 1-pEnter
        self.pwm = load(pwmFile,'r')
        self.background = load(backgroundFile,'r')
        self.length = len(self.pwm[0,:])
        
        for s in sum(self.pwm,0):
            if s!=1:
                raise("emission probabilities should sum up to one")
        if sum(self.background)!=1:
            raise("background probabilities should sum up to one")
               

    
    def __str__(self):
        rstr = 'Motif:\n'
        rstr += 'Probabilty to stay in BackGround state: %s\n' % str(self.pBG)
        rstr += 'Probabilty to enter the PWM states: %s\n' % str(self.pPWM) 
        rstr += 'Background:\n'   
        for i, p in enumerate(self.background):
            rstr += ("A", "C", "G", "T")[self.alphabet[i]] + '  '
            rstr += '%s' % p + '\n'
        rstr += 'Position Weight Matrix:\n'
        for i,letter in enumerate(self.alphabet):
            rstr += ("A", "C", "G", "T")[self.alphabet[i]] + '  '
            for p in self.pwm[i]:
                rstr += '%s' % p + '  '
            rstr += '\n'
        return rstr

"""
Implementation of the proposed fast HMM algorithm
"""
class FastHMM:
    def __init__(self,motif,index,cutoffWords=1,cutoffExtension = 1):
        if index.__class__ != Index: raise("Instance of Index class needed for scoring")
        self.index = index
        
        if motif.__class__ != Motif: raise("Instance of Motif class needed for scoring")
        self.motif = motif
        self.motiflength = self.motif.length
        
        #hits: [indexDB, indexPWM, scorePWM, scoreBG]
        self.hits = []
        self.wordlength = self.index.wordlength
        
        # cutoff ratios for p(hmm)/p(backround)
        self.cutoffWords = cutoffWords 
        self.cutoffExtension = cutoffExtension
        
        
        #self.standardThreshold = -log(self.motif.pBG**(self.wordlength+1) * (1.0/4)**self.wordlength)
        
    
    def run(self):
        db = self.index.db
        index = self.index
        print 'FAST HMM Mapping:    ',                  
        print '\t   [',                
        for indexWord,word in enumerate(index.allwords):
            if index.mappings[indexWord] != []:
                
                scoreWord = self.calcScoreWord(word)
                #print scoreWord                
                if scoreWord[0]!=inf and scoreWord[0]/scoreWord[1] < self.cutoffWords:
                                        
                    for indexDB in index.mappings[indexWord]:
                        
                        scoreExtension = self.calcScoreExtension(indexDB,scoreWord)
                        
                        if scoreExtension[0]!=inf and scoreExtension[0]/scoreExtension[1] < self.cutoffExtension:
                            
                            self.hits.append([indexDB, scoreExtension[0],scoreExtension[1]])
                                            
                                            
            #if i%(index.length/10) == 0: print '%3d%%' % (math.ceil(100.0*i/index.length))
            if indexWord%(index.length/10) == 0: print '|',                                                
        print ']'
        
        return        
        
    def calcScoreWord(self,word):
        """ The threshold is defined as the probablity that the hmm stayed in the background state """
        
        index = self.index
        
        windows = self.motiflength - self.wordlength + 1
        if windows < 1 : raise("Length of index word cannot be longer than the hmm")
        
        scorePWM = -log(self.motif.pPWM)
        scoreBG = -log(self.motif.pBG)
        
        for i,letter in enumerate(word):
            j = int(letter)
            scorePWM += -(log(self.motif.pwm[j,i]))
            scoreBG += -log(self.motif.background[j]*self.motif.pBG)
            if scorePWM == inf: break
 
        return [scorePWM,scoreBG]

    def calcScoreExtension(self,indexDB,score):        

        scorePWM = score[0]
        scoreBG = score[1]
        
        index = self.index        
        db = self.index.db


        for r in range(self.wordlength,self.motiflength):                        
            letter = int(db.db[indexDB+r])
            scorePWM += -(log(self.motif.pwm[letter,r]))
            scoreBG += -log(self.motif.background[letter]*self.motif.pBG)
            if scorePWM == inf: break

        return [scorePWM , scoreBG]    
    
    
    
class ForwardBackward:
    def __init__(self,motif,db):
        #if db.__class__ != DB: raise("Instance of DB class needed for index")        
        self.db = db
        self.obs = db.db
        self.ns = motif.length+1
            
        self.pE = concatenate((motif.background.reshape(-1,1),motif.pwm),axis=1)
        self.pE = self.pE.transpose()
        
        self.pT = zeros((self.ns,self.ns),float)        
        self.pT[0,0] = motif.pBG
        self.pT[0,1] = motif.pPWM
        self.pT[self.ns-1,0] = 1
        for i in range(1,self.ns-1):
            self.pT[i,i+1] = 1
        
        self.pT = self.pT.transpose()    
        
            
        #set first alpha , we assume we are first in the background state
        self.alpha =  zeros((self.ns,1),float)        
        self.alpha[0] = 1 * self.pE[0,0]        
        self.alphas = []
        self.alphas.append(self.alpha)
        
        self.beta = ones((self.ns,1),float)
        self.betas = []
        self.betas.append(self.beta)
        
        self.pSobs = []
       
    def run(self):
        print 'Forward Backward Algorithm:'
        self.forwardbackward()
        
        self.pObs = sum(self.alphas[-1])
        
        print 'Reconstruct :    ',                                            
        print '\t   [',                        
        l = len(self.obs) 
        for i in range(0,l):
            
            self.pSobs.append(self.alphas[i]*self.betas[-(i+1)]/self.pObs)
#            for t=1:200;
#                PsObs = [PsObs alphas(:,t).*betas(:,t)./Pobs];
#            end
            if i%(l/10) == 0: print '|',                                                
        print ']'         
             
        
    def forwardbackward(self):
        l = len(self.obs)
        suma = sum(self.alpha)
        
        print 'Alpha / Beta:    ',                                    
        print '\t   [', 
        for i in range(1,l):
            j = l-1-i
            #print suma
            self.alpha = self.pE[:,int(self.obs[i])].reshape(-1,1) *  dot(self.pT , self.alpha )
            
            self.beta = dot( self.pT.transpose(),  self.pE[:,int(self.obs[j+1])].reshape(-1,1) * self.beta)
            
            self.alpha = self.alpha / suma
            self.beta = self.beta / suma
            
            
            self.alphas.append(self.alpha)        
            self.betas.append(self.beta)
            
            suma = sum(self.alpha)
             
            if i%(l/10) == 0: print '|',                                                
        print ']'   
        
        
        
    


    
class Viterbi:
    
    def __init__(self,motif,db):
        #if db.__class__ != DB: raise("Instance of DB class needed for index")        
        self.db = db
        self.obs = db.db
        self.ns = motif.length+1
            
        self.pE = concatenate((motif.background.reshape(-1,1),motif.pwm),axis=1)
        self.pE = self.pE.transpose()   
        self.pElog = -numpy.log(self.pE) 
        
        self.pT = zeros((self.ns,self.ns),float)        
        self.pT[0,0] = motif.pBG
        self.pT[0,1] = motif.pPWM
        self.pT[self.ns-1,0] = 1
        for i in range(1,self.ns-1):
            self.pT[i,i+1] = 1
            
        self.pT = self.pT.transpose()
        self.pTlog = -numpy.log(self.pT)
        
        self.V = []
        
        self.nextV =  zeros((self.ns),float)
        #self.alpha[:] = 1e-300
        self.nextV[0] = 1 + self.pElog[0,0]    
        #self.alpha =  -numpy.log(self.alpha)        
        self.V.append(self.nextV)
        
    def run(self):
        
        for i in range(1,len(self.obs)):
            
            last = self.nextV
            
            self.pT
            
            
            
            #self.alpha = self.pElog[:,int(self.obs[i])] +  self.pT.sum(0) + self.alpha.sum()
                
#
#    Fast HMM Mapping
#    Fabian , dec08
#
#    hmm.py : A class representing the motif HMM that we are looking for in the db
#
#

#from pylab import *
#import numpy 
#import scipy 
import os

from hmm import *

#alphabet = [ 0, 1, 2, 3]
#alphabet = ['A','C','G','T']

class DB:
   
    def __init__(self,alphabet,name):
        self.alphabet = alphabet
        self.db = ''
        self.name = name
        self.length = 0
        
        if os.path.exists('db/'+self.name+'.txt'):
            print 'Loading Database %s' % self.name
            self.loadDB()
        else:
            print 'Database %s has to be created.' % self.name
            print '-> Call db.createDB(motif,length)' 
        

    def createDB(self,motif,length=1e2):
        """
        Method to create DB from given motif
        """
        #if motif.__class__ != Motif: raise("Motif needed to create DB")
        self.length = int(length)
        self.db=''
        count = 0
        while count < self.length:
            if random() < motif.pPWM:                
                #choose sequence drawn from pwm
                for position in motif.pwm.T:
                    #print position
                    r = random()
                    t = 0.0
                    for i,letter in enumerate(position):
                        if r < letter+t:
                            self.db += str(self.alphabet[i])
                            count +=1
                            break
                        else:
                            t += letter                          
                #increase count
                pass
            else:
                #choose nucleotide from background dist.
                r = random()
                t = 0.0
                for i,letter in enumerate(motif.background):
                    if r < letter+t:
                        self.db += str(self.alphabet[i])
                        count +=1
                        break
                    else:
                        t += letter  
        self.saveDB()                                    
        
    def saveDB(self):                
        file = open('db/'+self.name+'.txt','w') 
        file.write(self.db)
        file.close()
    def loadDB(self):                
        file = open('db/'+self.name+'.txt','r') 
        self.db = file.readline()
        self.length = len(self.db)
        file.close()  
        
    def stats(self):
        print 'Statistics for '+self.name
        count = [0 , 0 , 0, 0]
        
        for letter in self.db:
            #for i in range(0,len(self.alphabet)):
            #   if(letter==self.alphabet[i]):
            count[int(letter)]+=1
                    #break
        for i in range(0,len(self.alphabet)):
            print '\t%s : %s (%f%%)' %(("A", "C", "G", "T")[self.alphabet[i]],count[i],(100.0*count[i]/self.length))
            
    def __str__(self):
        return self.db      
    
   
class Index:
    def __init__(self,db,wordlength=8):
        
        if db.__class__ != DB: raise("Instance of DB class needed for index")
        
        self.db = db
        self.alphabet = self.db.alphabet
        self.name = db.name+'.Index'
        self.wordlength =  wordlength
        self.length = len(self.alphabet)**self.wordlength
        
        self.allwords = []
        self.mappings = []
        
        if os.path.exists('index/allwords'+str(self.wordlength)):
            print 'Loading List of all words %s' % self.name
            self.loadWords()
        else:
            print 'Database %s has to be created.' % self.name             
            self.createAllSeq(self.alphabet,[],self.wordlength)
            self.saveWords()

        if os.path.exists('index/'+self.name+'.'+str(self.wordlength)):
            print 'Loading Index %s' % self.name
            self.loadIndex()
        else:
            print 'Index %s has to be created.' % self.name
            self.createMappings()
            self.saveIndex()
        
    #def create(self):
        #create all possible words for index
        #self.createAllSeq(self.alphabet,[],self.wordlength)
        #make the mapping
        #self.createMappings()
        #return
    
    def createAllSeq(self,alphabet,cs=[],length=6):
        if length==0:
            rstr = ''
            for e in cs:
                rstr += str(e)
            self.allwords.append(rstr)            
            pass
        else:
            a = alphabet[:]                
            for i,letter in enumerate(a):
                ct = cs[:]                 
                ct.append(letter)                
                anew = a[:]
                anew.remove(letter)                
                self.createAllSeq(a, ct, length-1)      
                ct = cs         
                    
    def createMappings(self):
        #for wordlength 8
        #weights = [16384,4096,1024,256,64,24,4,0] 
        weights = []
        pw = 0
        for i in range(self.wordlength-1,-1,-1):
            weights.append(2**(2*i))
                       
        for i in range(0,len(self.alphabet)**self.wordlength):            
            self.mappings.append([])
        #self.mappings = [[]] * self.length 
        
        print 'Creating Index Mapping (each | is 10%)'    
        for i in range(0,(self.db.length-self.wordlength+1)):
            cI = 0
            seq = self.db.db[i:(i+self.wordlength)]
            for j,l in enumerate(seq):
                cI += int(l) * weights[j]
            self.mappings[cI].append(i)
            if i%(self.db.length/10) == 0: print '|',
        print 
                   
    def saveWords(self):
        file = open('index/allwords'+str(self.wordlength),'w') 
        for word in self.allwords:
            file.write(word+'\n')
        file.close()
        
    def saveIndex(self):    
        file = open('index/'+self.name+'.'+str(self.wordlength),'w')
        for mapping in self.mappings:
            if mapping == []:
                file.write('none\n')
            else:
                for m in mapping:
                    file.write(str(m)+' ')
                file.write('\n')
        file.close()
        
    def loadWords(self):
        self.allwords = []
        self.mappings = []

        file = open('index/allwords'+str(self.wordlength),'r') 

        for line in file.readlines(): 
            self.allwords.append(line[:-1])      
        file.close()
        
    def loadIndex(self): 
        file = open('index/'+self.name+'.'+str(self.wordlength),'r')
        for line in file.readlines():            
            if line[0]!='n': 
                t = []
                for n in line[:-1].split():
                    t.append(int(n))
                self.mappings.append(t)   
            else: self.mappings.append([])
            
        file.close()
        
                    
                               
    
    

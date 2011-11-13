
realalphabet = ['A','C','G','T']
def printWord(word):
    s=''
    for w in word:
        s+=realalphabet[int(w)]
    return s
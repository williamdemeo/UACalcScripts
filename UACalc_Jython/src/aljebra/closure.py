class Closure(object):
    
    
    def sd_embedding(self):
        '''Return an array with (i,j) entry equal to the index of the block of pars[j] containing i.'''
        answer = []
        n = pars[0].universeSize()
        for i in range(n):
            temp = []
            for j in range(len(pars)):
                temp.append(pars[j].blockIndex(i))
            answer.append(temp)
            return answer

    
    




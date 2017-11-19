'''
RNA structure
a single RNA structure
'''
class RNAstructure:

    def __init__(self):
        self.sequence = ''
        self.pair = []   # base number of the paired base

    def CTRead(self,filename):
        '''
        Read in RNAstructure CT file
        format example:
          240  ENERGY = -49.3   mr_s129
            1 A       0    2    0    1
            2 A       1    3    0    2
            3 C       2    4    0    3
            4 C       3    5    0    4
            5 A       4    6    0    5

            base_number base previous_base next_base paired_base base_number

        usage
            rna.CTRead(filename)
        '''
        with open(filename, 'r') as ct:
            line = ct.readline()
            print('firstline:', line)                
            #print('field:',field)

            if line.find('ENERGY')>=0:
                # 
                pass
            else:
                # probknot file
                field = line.split()
                self.length = int(field[0])
                self.id = field[1]
            
            self.pair = [0] * (self.length+1)
            nline = 0
            nbase = 0
            for line in ct:
                n,base,prev,next,pair,n2 = line.split()
                print( 'n:',n,'base:',base,'pref:',prev,'next:',next,'pair:',pair,'n2:',n2)
                self.sequence += base
                if pair != '0':
                    self.pair[int(pair)] = int(n)
                    self.pair[int(n)] = int(pair)
                nbase += 1

    def __str__(self):
        str = 'RNA:\n'
        for key in self.__dict__:
            str += '{0} = {1}\n'.format(key, self.__dict__[key])

        return str

    def getStemList(self):
        unpaired = 2
        instem = False
        lpos = 0
        rpos = 0
        for pos in range(0, len(self.pair)-1):
            if self.pair[pos] and self.pair[pos] < pos:
                # check that left stem is smaller than right stem
                continue

            term = pos-lpos > unpaired or rpos-self.pair[lpos] > unpaired
            #print( 'pos:',pos, 'pair:',self.pair[pos], 'instem:', instem, 'term:', term, 'lpos:',lpos, 'rpos:', rpos)
            if instem:
                # currently in a stem
                if term:
                    # end old stem
                    print( 'stem:',lstart,lpos,self.pair[lpos],rstop )
                    lstart = pos
                    rstop = self.pair[pos]
                    instem = False
               
                if self.pair[pos]:
                        lpos = pos
                        rpos = self.pair[pos]
                        instem = True
                        
            else:
                # not in a stem
                if self.pair[pos]:
                    # start a new stem
                    lstart = pos
                    rstop = self.pair[pos]
                    lpos = pos
                    rpos = self.pair[pos]
                    instem = True

class Stem:
    def __init__(self):
        self.lbegin = 0
        self.lend   = 0
        self.rbegin = 0
        self.rend   = 0
        self.lvienna = ''
        self.rvienna = ''
        
if __name__ == '__main__':
    rna = RNAstructure()
    rna.CTRead('data/mr_s129.probknot.ct')
    print(rna )
    rna.getStemList()
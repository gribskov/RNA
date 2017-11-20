'''
RNA structure
a single RNA structure
'''
class RNAstructure:

    def __init__(self):
        self.sequence = ''
        self.pair = []   # base number of the paired base
        self.stemlist = []

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

    def stemListGet(self):
        stem = Stem()
        unpaired = 2
        instem = False
        for pos in range(0, len(self.pair)-1):
            if self.pair[pos] and self.pair[pos] < pos:
                # check that left stem is smaller than right stem
                continue

            term = pos-stem.lend > unpaired or stem.rend-self.pair[stem.lend] > unpaired
            #print( 'pos:',pos, 'pair:',self.pair[pos], 'instem:', instem, 'term:', term, 'lpos:',lpos, 'rpos:', rpos)
            if instem:
                # currently in a stem
                if term:
                    # end old stem
                    self.stemlist.append(stem)
                    stem = Stem()
                    stem.lbegin = pos
                    stem.rend = self.pair[pos]
                    instem = False
               
                if self.pair[pos]:
                        stem.lend = pos
                        stem.rbegin = self.pair[pos]
                        instem = True
                        
            else:
                # not in a stem
                if self.pair[pos]:
                    # start a new stem
                    stem = Stem()
                    stem.lbegin = pos
                    stem.rend = self.pair[pos]
                    stem.lend = pos
                    stem.rbegin = self.pair[pos]
                    instem = True

        if instem:
            # if you  end in a stem, save it before closing
            self.stemlist.append(stem)

    def stemlistFormat(self):
        n = 0
        str = ''
        for stem in self.stemlist:
            n += 1
            str += '{0}\t{1}\n'.format(n, stem.formatted())
        return str
 
class Stem:
    def __init__(self):
        self.lbegin = 0
        self.lend   = 0
        self.rbegin = 0
        self.rend   = 0
        self.lvienna = ''
        self.rvienna = ''

    def formatted(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(self.lbegin,self.lend,self.rbegin,self.rend,self.lvienna,self.rvienna)

if __name__ == '__main__':
    rna = RNAstructure()
    rna.CTRead('data/mr_s129.probknot.ct')
    rna.stemListGet()
    print(rna )
    print('Stemlist\n')
    print( rna.stemlistFormat() )
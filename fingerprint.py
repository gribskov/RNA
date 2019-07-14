import sys
import json
import yaml

class Fingerprint(dict):
    """=============================================================================================
    A finger print is a spectrum of fixed size motifs identified in a structure.  The motifs are
    referenced to a motif dictionary to enable simple identification of parents.

    Sorting the results

    # sort by count
    for motif,count in sorted(fingerprint.motif.items(),key=lambda x:x[1],reverse=True):

    # sort alphabetically
    for motif,count in sorted(fingerprint.motif.items(),key=lambda x:x[0]):

    10 July 2019     Michael Gribskov
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        super().__init__(self)
        self.information = { 'Date': ''
                            }
        self.motif = {}
        self.nmotif = 0
        self.total = 0      # total count of added motifs


    def add(self, string, count=1):
        """-----------------------------------------------------------------------------------------
        Adss a motif to the fingerprint, creating a new entry if needed.  Agnostic as to the
        encoding of the motif string

        :param string: str, encoded motifs string
        :return: int, number of motifs in fingerprint
        -----------------------------------------------------------------------------------------"""
        self.total += count
        if string in self.motif:
            self.motif[string] += count
        else:
            self.motif[string] = count
            self.nmotif += count

        return self.nmotif
    
    def mincount(self):
        """-----------------------------------------------------------------------------------------
        Return the count of the motif with the smallest count

        :return: int
        -----------------------------------------------------------------------------------------"""
        return min(self.motif.values(), default=0)

    def minkey(self):
        """-----------------------------------------------------------------------------------------
        Return the motif with the lowest count

        :return: str
        -----------------------------------------------------------------------------------------"""
        m = self.motif
        m = min(m.items(), default=0, key= lambda x:x[1])
        # print('smallest', m)

        return m[0]

    def toJSON(self):
        """-----------------------------------------------------------------------------------------
        Convert fingerprint to JSON string

        :return: str
        -----------------------------------------------------------------------------------------"""
        fields = ['information', 'total', 'nmotif', 'motif']

        root = []
        for i in range(len(fields)):
            # look up the values of each field and add to dispatch list
            root.append(getattr(self, fields[i]))

        return json.dumps({fields[i]: root[i] for i in range(len(fields))}, indent=4)

    def toYAML(self):
        """-----------------------------------------------------------------------------------------
        Convert fingerprint to JSON string

        :return: str
        -----------------------------------------------------------------------------------------"""
        fields = ['information', 'total', 'nmotif', 'motif']

        root = [ { 'fingerprint': [ {'information':self.information},
                                    {'total':self.total},
                                    {'nmotif':self.nmotif},
                                    {'motif':self.motif} ] } ]

        return yaml.dump(root, indent=2)

    def writeYAML(self, file):
        """-----------------------------------------------------------------------------------------
        Write the fingerprint to a file in YAML format

        :param file: fp/str, either a string or an open file
        :return: True
        -----------------------------------------------------------------------------------------"""
        if isinstance(file, str):
            # file argument is string, try to open
            try:
                fp = open(file, 'w')
            except OSError:
                sys.stderr.write('fingerprint.writeYAML - error opening file ({})'.format(file))
                exit(1)
        else:
                # file is not str, assume it is a file pointer
                fp = file

        fp.write(self.toYAML())
        # fp.close()    if fp is stdout, not good

        return True


    def readYAML(self, file):
        """-----------------------------------------------------------------------------------------
        read the fingerprint from a file in YAML format

        :param file: fp/str, either a string or an open file
        :return: int, number of motfs
        -----------------------------------------------------------------------------------------"""
        if isinstance(file, str):
            # file argument is string, try to open
            try:
                fp = open(file, 'r')
            except OSError:
                sys.stderr.write('fingerprint.readYAML - error opening file ({})'.format(file))
                exit(1)
        else:
            # file is not str, assume it is a file pointer
            fp = file

        f = yaml.load(fp)
        root = f[0]['fingerprint']

        # fields = ['information', 'total', 'nmotif', 'motif']

        self.information = root[0]['information']
        self.total = root[1]['total']
        self.nmotif = root[2]['nmotif']
        self.motif = root[3]['motif']


        return nmotif


    def add_parents(self, motifdb):
        """-----------------------------------------------------------------------------------------
        Look up the parents in the provided motif database and add them to the fingerprint.  The
        count of each parent is incremented for each of its children.  for nested children this may
        not make perfect sense

        :param motifdb:
        :return: int, total number of motifs
        -----------------------------------------------------------------------------------------"""
        children = list(self.motif.keys())
        for child in children:
            for parent in motifdb.parent[child]:
                # print('child {}\t\tparent {} {}'.format(child, parent, self.motif[child]))
                self.add(parent, count=self.motif[child])

        return self.nmotif

# ##################################################################################################
# Testing
# ##################################################################################################
if __name__ == '__main__':

    finger = Fingerprint()
    finger.append('dummy')
    finger.append('dummier')
    finger.append('dummiest')
    print(finger[0])
    print(finger)


    exit(0)

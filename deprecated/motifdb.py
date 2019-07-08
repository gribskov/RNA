"""=================================================================================================
DEPRECATED - merged into Xios
====================================================================================================
MotifDB class is for creating and using XIOS graph dictionaries

011220

0,5  1,2  3,4

00 -> 11 *
    0011 -> 1122
        010122
        011022
        011202 *
        011220 *
    0101 -> 1212  *
        010212 *
        012012 *
        012102 *
        012120 *
    0110 -> 1221  *
        010221 *
        012021 *
        012201 *
        012210 *

01
    01 23

    02 13
    03 12


Michael Gribskov     15 June 2019
================================================================================================="""
import json
import datetime
from xios import Xios


class MotifDB():
    """=============================================================================================
    I keep going back and forth over whether the object should be a dict, or whether it should have
    defined fields.  Defined fields requires fewere accessors, but a general accessor with getattr
    could be used.
        motifdb.fields = list of fields in order for serialization
        motifdb.information = metadata, standard info has accessors, but any notes can be added
        motifdb.db[nstem] = [dfshex, dfshex, dfshex ...]
        motifdb.lenidx = list, index is number of stems
    ============================================================================================="""

    def __init__(self):
        self.fields = ['information', 'n', 'db', 'lenidx', 'parent']
        self.information = {}  # for metadata
        self.n = 0
        self.db = []
        self.parent = {}
        self.lenidx = []  # lists of motifs indexed by number of stems (motif length)

    def add_with_len(self, motif):
        """-----------------------------------------------------------------------------------------
        Add a single motif to the database

        :param motif:
        :return: int, number of motifs
        -----------------------------------------------------------------------------------------"""
        db = self.db
        n = len(db)
        nstem = len(motif)//2
        if nstem > len(self.lenidx):
            for i in range(len(self.lenidx), nstem):
                self.lenidx.append([])

        db.append(motif)
        self.lenidx[nstem - 1].append(db[n])
        self.n = len(db)

        return self.n

    def add_parent(self, child, parent):
        """-----------------------------------------------------------------------------------------
        add parent to the parent list of child, and add all the parents of parent to the parent
        list of child

        :param parent:
        :return:
        -----------------------------------------------------------------------------------------"""
        # if child not in self.parent:
        #     self.parent[child] = []

        if parent not in self.parent[child]:
            self.parent[child].append(parent)

        for p in self.parent[parent]:
            if p not in self.parent[child]:
                self.parent[child].append(p)

        return len(self.parent)


    def setdate(self):
        """-----------------------------------------------------------------------------------------
        set date in information
        :return: datetime
        -----------------------------------------------------------------------------------------"""
        daytime = datetime.datetime.now()
        self.information['date'] = daytime.strftime('%Y-%m-%d %H:%M:%S')

        return daytime

    def setname(self, string):
        """-----------------------------------------------------------------------------------------
        set database name in information

        :return: int, length of new name field
        -----------------------------------------------------------------------------------------"""
        self.information['name'] = string

        return len(string)

    def setsource(self, source):
        """-----------------------------------------------------------------------------------------
        Name of the program creating the database
        :return:
        -----------------------------------------------------------------------------------------"""
        self.information['source'] = source

        return len(source)

    def toJSON(self):
        """-----------------------------------------------------------------------------------------
        Convert database to JSON string

        :return: str
        -----------------------------------------------------------------------------------------"""
        fields = self.fields
        dispatch = []
        for i in range(len(fields)):
            dispatch.append(getattr(self, fields[i]))

        return json.dumps({fields[i]: dispatch[i] for i in range(len(fields))},indent=4)

    def toFile(self,fp):
        """-----------------------------------------------------------------------------------------
        Write a formatted version of the motif database to a file. use sys.stdout to as the file if
        you want it on STDOUT

        :param fp:
        :return:
        -----------------------------------------------------------------------------------------"""
        fields = ['information', 'n', 'db' ]

        for i in range(len(fields)):
            data = getattr(self, fields[i])
            if fields[i] == 'information':
                for tag in data:
                    fp.write('{}\t{}\n'.format(tag, data[tag]))
            elif fields[i] == 'n':
                fp.write('{}\t{}\n'.format('motifs', data))
            elif fields[i] == 'db':
                for m in data:
                    x = Xios()
                    x.ascii_decode(m)
                    fp.write('\t{}\n'.format(x))
                    for p in self.parent[m]:
                        x = Xios()
                        x.ascii_decode(p)
                        fp.write('\t\t{}\n'.format(x))

        return

# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # motifdb
    motif = MotifDB()
    motif.setdate()
    motif.setname('test database')

    motif.add_with_len('aabbcc', 1)
    motif.add_with_len('bbccdd', 2)
    motif.add_with_len('ccddee', 2)
    print(motif.toJSON())

    exit(1)

    # SerialRNA
    rnas = [[0, 0, 1, 1, 2, 2],
            [0, 1, 0, 1, 2, 2],
            [0, 1, 1, 2, 2, 0],
            [0, 1, 2, 1, 2, 0],
            [0, 0], []
            ]

    print('canonical form')
    noncanonical = [[3, 3, 0, 0, 1, 1], [1, 1, 2, 2, 3, 3], [1, 1, 2, 2, 4, 4], [3, 2, 0, 2, 0, 3]]
    for testcase in noncanonical:
        rna = SerialRNA(testcase)
        print('RNA {}'.format(rna))
        rna.canonical()
        print('\tcanonical {}'.format(rna))

    print('\nConnected components')
    for testcase in rnas:
        rna = SerialRNA(testcase)
        print('RNA {}'.format(rna))
        connected = rna.connected()
        if len(connected) > 1:
            for i in range(len(connected)):
                print('\tcomponent {}: {}'.format(i, connected[i]))

    print('\nExtension')
    for testcase in rnas:
        rna = SerialRNA(testcase)
        print('RNA {}'.format(rna))
        for new in rna.addstemzero():
            print('\t{} {}'.format(new, len(new.connected())))

    exit(0)

"""=====================================================================================================================
command.py

test version to improve command processing

Michael Gribskov 2/16/2026
====================================================================================================================="""
import yaml

class Command:
    """#############################################################################################
    convert the yaml file describing the workflow to a list of executable commands
    #############################################################################################"""

    def __init__(self, filename='workflow1.yaml'):
        """-----------------------------------------------------------------------------------------
        file        path to yaml file containing the workflow
        plan        python structure version of the plan
        stage       current stage
        def_main    dict, global definitions
        def_stage   dict, definitions for current stage
        -----------------------------------------------------------------------------------------"""
        self.file = filename
        self.parsed = None
        self.stage = ''
        self.stage_dir = ''
        self.def_main = {}
        self.def_stage = {}
        self.command = ''
        self.mult = {}
        self.late = {}
        self.defs = []

    def read(self):
        """-----------------------------------------------------------------------------------------
        read the workflow description from the yaml file. File is stored in Command.file

        :return: dict               workflow as a python dictionary
        -----------------------------------------------------------------------------------------"""
        fp = open(self.file, 'r')
        self.parsed = yaml.load(fp, Loader=yaml.FullLoader)
        fp.close()
        return self.parsed

    def main(self):
        """-----------------------------------------------------------------------------------------
        Set up the main definitions for the plan
        :return:
        -----------------------------------------------------------------------------------------"""
        pass

        def expand(self):
        """-----------------------------------------------------------------------------------------
        recursively substitute symbols in defs with their values. Definitions to expand are in
        self.parsed['definitions']

        :return: dict       definition dict after substitution
        -----------------------------------------------------------------------------------------"""
        stack = []
        symbol = {}
        defs = self.parsed['definitions']
        for d in defs:
            # add $ to the definition keys and add definitions that have $ on the stack,
            # these need expansion
            symbol['$'+d] = defs[d]
            if defs[d].find('$') > -1:
                stack.append(d)

        while stack:
            d = stack.pop()
            for s in symbol:
                if defs[d].find(s) == -1:
                    # symbol not present, skip to next symbol
                    continue

                # symbol is in this definition, expand it
                defs[d] = defs[d].replace(s, symbol[s])
                if defs[d].find('$') == -1:
                    # no additional $ in this definition, done with this definition, back to stack
                    break

                # there are still $, push definition back on the stack, and start from beginning of
                # symbol list
                stack.append(d)
                break

        return defs

    def command_generate(self):
        """-----------------------------------------------------------------------------------------
        generate a list of commands from the workflow for a specific stage
        1. split the command into tokens
        2. for each token, match to the keywords in plan['stage'][stage] and replace
           $keyword with the value from the yaml.
        3. for each token, match to keywords in plan['definitions'] and replace $keyword with the
           value from th yaml
        4. options in <> are processed using plan['stage'][stage]['rule'] at runtime

        break command lines from workflow yaml into
        commands
        mult - symbols generated with wild cards
        late - dependent symbols created from other symbols. they are late because they geneerally
                can't be resolved unitl the previous stage finishes

        :return: string                 command string
        -----------------------------------------------------------------------------------------"""
        current = self.parsed['stage'][stage]

        # merge stage definitions with global definitions and expand the command
        # separate the definitions into the command, simple symbols, mutltiple processing
        # symbols, and late symbols
        self.def_stage = {k: self.def_main[k] for k in self.def_main}
        for item in current:
            if item == 'command':
                self.command = current[item]
            elif current[item].find('*') > -1:
                self.mult[item] = current[item]
            elif current[item].find(':') == 0:
                self.late[item] = current[item]
            else:
                self.def_stage[item] = current[item]

        self.def_stage = self.expand(self.def_stage)
        for d in self.def_stage:
            self.command = self.command.replace(f'${d}', self.def_stage[d])

        for m in self.mult:
            for d in self.def_stage:
                self.mult[m] = self.mult[m].replace(f'${d}', self.def_stage[d])
        for litem in self.late:
            for d in self.def_stage:
                self.late[litem] = self.late[litem].replace(f'${d}', self.def_stage[d])

        return self.multiple()

    def multiple(self):
        """-----------------------------------------------------------------------------------------
        based on the entries in self.mult, that is, symbols that contain wildcards, generate
        lists of matching file names to be used in individual commands

        :return:
        -----------------------------------------------------------------------------------------"""
        filelist = {}
        for symbol in self.mult:
            filelist[symbol] = glob.glob(self.mult[symbol])
            # sys.stderr.write(f'filelist:{filelist}\tsymbol[{symbol}]:{self.mult[symbol]}\n')
            if not filelist[symbol]:
                # symbol expansion  produced an empty list
                # self.log.add('stage', f'Error: expansion of symbol ${symbol} produced and empty file list')
                sys.stderr.write(f'Error: expansion of symbol ${symbol} produced and empty file list\n')

        commandlist = []
        expand = self.command
        # sys.stderr.write(f'command:{self.command}\n')
        for m in self.multiple_gen(filelist):
            # m is a dict with a value for every multiple symbol
            # sys.stderr.write(f'm:{m}\n')
            expand = self.command
            for symbol in m:
                expand = expand.replace(f'${symbol}', m[symbol])
                for litem in self.late:
                    lcom = self.late[litem][1:]
                    for symbol in m:
                        # TODO not sure why the for symbol in m loop is doubly nested, it think its so that symbols get
                        # expanded both before and after a late function runs. may not be necessary but currently works
                        # so i'm leaving as is
                        file = os.path.basename(m[symbol])
                        # sys.stderr.write(f'base:{file}\tsymbol:{symbol}\tlcom:{lcom}\n')
                        lcom = lcom.replace(f'%{symbol}', f'{file}"')
                        lcom = f'"{lcom}'
                        # sys.stderr.write(f'lcom:{lcom}\n')
                        t = eval(lcom)
                        # sys.stderr.write(f't:{t}\tl:{l}\n\n')
                        expand = expand.replace(f'${litem}', t)

            commandlist.append(expand)

        if not commandlist:
            # in case there were no symbols in command
            commandlist.append(expand)

        return commandlist

    def multiple_gen(self, filelist):
        """-----------------------------------------------------------------------------------------
        generator to create all combinations of multiple values

        :return: dict   value for each entry in self.mult
        -----------------------------------------------------------------------------------------"""
        pos = {symbol: 0 for symbol in filelist}
        while True:
            combo = {symbol: filelist[symbol][pos[symbol]] for symbol in filelist}
            yield combo
            carry = 1
            for symbol in pos:
                pos[symbol] += carry
                carry = 0
                if pos[symbol] >= len(filelist[symbol]):
                    carry = 1
                    pos[symbol] = 0

            if carry:
                break

        return


####################################################################################################
# end of class Command
####################################################################################################
# ======================================================================================================================
# testing
# ======================================================================================================================
if __name__ == '__main__':
    workflow_file = 'workflow1.yaml'
    cmd = Command(filename=workflow_file)
    cmd.read()
    # cmd.parsed['definitions']
    cmd.def_main = cmd.expand()
    cmd.command_generate()

    exit(0)

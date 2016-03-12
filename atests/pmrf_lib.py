import subprocess
import os

PMRF_EXEC = '%s/../pmrf'%os.path.dirname(__file__)

class pmrf_lib(object):

    def __init__(self):
        self._output = ''
        pass

    def run_pmrf(self, opts):
        cmd = '%s %s'%(PMRF_EXEC, opts)
        self._output = subprocess.check_output(cmd.split())
        return self._output

    def read_file(self, filename):
        return open(filename).read()

    def should_exist(self, filename):
        if not os.path.exists(filename):
            raise AssertionError("'%s' should have existed"%filename)

    def remove_file(self, filename):
        os.remove(filename)

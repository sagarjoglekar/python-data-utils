import unittest
from figs_base import *
from tests import *
import os

class testfigsbase(unittest.TestCase):

    def test_read_write_line(self):
        write_line("test.csv",[(1,2),(3,4),(5,6)],xlbl="X",ylbl="Y",comment=None)
        if os.path.exists("test.csv"):
            self.assertTrue(True)
            x,y = read_line("test.csv")
            if x[0] == '1' and y[0] == '2':
                self.assertTrue(True)
            else:
                print "File readLine failed " + "1 not " + str(x) + "2 not " + str(y)
                self.assertTrue(False)

        else:
            self.assertTrue(False)




if __name__ == '__main__':
    unittest.main()

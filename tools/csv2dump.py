import sys
import lammpsUtils as LU

HELP_MESSAGE="""Program usage:

    python csv2dump.py -i csvfilename.csv -o dumpfilename
    
"""

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", help="Input file CSV file with atoms to be save into LAMMPS dump format.")
parser.add_option("-o", "--output", action="store", type="string", help="Name of output LAMMPS dump file")

(options, args) = parser.parse_args()
opt_dict = vars(options)
if options.input==None:
    sys.exit("ERROR: Please, specify input CSV file.")
if options.output==None:
    sys.exit("ERROR: Please, specify output LAMMPS dump file name.")
opt_dict = vars(options)

#if options.F and options.coordination:
LU.csv2dump(options.input,options.output)

#===============================================#
# BWT Implementation                            #
# Michael Ting					#
# 10 October 2013				#
#===============================================#

The implementation of the Burrows-Wheeler Transform
can be found in bwt.py.

The script was written in Python and tested and run 
with Python2.7.

To utilize the script, run:

BWT encoding:

  $ python bwt.py -bwt seq.fasta seq.bwt.fasta

where seq.fasta is your input sequence in FASTA format, and
seq.bwt.fasta is your output sequence file name, or

BWT Decoding:

  $ python bwt.py -ibwt seq.bwt.fasta seq.ibwt.fasta

where seq.bwt.fasta is your input BWT in FASTA format, and
seq.ibwt.fasta is your output sequence file name.

The input to bwt.py assumes your FASTA sequence has
a termination symbol such as "$" in "hello$".

#===========#
# Test Code #
#===========#

I also wrote a script to test functions in bwt.py, called 
test_bwt.py. This script is run using the nose 1.3.0 package
for Python2.7.

To run the tests, install nose 1.3.0 for Python2.7 and run:

$ nosetests test_bwt.py -v

Which should indicate whether or not the functions in
bwt.py pass the test code in test_bwt.py


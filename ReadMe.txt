The convert_nef_nmrpipe.py script converts a nef formatted list of protien assignments to the NMRPipe Table Format, for example for use 
in TALOS. First, copy the assignments section of your nef file into a txt file and give this a name. See CESA8_Assignment.txt for an 
example. Second, in convert_nef_nmrpipe.py, update the PATH and Sequence values to match your protein. Also, in #Read and extract data
adjust the number of NMRpipe_table += inputs so that all of your Sequence1 - SequenceX lines are added in order.

This script works for protein backbone assignments. Some manual adjustments may be required with the output table format. Compare the 
output to an example table format file that has worked in TALOS before, for example CesA8_fullassignment_NMRPipe_table_format.tab.
"""
The script converts .nef format text files into NMRPipe Table Format. The nef format is:

         A      9-1     MET    C   173.9428312  .             C  13  1  false  .

Some manual editing may be  required after. See a previously uploaded file and match the corrections. Also, for Gly, rename protons to HA2 and HA3 if required



"""

#Import

#Constant

PATH = 'Paxx_assignment.txt'
Sequence1 = 'DATA SEQUENCE MDPLSPPLCT LPPGPEPPRF VCYCEGEESG EGDRGGFNLY VTDAAELWST\n' #NMRPipe wants the sequence split up like this, in 50 AA sequences.
Sequence2 = 'DATA SEQUENCE CFTPDSLAAL KARFGLSAAE DITPRFRAAC EQQAVALTLQ EDRASLTLSG\n'
Sequence3 = 'DATA SEQUENCE GPSALAFDLS KVPGPEAAPR LRALTLGLAK RVWSLERRLA AAEETAVSPR\n'
Sequence4 = 'DATA SEQUENCE KSPRPAGPQL FLPDPDPQRG GPGPGVRRRC PGESLINPGF KSKKPAGGVD\n'
Sequence5 = 'DATA SEQUENCE FDET\n\n'
FIRST_RESID = 'DATA FIRST_RESID 1\n\n' #M is Met1 for example.
REMARK = 'REMARK Chemical Shift Table for Paxx\n\n'
filename = 'Paxx_assignment_NMRPipe_table_format.tab'
full_sequence = 'MDPLSPPLCTLPPGPEPPRFVCYCEGEESGEGDRGGFNLYVTDAAELWSTCFTPDSLAALKARFGLSAAEDITPRFRAACEQQAVALTLQEDRASLTLSGGPSALAFDLSKVPGPEAAPRLRALTLGLAKRVWSLERRLAAAEETAVSPRKSPRPAGPQLFLPDPDPQRGGPGPGVRRRCPGESLINPGFKSKKPAGGVDFDET'
pdb_shift = 0 #How much is the mismatch between my chain and the pdb chain indexes?

#Amino acid dictionary

AADict = dict()
single_letter = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
three_letters = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
for letter in single_letter:
    index = single_letter.index(letter)
    three_letter = three_letters[index]
    AADict[letter] = three_letter

#Load data

with open(PATH, mode = 'r') as datafile:
    nef_data = datafile.read()

#Read and extract data

NMRpipe_table = ''
NMRpipe_table += REMARK

NMRpipe_table += FIRST_RESID
NMRpipe_table += Sequence1
NMRpipe_table += Sequence2
NMRpipe_table += Sequence3
NMRpipe_table += Sequence4
NMRpipe_table += Sequence5
NMRpipe_table += 'VARS   RESID RESNAME ATOMNAME SHIFT\n'
NMRpipe_table += 'FORMAT %4d   %1s     %4s      %8.3f\n\n'

x = ['1','2','3','4','5','6','7','8','9']

nef_data_rows = nef_data.splitlines()
for row in nef_data_rows: #all of these operations, The script reads each row e.g. '         A      9-1     MET    C   173.9428312  .             C  13  1  false  .' from the left and extracts details.
    row = row.strip(' ')
    chain = row[0] # E.g. Chain A
    row = row.strip(chain)
    row = row.strip(' ')
    RESID = row.split(' ')[0]
    if '-1' in RESID: #Added this line to adjust 150-1 back to 149.
        RESID_plusone = RESID[:-2] #this removes the last two values, being the '-1' in the string. Its still plus one of the value. 150-1 is 150 here.
        row = row.strip(RESID)
        row = row.strip('-')
        row = row.strip('1')
        RESID = str((int(RESID_plusone)-1)) #Now 150-1 has been written back to 149.   
    row = row.strip(RESID)
    row = row.strip(' ')
    RESNAME = row.split(' ')[0]
    if RESNAME == '.': #In the 150-1 rows, there is no RESNAME. Need to determine the RESID and get RESNAME from the protein sequence.
        single_aa = full_sequence[int(RESID)-1] #-1 because in python, 0 refers to the first residue. So index of amino acid 1 is 0 in the string
        row = row.strip(RESNAME)
        RESID=(AADict[single_aa])
    row = row.strip(RESNAME)
    row = row.strip(' ')
    ATOMNAME = row.split(' ')[0]
    if ATOMNAME == 'H': #NMRPipe requires HN for amide proton
        ATOMNAME = 'HN'
    row = row.strip(ATOMNAME)
    row = row.strip(' ')
    SHIFT = row.split(' ')[0]
    if len(SHIFT) != 11: #Basically, all shifts are strings of len=11 e.g. '173.9428312'. If it is less, e.g. '29.4683351', then we add spaces to make it '  29.4683351'. TALOS wants it like this.
        SHIFT = ' '+SHIFT
    if len(SHIFT) != 11:
        SHIFT = ' '+SHIFT
    if len(SHIFT) != 11:
        SHIFT = ' '+SHIFT
    if len(SHIFT) != 11:
        SHIFT = ' '+SHIFT
    if len(SHIFT) != 11:
        SHIFT = ' '+SHIFT
    if len(SHIFT) != 11: #These ones contain an extra shift. Thus, want rid of first 5 spaces that we just added
        SHIFT = SHIFT[5:]
    if RESNAME != '.': #Quality control, don't want any unassigned peaks in here.
        resid = int(RESID)
        correctresid = resid + pdb_shift #This adjusts the resid in the output by the amount in pdb_shift, if required. Normally set pdb_shift to 0.
        RESID = str(correctresid)
        #All of the below are to ensure the string length of the text is 31. This can change based on if its HBx vs HA, or on RESID of 6 vs 158.
        if RESID in x and len(ATOMNAME) == 2: #CA, CB, NH, HA 
            line = '   '+RESID+' '+RESNAME+'     '+ATOMNAME+'    '+SHIFT+'\n'
            NMRpipe_table += line
        if RESID in x and len(ATOMNAME) == 1: #N
            line = '   '+RESID+' '+RESNAME+'     '+ATOMNAME+'     '+SHIFT+'\n'
            NMRpipe_table += line
        if RESID in x and len(ATOMNAME) == 3: #e.g. HBx, HBy
            line = '   '+RESID+' '+RESNAME+'     '+ATOMNAME+'   '+SHIFT+'\n'
            NMRpipe_table += line
        if RESID in x and len(ATOMNAME) == 4: #e.g. HD1x
            line = '   '+RESID+' '+RESNAME+'     '+ATOMNAME+'  '+SHIFT+'\n'
            NMRpipe_table += line
        if RESID not in x and len(ATOMNAME) == 2: 
            line = '  '+RESID+' '+RESNAME+'     '+ATOMNAME+'    '+SHIFT+'\n'
            if len(line) == 32: #If RESID is three digits, len line is 32. if RESID is four digits, len line is 33. # This trims some whitespace to make all lines 31 digits lone (including whitespace)
                line = line[1:]
            NMRpipe_table += line
        if RESID not in x and len(ATOMNAME) == 1:
            line = '  '+RESID+' '+RESNAME+'     '+ATOMNAME+'     '+SHIFT+'\n'
            if len(line) == 32: #If RESID is three digits, len line is 32. if RESID is four digits, len line is 33. # This trims some whitespace to make all lines 31 digits lone (including whitespace)
                line = line[1:]
            NMRpipe_table += line
        if RESID not in x and len(ATOMNAME) == 3: #e.g. HBx, HBy
            line = '  '+RESID+' '+RESNAME+'     '+ATOMNAME+'   '+SHIFT+'\n'
            if len(line) == 32: #If RESID is three digits, len line is 32. if RESID is four digits, len line is 33. # This trims some whitespace to make all lines 31 digits lone (including whitespace)
                line = line[1:]
            NMRpipe_table += line
        if RESID not in x and len(ATOMNAME) == 4: #e.g. HD1x
            line = '  '+RESID+' '+RESNAME+'     '+ATOMNAME+'  '+SHIFT+'\n'
            if len(line) == 32: #If RESID is three digits, len line is 32. if RESID is four digits, len line is 33. # This trims some whitespace to make all lines 31 digits lone (including whitespace)
                line = line[1:]
            NMRpipe_table += line

NMRpipe_table += '...'

#Save output

file = open(filename, 'w')
file.write(NMRpipe_table)
file.close()

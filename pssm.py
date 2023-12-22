# J Ludwig 2017
# My first python script

#/usr/bin/python

import sys #reading files
import re #matching the header line
import math #log
a=[] #contains mottif headers
b=[] #contains mottif sequences
header_counter=-1
with open(sys.argv[1], 'r') as mottif_file: #readfile from input argument 1 as mottif data
	for line in mottif_file.readlines():
		line = line.rstrip('\n') #remove "\n" 
		if re.match('>',line): #detecting headers from fasta file
			a.append(line) #adding mottif header line data
			flag=0 #flag used to identify the next line of data as sequence data
			header_counter = header_counter + 1 #counting the number of headers in mottif file
		elif flag==0:
			b.append(line) #adding sequence data
			flag=1 #flag to continue adding sequence data until next ">"
		else:
			b[header_counter] = b[header_counter] + line #adding sequence data; this method was used over append as append added "," between entries
maxb = len(max(b)) #variable used later in loop; set to the max length of the longest detected mottif sequence
c=[] #contains sequence headers
d=[] #contains sequence sequences/data
Gcount=0 #variable containing count of G nucleotide
Ccount=0 #variable containing count of C nucleotide
Acount=0 #variable containing count of A nucleotide
Tcount=0 #variable containing count of T nucleotide 
Totalcount=0 #variable containing count of all nucleotide
header_counter=-1
with open(sys.argv[2], 'r') as sequence_file: #readfile from input argument 2 as sequence data
	for line in sequence_file.readlines():
		line=line.rstrip() #remove "\n" 
		if re.match('>',line): #detecting headers from fasta file
			c.append(line) #adding sequence header line
			flag=0 #flag used to identify the next line of data as sequence data
			header_counter = header_counter + 1 #counting the number of headers in sequence file
		elif flag==0:
			d.append(line) #adding sequence data
			flag=1
			Gcount = Gcount + (line.count('G')) #counting G nucleotide
			Ccount = Ccount + (line.count('C')) #counting C nucleotide
			Acount = Acount + (line.count('A')) #counting A nucleotide
			Tcount = Tcount + (line.count('T')) #counting T nucleotide
			Totalcount = Totalcount + (line.count('')) #counting Total nucleotide
		else:
			d[header_counter] = d[header_counter] + line #adding sequence data; this method was used over append as append added "," between entries
			Gcount = Gcount + (line.count('G')) #counting G nucleotide
			Ccount = Ccount + (line.count('C')) #counting C nucleotide
			Acount = Acount + (line.count('A')) #counting A nucleotide
			Tcount = Tcount + (line.count('T')) #counting T nucleotide
			Totalcount = Totalcount + (line.count('')) #counting Total nucleotide
if re.match('\d+%',sys.argv[3]):
	prescore=sys.argv[3]
	score=prescore.rstrip(("%"))
else:
	score = int(sys.argv[3]) #score from input argument 3
maxd = len(max(d)) #variable used later in loop; set to the max length of the longest detected mottif sequence
Gpercent = Gcount/Totalcount #calculate % of Genome consisting of G
Cpercent = Ccount/Totalcount #calculate % of Genome consisting of C
Apercent = Acount/Totalcount #calculate % of Genome consisting of A
Tpercent = Tcount/Totalcount #calculate % of Genome consisting of T
x=0 #variable for next loop
y=0 #variable for next loop
Matrix = [[0.00 for i in range(4)] for i in range(maxb)] #create array with 4 columns and number of rows = longest mottif sequence length; 0.00 is value for all entry
while (x< maxb):
	Matrix[y][0] = Apercent #set value of first column to background % of A in test sequence
	Matrix[y][1] = Cpercent #set value of first column to background % of C in test sequence
	Matrix[y][2] = Gpercent #set value of first column to background % of G in test sequence
	Matrix[y][3] = Tpercent #set value of first column to background % of T in test sequence
	x = x + 1
	y = y + 1 
x=0 #variable for next loop
while (x < len(b)): #note currently only lower case
	y = 0 #variable for loop
	for sequence in b[x]: #note these are currently set to lower case
		if sequence == 'a':
			Matrix[y][0] = Matrix[y][0] + 1 #counting number of A nucleotide in motiff sequence
		if sequence == 'c':
			Matrix[y][1] = Matrix[y][1] + 1 #counting number of C nucleotide in motiff sequence
		if sequence == 'g':
			Matrix[y][2] = Matrix[y][2] + 1 #counting number of G nucleotide in motiff sequence
		if sequence == 't':
			Matrix[y][3] = Matrix[y][3] + 1 #counting number of T nucleotide in motiff sequence
		y = y + 1
	x = x + 1
i = 0
print ("@@Count of nucleotides by position in motiff file")
print ("{:<7}{:<7}{:<7}{:<7}{:<7}".format("POS","A","C","G","T"))
while (i < maxb):	
	print ("{:3}{!s:>7.5}{!s:>7.5}{!s:>7.5}{!s:>7.5}".format(i+1,Matrix[i][0],Matrix[i][1],Matrix[i][2],Matrix[i][3])) #print count of nucleotides in mottif
	i = i + 1
f = 0
Matrix2 = [[0 for i in range(4)] for i in range(maxb)] 
while (f < maxb): #calculate frequency of test sequence nucleotides
	Matrix2[f][0] = Matrix[f][0]/(Matrix[f][0]+Matrix[f][1]+Matrix[f][2]+Matrix[f][3]) 
	Matrix2[f][1] = Matrix[f][1]/(Matrix[f][0]+Matrix[f][1]+Matrix[f][2]+Matrix[f][3]) 
	Matrix2[f][2] = Matrix[f][2]/(Matrix[f][0]+Matrix[f][1]+Matrix[f][2]+Matrix[f][3]) 
	Matrix2[f][3] = Matrix[f][3]/(Matrix[f][0]+Matrix[f][1]+Matrix[f][2]+Matrix[f][3]) 
	f = f + 1
f = maxb - 1
g=0
Matrix3 = [[0 for i in range(4)] for i in range(maxb)] 
while (f > -1): #invert instead of reverse compliment the test sequence to test for both directions
	Matrix3[f][0] = Matrix[g][0]/(Matrix[g][0]+Matrix[g][1]+Matrix[g][2]+Matrix[g][3]) 
	Matrix3[f][1] = Matrix[g][1]/(Matrix[g][0]+Matrix[g][1]+Matrix[g][2]+Matrix[g][3]) 
	Matrix3[f][2] = Matrix[g][2]/(Matrix[g][0]+Matrix[g][1]+Matrix[g][2]+Matrix[g][3]) 
	Matrix3[f][3] = Matrix[g][3]/(Matrix[g][0]+Matrix[g][1]+Matrix[g][2]+Matrix[g][3]) 
	g = g + 1
	f = f - 1
	
print ("@@ Frequency Matrix")
print ("{:<7}{:<7}{:<7}{:<7}{:<7}".format("POS","A","C","G","T"))
h = 0
while (h < maxb):	
	print ("{:3}{!s:>7.5}{!s:>7.5}{!s:>7.5}{!s:>7.5}".format(h+1,Matrix2[h][0],Matrix2[h][1],Matrix2[h][2],Matrix2[h][3])) #print frequency of nucleotide of test sequence

	h = h + 1
z=0
print ("@@PSSM")
print ("{:<7}{:<7}{:<7}{:<7}{:<7}".format("POS","A","C","G","T"))
while (z < maxb): #print PSSM
	print ("{:3}{!s:>7.5}{!s:>7.5}{!s:>7.5}{!s:>7.5}".format(z + 1,(math.log2((Matrix2[z][0]) / Apercent)),(math.log2((Matrix2[z][1]) / Cpercent)),(math.log2((Matrix2[z][2]) / Gpercent)),(math.log2((Matrix2[z][3]) / Tpercent))))
	z = z+1


start=0
finish=maxb #sliding window set to length of longest mottif sequence
while (finish<maxd): 
	m=0
	thing = (d[0][start:finish]) #get the nucleotides from the test sequence between start position and end position
	x=0
	x2 = 0 
	for thing2 in thing:
		if thing2 == 'A':
			s =(math.log2((Matrix2[m][0]) / Apercent))
			s2 = (math.log2((Matrix3[m][3]) / Tpercent))
		if thing2 == 'C':
			s =(math.log2((Matrix2[m][1]) / Cpercent))
			s2 = (math.log2((Matrix3[m][2]) / Gpercent))
		if thing2 == 'G':
			s =(math.log2((Matrix2[m][2]) / Gpercent))
			s2=(math.log2((Matrix3[m][1]) / Cpercent))
		if thing2 == 'T':
			s =(math.log2((Matrix2[m][3]) / Tpercent))
			s2=(math.log2((Matrix3[m][0]) / Apercent))
		m = m + 1
		x = x + s
		x2 = x2 + s2
	
	if x > score: #{:18} 18 is length to display, {:^5d} center align 5 digits
		print ("{:18}{:^8d}{:^11}{:8d}{:^7}{:5.1f}{:^10}{:}".format("Match starting at",start+1, "Ending at", finish, "Score", x, "Sequence", d[0][start:finish]))
	if x2 > score:
		print ("{:18}{:^8d}{:^11}{:8d}{:^7}{:5.1f}{:^10}{:}".format("Match starting at",start+1, "Ending at", finish, "Score", x2, "Sequence", d[0][start:finish]))
	start = start + 1
	finish = finish + 1
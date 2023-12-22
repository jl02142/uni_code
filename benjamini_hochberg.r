# J Ludwig 2019
# Benjamini Hochberg multiple test correction
# To show students how easy it is to write their own code instead of spending more time installing a library and reading the documentation

#read data
# need to update code with name of file below
GeneExpression <- read.csv("") 

#Set normal and cancer data into their own
norm = GeneExpression[1:20,2:101]
canc = GeneExpression[21:40,2:101]

#Init t & loop through data conducting t-test and storing the p-value into t
t=list()
for(i in 1:100){ t[i] = (t.test(norm[i], canc[i])$p.value)}

#The function unlist is used to fix dumb list errors (I'm being lazy and should have stored the data differently but this works)
#BH requires the data to be sorted from smallest to largest so use sort function
sort_t=(sort(unlist(t)))

#Init bh_t to store data in
bh_t=list()

#This is the Benjamini�Hochberg procedure to calculate the new crit values
# https://www.statisticshowto.com/benjamini-hochberg-procedure/
for(i in 1:100){bh_t[i]=(i/100)*0.05}

#This is part of the Benjamini�Hochberg procedure to determine where the cut off is. 
#The 13th element is the cutoff, the first 12 pass from the class data
for(i in 1:100){if(sort_t[i] >= bh_t[i]){print(i) & break}}

#This tells you what gene numbers in the original data set have small enough p-values to be considered significant post BH
for(i in 1:12){for(j in 1:100){if(sort_t[i] == t[j]){
	cat(sprintf("%s\ %e\ \n", j, t[[j]]))
}}}

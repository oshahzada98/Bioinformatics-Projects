import math 

'''
    Owais A. Shahzada 
    
    October 22, 2018
    
    Objective: Determine whether the protein is "extended" or "collapsed" 
    
    Description: Each amino acid has its own physical properties. Proteins 
    physical properties can determine whether the protein sequence will be 
    collapsed or extended. Using the Average Net Charge and the Average 
    Hydrophobicity table we must predict whether or not the protein will 
    be extended or collapseed. 
    
    Equation: 
    raw_prediction = 0.343 x abs(ANC) - 0.940 x AH + 0.380 
    prediction = 1/(1 + exp(-20 * raw_prediction))

'''

 
file_FastaFile = "mt_fasta_nn_2.txt"
file_HydroTable = "hydrophobicity_table.txt"
file_Translation = "translation_table.txt"
print("\n")

file1 = open(file_FastaFile)
file2 = open(file_HydroTable)
file3 = open(file_Translation)

seqID = 0    
dnaSeq = 1    
char = 0     
hydroVal = 1 
rnaCod = 0   
trans = 1    

fastaDict = {}       
hydroDict = {}        
translationDict = {} 


for line in file2: 
    lineSplit2 = line.strip().split("\t") 
    hydroDict[lineSplit2[char]] = float(lineSplit2[hydroVal]) 


for line in file3:
    lineSplit3 = line.strip().split("\t") 
    translationDict[lineSplit3[rnaCod]] = lineSplit3[trans]

def calculateValues(lineSplit): 
    seqIndex = 0     
    dnaIndex = 1     
    seqID = lineSplit[seqIndex] 
    dnaSeq = lineSplit[dnaIndex] 

  
    rnaNuc = ""
    for character in dnaSeq:  
        if character == "C":    
            rnaNuc += "C"       
        if character == "G":   
            rnaNuc += "G"       
        if character == "A":
            rnaNuc += "A"
        if character == "T":
            rnaNuc += "U"
    
      
    transCodon = ""
    start = rnaNuc.find("AUG") 
    for i in range(start, len(rnaNuc), 3): 
        codon = rnaNuc[i: (i+3)] 
        print(translationDict)
        if translationDict[codon] == "*": 
            break
        else:
           transCodon += translationDict[codon]  

    count_R = transCodon.count("R")
    count_D = transCodon.count("D")
    count_H = transCodon.count("H")
    count_K = transCodon.count("K")
    count_E = transCodon.count("E")
    protein_Len = len(transCodon) 


    ancVal = (count_K + count_R + 0.5 * count_H - count_D - count_E) / protein_Len
    ancVal = round(ancVal, 4)
    
    sum = 0 
    for aminoAcid in transCodon:  
      sum += hydroDict[aminoAcid] 

    avgHydro = sum / protein_Len 
    avgHydro = round(avgHydro, 4)

    raw_Predict = 0.343 * abs(ancVal) - 0.940 * avgHydro + 0.380
    prediction = 1/(1+ math.exp(-20 * raw_Predict))
    prediction = round(prediction, 4)
    if (prediction >= 0.5):
        proteinClass = "extended"
    else:
        proteinClass = "collapsed"
    


print("seq_id", "hydro", "charge", "value", "class")
sequenceText = file1.readlines() 
newSeqText = []  

for line in sequenceText:
    
    if ">" in line:    
        newSeqText.append(line) 
    else:
        newLine = line.replace("\n", "")   
        newSeqText.append(newLine)        

joinedText = "".join(newSeqText)            
replacedText = joinedText.replace("\n", " ")
replacedText = replacedText.split(">")      
del replacedText[0]                        


for seq in replacedText:                    
    lineSplit = seq.split(" ")
    print(lineSplit)
    fastaDict[lineSplit[seqId]] = lineSplit[dnaSeq]
    print(fastaDict)
    calculateValues(lineSplit) 

from textwrap import dedent
import numpy as np 
import re 
import sys 


'''
    Owais A. SHahzada 
    02/19/2019
    Objective: Build the Needleman-Wunsch Algorithm 
    Description: Needleman-Wunsch is an algorithm to align 
    nucleotide or protein sequences. 
    
'''

#reverse sequences function
def reverse_sequence_B(alignB):
    reversed_B = alignB[::-1]

    return reversed_B

def reverse_sequence_A(alignA):
    reversed_A = alignA[::-1]
    
    return reversed_A

#Traceback function 
def traceback(seq1,seq2, matrix, score_table, gap,score1):

    i = len(seq1)
    j = len(seq2)

    alignA = ""
    alignB = ""

    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i,j] == matrix[i-1,j-1] + score_table[i,j]:
            alignA += seq1[i-1]
            alignB += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0 and matrix[i,j] == matrix[i-1][j] + gap:
            alignA += seq1[i-1]
            alignB += "-"
            i -= 1
        else:
            alignA += "-"
            alignB += seq2[j-1]
            j -= 1
    aligned_sequence_output_A = reverse_sequence_A(alignA)
    aligned_sequence_output_B = reverse_sequence_B(alignB)
    
    print(dedent(f"""
{matrix}
Max score: {score1}
{aligned_sequence_output_A}
{aligned_sequence_output_B}
    """))
    
# Creates matrix 
def needleman_matrix(seq1, seq2):
    
    #Penalties 
    match = 1
    mismatch = -1
    gap = -2

    N,M = len(seq1), len(seq2)

    #2d array of zeros 
    array_of_zeros = np.zeros((N+1,M+1))    
    matrix = array_of_zeros           # Keeps track of all the paths taken
    score_table = np.zeros((N+1,M+1)) # Keeps track of all the scores 
    
    for i in range(N+1):
        matrix[i][0] = gap * i
    for j in range(M+1):
        matrix[0][j] = gap * j
        
    for i in range(1,N+1):
        for j in range(1,M+1):
            if seq1[i-1] == seq2[j-1]: # Condtion checked if the nucleotides of the two sequeces match
                score = match
            else:
                score = mismatch 
            score1 = matrix[i-1][j-1] + score # Depending on the match or mismatch score1 will store that
            score2 = matrix[i][j-1] + gap     # Score2 reports back the score in an upwards direction
            score3 = matrix[i-1][j] + gap     # Score3 reprots back a score in the leftward direction

            matrix[i,j] = max(score1,score2,score3) #reports back the max score of the three 
            score_table[i,j] = score
    traceback(seq1,seq2,matrix,score_table,gap,score1)


def main():
    #Catches the seqeunces and stores into a list
    with open(sys.argv[1], 'r') as fasta: 
        text = fasta.read()         
        replaced_text = text.replace("\n", " ")
        pattern = r"(\b[GATCRY\n]+)\b" #Regex pattern
        sequences = re.findall(pattern, replaced_text)
        seq1 = sequences[0] 
        seq2 = sequences[1]

    needleman_matrix(seq1, seq2)
    

if __name__ == "__main__":
    main()
        

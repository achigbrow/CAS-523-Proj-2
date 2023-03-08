# A_problems.py
 
import math

# permutation = combinations * 3^k
def choose(n,k):
  comb = math.factorial(n)/(math.factorial(k)*math.factorial(n-k))
  return comb * pow(3,k)

# A-1-1
print('# of permutations after genome is mutated once')
print('full=',perms(29903,1))
print('spike=',perms(3138,1))

# A-1-3 
print('# of permutations after genome is mutated thrice')
print('full=',perms(29903,3))
print('spike=',perms(3138,3))

# A-1-2

print('# of permutations with different amino acid after genome is mutated once')

# triplet -> amino acid
'''
UUU	F		UCU	S		UAU	Y		UGU	C
UUC	F		UCC	S		UAC	Y		UGC	C
UUA	L		UCA	S		UAA	$		UGA	$
UUG	L		UCG	S		UAG	$		UGG	W
CUU	L		CCU	P		CAU	H		CGU	R
CUC	L		CCC	P		CAC	H		CGC	R
CUA	L		CCA	P		CAA	Q		CGA	R
CUG	L		CCG	P		CAG	Q		CGG	R
AUU	I		ACU	T		AAU	N		AGU	S
AUC	I		ACC	T		AAC	N		AGC	S
AUA	I		ACA	T		AAA	K		AGA	R
AUG	M		ACG	T		AAG	K		AGG	R
GUU	V		GCU	A		GAU	D		GGU	G
GUC	V		GCC	A		GAC	D		GGC	G
GUA	V		GCA	A		GAA	E		GGA	G
GUG	V		GCG	A		GAG	E		GGG	G
'''


amino_acids = {
'UUU' : 'F', 'UCU' : 'S', 'UAU' : 'Y', 'UGU' : 'C', 'UUC' : 'F', 'UCC' : 'S', 'UAC' : 'Y', 'UGC' : 'C', 'UUA' : 'L', 'UCA' : 'S', 'UAA' : '$', 'UGA' : '$', 'UUG' : 'L', 'UCG' : 'S', 'UAG' : '$', 'UGG' : 'W', 'CUU' : 'L', 'CCU' : 'P', 'CAU' : 'H', 'CGU' : 'R', 'CUC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CUA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CUG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'AUU' : 'I', 'ACU' : 'T', 'AAU' : 'N', 'AGU' : 'S', 'AUC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'AUA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'AUG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GUU' : 'V', 'GCU' : 'A', 'GAU' : 'D', 'GGU' : 'G', 'GUC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GUA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GUG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'
}
    
S = 'ACDEFGHIKLMNPQRSTVWY$' # $ is stop

hist = [0] * 21

for k,v in amino_acids.items():
  idx = S.find(k)
  print(idx)
  hist[idx] = hist[idx]+1 
    
print(hist)

for i in range(1,21):
  hist[i] = hist[i]/sum(hist)

print('probability of amino acid occurring')
print(hist)

L = 3138/3

#M = pow(20,L) 






import math

# permutation = combinations * 3^k
def perms(n,k):
  comb = math.factorial(n)/(math.factorial(k)*math.factorial(n-k))
  return comb * pow(3,k)

# A-1-1
print('# of permutations after genome is mutated once')
print('full=',perms(29903,1))
print('spike=',perms(3138,1))

# A-1-2

print('# of permutations with different amino acid after genome is mutated once')

# amino acid table's probabilities ...
2/64 # Phe
2/64 # 
...
4/64 # 



# A-1-3
print('# of permutations after genome is mutated thrice')
print('full=',perms(29903,3))
print('spike=',perms(3138,3))

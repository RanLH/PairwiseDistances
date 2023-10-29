# ed.py 
# Ran Libeskind-Hadas
# July 15, 2023
# Functions using memoization to...
#   1. Compute the Levenstein edit distance of a pair of strings and the number
#       of optimal solutions.
#   2. Compute the diameter of the space of optimal solutions.
#   3. Compute the pairwise-distance vector of the space of optimal solutions.

from Vector import *
import random
import sys
import matplotlib.pyplot as plt

sys.setrecursionlimit(10000)
# Directions of move in the dynamic programming (or memoization) table.
LEFT = (0, -1)
UP = (-1, 0)
DIAGONAL = (-1, -1)

def EDcount(S1, S2, memo = {}):
    """ Takes as input two strings and returns an ordered pair comprising 
        the optimal edit distance and the number of solutions with that optimal score. """
    if (S1, S2) in memo: return memo[(S1, S2)]
    if S1 == "": return (len(S2), 1)
    elif S2 == "": return (len(S1), 1)
    elif S1[-1] == S2[-1]:
        solution = EDcount(S1[:-1], S2[:-1], memo)
        memo[(S1, S2)] = solution
        return solution
    else:
        num = 0
        option1, num1 = EDcount(S1, S2[:-1])
        option2, num2 = EDcount(S1[:-1], S2, memo)
        option3, num3 = EDcount(S1[:-1], S2[:-1], memo)
        optimal = min(option1, option2, option3)
        if option1 == optimal:
            num += num1
        if option2 == optimal:
            num += num2
        if option3 == optimal:
            num += num3
        solution = (1 + optimal, num)
        memo[(S1, S2)] = solution
        return solution


def EDalignments(S1, S2, memo = {}):
    """ Takes two strings as input and returns an ordered pair comprising the 
        optimal edit distance and a list of all of alignments with that score. """
    if (S1, S2) in memo: return memo[(S1, S2)]
    if S1 == "": return (len(S2), [("_"*len(S2), S2)])
    elif S2 == "": return (len(S1), [(S1, "_"*len(S1))])
    elif S1[-1] == S2[-1]:
        score, subsolutions = EDalignments(S1[:-1], S2[:-1])
        solutions = [(x + S1[-1], y + S2[-1]) for (x, y) in subsolutions]
        package = (score, solutions)
        memo[(S1, S2)] = package
        return package
    else:
        alignments = []
        option1, alignments1 = EDalignments(S1, S2[:-1])
        option2, alignments2 = EDalignments(S1[:-1], S2)
        option3, alignments3 = EDalignments(S1[:-1], S2[:-1])
        optimal = min(option1, option2, option3)
        if option1 == optimal:
            alignments.extend([(x+"_", y+S2[-1]) for (x, y) in alignments1])
        if option2 == optimal:
            alignments.extend([(x+S1[-1], y+"_") for (x, y) in alignments2])
        if option3 == optimal:
            alignments.extend([(x + S1[-1], y + S2[-1]) for (x, y) in alignments3])
        package = (optimal+1, alignments)
        memo[(S1, S2)] = package
        return package

def EDannotate(S1, S2, memo):
    """ Takes as input two strings and returns an ordered pair comprising the 
        optimal edit distance.
        As a side-effect, populates the memo dictionary with 
        a tuple (d, numSolutions, L) where d is the edit distance, numSolutions
        is the number of optimal solutions, and L is the list of pointers
        for optimal solutions, where pointer can be
        LEFT, UP, or DIAGONAL.
    """
    if (S1, S2) in memo: return memo[(S1, S2)]
    if S1 == "" and S2 == "": return [0, 1, []]
    elif S1 == "" and S2 != "": 
        subProblemSolution = EDannotate(S1, S2[:-1], memo)
        score, numSolutions, _ = subProblemSolution
        solution = [score + 1, numSolutions, [LEFT]]
        memo[(S1, S2)] = solution
        return solution
    elif S1 != "" and S2 == "": 
        subProblemSolution = EDannotate(S1[:-1], S2, memo)
        score, numSolutions, _ = subProblemSolution
        solution = [score + 1, numSolutions, [UP]]
        memo[(S1, S2)] = solution
        return solution
    elif S1[-1] == S2[-1]:
        subProblemSolution = EDannotate(S1[:-1], S2[:-1], memo)
        score, numSolutions, _ = subProblemSolution
        solution = [score, numSolutions, [DIAGONAL]]
        memo[(S1, S2)] = solution
        return solution
    else:
        num = 0
        option1, num1, _ = EDannotate(S1, S2[:-1], memo)
        option2, num2, _ = EDannotate(S1[:-1], S2, memo)
        option3, num3, _ = EDannotate(S1[:-1], S2[:-1], memo)
        optimal = min(option1, option2, option3)
        directions = []
        if option1 == optimal:
            num += num1
            directions.append("LEFT")
        if option2 == optimal:
            num += num2
            directions.append("UP")
        if option3 == optimal:
            num += num3
            directions.append("DIAGONAL")
        solution = [1 + optimal, num, directions]
        memo[(S1, S2)] = solution
        return solution
    
def subsumes(i, j, k, l):
    # Returns true iff ordered pair (i, j) subsumes ordered pair (k, l)
    return (i > k and j >= l) or (i >=k and j > l)
    
def diameter(i, j, k, l, DP1, DP2 = {}):
    """ Takes as input ordered pair i, j and ordered pair k, l representing
    the indices in the DP table (DP1) for the alignment of two strings.
    Uses DP2 as the memoization table to compute the diameter of the solution
    space - that is - the largest number of differences in edges between a
    path starting at location (i, j) and a path starting at location (k, l)
    in DP table DP1.
    """
    if i == 0 and j == 0 and k == 0 and l == 0: return 0
    if (i, j, k, l) in DP2: return DP2[(i, j, k, l)]
    else:
        if i == 0 and j == 0: listOfDirections1 = []
        else:
            _, _, listOfDirections1 = DP1[(i, j)]
        if k == 0 and l == 0: listOfDirections2 = []
        else: _, _, listOfDirections2 = DP1[(k, l)]
        options = []
        if subsumes(i, j, k, l):
            for direction in listOfDirections1:
                d1, d2 = direction
                q = i + d1
                r = j + d2
                option = 1 + diameter(q, r, k, l, DP1, DP2)
                options.append(option)
            solution = max(options)
            DP2[(i, j, k, l)] = solution
            return solution
        elif subsumes(k, l, i, j) or (i != k and j != l):
            for direction in listOfDirections2:
                d1, d2 = direction
                q = k + d1
                r = l + d2
                option = 1 + diameter(i, j, q, r, DP1, DP2)
                options.append(option)
            solution = max(options)
            DP2[(i, j, k, l)] = solution
            return solution
        else: # i == k and j == l
            for direction1 in listOfDirections1:
                dx1, dy1 = direction1
                for direction2 in listOfDirections2:
                    dx2, dy2 = direction2
                    if direction1 == direction2: add = 0
                    else: add = 2
                    option = add + diameter(i+dx1, j+dy1, k+dx2, l+dy2, DP1, DP2)
                    options.append(option)
            solution = max(options)
            DP2[(i, j, k, l)] = solution
            return solution

def EDannotate1(S1, S2, len1, len2, memo):
    """ Takes as input two strings, S1 and S2, and the lengths of the 
        prefixes of those strings to align (len1 for S1 and len2 for S2).
        Note that S1 and S2 never change in this function. Instead, we 
        vary the values of len1 and len2 to indicate the lengths of the prefixes
        of those strings that we seek to align.
        Returns an ordered pair comprising the optimal edit distance.
        As a side-effect, populates the memo dictionary with 
        a tuple (d, numSolutions, L) where d is the edit distance, numSolutions
        is the number of optimal solutions, and L is the list of pointers
        for optimal solutions, where pointer can be
        LEFT, UP, or DIAGONAL
    """
    if (len1, len2) in memo: return memo[(len1, len2)]
    if len1 == 0 and len2==0: return [0, 1, []]
    elif len1 == 0 and len2 > 0: 
        subProblemSolution = EDannotate1(S1, S2, len1, len2-1, memo)
        score, numSolutions, _ = subProblemSolution
        solution = [score + 1, numSolutions, [LEFT]]
        memo[(len1, len2)] = solution
        return solution
    elif len1 > 0 and len2 == 0: 
        subProblemSolution = EDannotate1(S1, S2, len1-1, len2, memo)
        score, numSolutions, _ = subProblemSolution
        solution = [score + 1, numSolutions, [UP]]
        memo[(len1, len2)] = solution
        return solution
    elif S1[len1-1] == S2[len2-1]:
        subProblemSolution = EDannotate1(S1, S2, len1-1, len2-1, memo)
        score, numSolutions, _ = subProblemSolution
        solution = [score, numSolutions, [DIAGONAL]]
        memo[(len1, len2)] = solution
        return solution
    else:
        num = 0
        option1, num1, _ = EDannotate1(S1, S2, len1, len2-1, memo)
        option2, num2, _ = EDannotate1(S1, S2, len1-1, len2, memo)
        option3, num3, _ = EDannotate1(S1, S2, len1-1, len2-1, memo)
        optimal = min(option1, option2, option3)
        directions = []
        if option1 == optimal:
            num += num1
            directions.append(LEFT)
        if option2 == optimal:
            num += num2
            directions.append(UP)
        if option3 == optimal:
            num += num3
            directions.append(DIAGONAL)
        solution = [1 + optimal, num, directions]
        memo[(len1, len2)] = solution
        return solution
    
def pdv(i, j, k, l, dim, DP1, DP2 = {}):
    """ Takes as input the ordered pair i, j and ordered pair k, l,
    the DP table (DP1) for the edit distance between two strings, the dimension
    of the output vector (which should be at least one more than the
    diameter of the solution space, and uses DP2 to as the memoization table
    to compute the pairwise distance vector. """
    if i == 0 and j == 0 and k == 0 and l == 0: 
        v = Vector(dim)
        v.set(0, 1)
        return v
    if (i, j, k, l) in DP2: return DP2[(i, j, k, l)]
    else:
        if i == 0 and j == 0: listOfDirections1 = []
        else:
            _, _, listOfDirections1 = DP1[(i, j)]
        if k == 0 and l == 0: listOfDirections2 = []
        else: _, _, listOfDirections2 = DP1[(k, l)]
        solution = Vector(dim)
        if subsumes(i, j, k, l):
            for direction in listOfDirections1:
                d1, d2 = direction
                q = i + d1
                r = j + d2
                option = pdv(q, r, k, l, dim, DP1, DP2)
                option = option.shift(1)
                solution = solution + option
            DP2[(i, j, k, l)] = solution
            return solution
        elif subsumes(k, l, i, j) or (i != k and j != l):
            for direction in listOfDirections2:
                d1, d2 = direction
                q = k + d1
                r = l + d2
                option = pdv(i, j, q, r, dim, DP1, DP2)
                option = option.shift(1)
                solution = solution + option
            DP2[(i, j, k, l)] = solution
            return solution
        else: # i == k and j == l
            for dir1 in range(len(listOfDirections1)):
                dx1, dy1 = listOfDirections1[dir1]
                for dir2 in range(dir1, len(listOfDirections2)):
                    dx2, dy2 = listOfDirections2[dir2]
                    if dir1 == dir2: add = 0
                    else: add = 2
                    option = pdv(i+dx1, j+dy1, k+dx2, l+dy2, dim, DP1, DP2)
                    option = option.shift(add)
                    solution = solution + option
            DP2[(i, j, k, l)] = solution
            return solution
    
def formatAlignments(alignmentList):
    for (s1, s2) in alignmentList:
        print(s1)
        print(s2)
        print()

def main(S1, S2):
    score, alignmentList = EDalignments(S1, S2)
    print("Alignment score = ", score)
    formatAlignments(alignmentList)

def experiment(numTrials, seqLength):
    symbols = ['A', 'T', 'C', 'G']
    numSolutions = []
    for _ in range(numTrials):
        string1 = ""
        string2 = ""
        for _ in range(seqLength):
            string1 = string1 + random.choice(symbols)
            string2 = string2 + random.choice(symbols)
        dist, num = EDcount(string1, string2)
        numSolutions.append(num)
    return min(numSolutions), max(numSolutions)

length = 100
bases = ["A", "T", "C", "G"]
S1 = "".join(random.choices(bases, k=length))
S2 = "".join(random.choices(bases, k=length))
dim = 2*(len(S1) + len(S2)) # Upper-bound on diameter
DP1 = {}
EDannotate1(S1, S2, len(S1), len(S2), DP1)
print("S1 = ", S1)
print("S2 = ", S2)
print("-----------------")
print("DP table")
print(DP1)
print("-----------------")
diam = diameter(len(S1), len(S2), len(S1), len(S2), DP1, {})
print("Diameter: ", diam)
print("-----------------")
print("Pairwise Distance Vector")
vector = pdv(len(S1), len(S2), len(S1), len(S2), dim, DP1, {})
print(vector)
plt.bar(x = list(range(0, diam+1)), height = vector.v[:diam+1])
plt.show()

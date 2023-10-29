# trees.py 
# Ran Libeskind-Hadas
# October 2023
# Functions using memoization to...
#   1. Compute solve the small parsimony problem.
#   2. Compute the diameter of the space of optimal solutions.
#   3. Compute the pairwise-distance vector of the space of optimal solutions.

import math
import random
from Vector import *
import matplotlib.pyplot as plt

# A phylogenetic tree is represented as dictionary where each key is a node where
# is a non-negative integer, with 0 designating the root. The value associated with each
# node key is a list of the children of that node. That list has either two entries (in
# the case of an internal node) or zero entries (in the case of a leaf).

tree1 = {0: [1, 2], 1: [3, 4], 2: [5, 6], 3: [], 4: [], 5: [], 6: []}

# The characters under consideration are represented by a set.

characters1 = {"A", "T", "C", "G"}

# Each leaf in the tree has an associated character label, represented by dictionary
# in which the key is the name of the leaf and the value is the character label

leaves1 = {3: "A", 4: "A", 5: "T", 6: "G"}

# Another test set...

tree2 = {0: [1, 2], 1: [3, 4], 3: [5, 6], 5: [7, 8], 2: [], 4: [], 6: [], 7: [], 8: []}
characters2 = [0, 1, 2]
leaves2 = {2: 1, 4: 0, 6: 0, 7: 1, 8: 2}

# We seek to compute a maximum parsimony labeling for the tree. Once computed, this is stored
# in a dictionary in which the keys are internal nodes and the values are themselves
# character dictionaries.
# In a character dictionary, the key is a character from the character set and the value
# is dictionary with keys: "score", "num solutions", "left set", "right set"
# where score is the optimal score, num solutions is the number of solutions with this score
# left set is the set of characters that can be found on the left child in such optimal solutions
# and right set is the set of characters that can be found the right child in such optimal solutions.

def leaf(node, tree):
    """ Takes a node and tree as input and returns True if the node is a leaf and false
        otherwise. """
    return tree[node] == []

def diff(x, y):
    """ Returns 0 if x == y and 1 otherwise. """
    if x == y: return 0
    else: return 1

def maxparsimony(node, tree, characters, leaves, labeling):
    """ Takes as input a node, a phylogenetic tree, a list of characters, a leaves dictionary 
        and returns a labeling of the internal nodes. 
        Returns nothing, but updates the labeling as a side-effecct. """
    if leaf(node, tree): 
        labeling[node] = {}
        for char in characters:
            if char == leaves[node]: 
                labeling[node][char] = {"score": 0, "num solutions": 1, "left set": {}, "right set": {}}
            else: 
                labeling[node][char] = {"score": math.inf, "num solutions": 0, "left set": {}, "right set":{}}
    else:
        leftChild, rightChild = tree[node]
        maxparsimony(leftChild, tree, characters, leaves, labeling)
        maxparsimony(rightChild, tree, characters, leaves, labeling)
        labeling[node] = {}
        for character in characters:
            candidates = []
            for leftLabel in characters:
                leftDict = labeling[leftChild][leftLabel]
                leftScore = leftDict["score"]
                leftNumSolutions = leftDict["num solutions"]
                for rightLabel in characters:
                    rightDict = labeling[rightChild][rightLabel]
                    rightScore = rightDict["score"]
                    rightNumSolutions = rightDict["num solutions"]
                    score = diff(character, leftLabel) + diff(character, rightLabel) + leftScore + rightScore
                    candidates.append((score, leftLabel, rightLabel, leftNumSolutions * rightNumSolutions))
            optimalScoreForChar = min(candidates)[0]    # Optimal score for the current character
            leftChildSet = set()
            rightChildSet = set()
            numSolutions = 0
            for candidate in candidates:
                if candidate[0] == optimalScoreForChar:
                    leftChildSet.add(candidate[1])
                    rightChildSet.add(candidate[2])
                    numSolutions += candidate[3]
            labeling[node][character] = {"score": optimalScoreForChar, "num solutions": numSolutions, "left set": leftChildSet, "right set": rightChildSet}

def diameter(tree, labeling):
    """ Takes a tree and labeling (as computed by the maxparsimony function) and returns
        the diameter of the space of optimal solutions. """
    optimum = min([labeling[0][c]["score"] for c in labeling[0]])
    optimumRootChars = [c for c in labeling[0] if labeling[0][c]["score"] == optimum]
    candidates = []
    for char1 in optimumRootChars:
        for char2 in optimumRootChars:
            candidates.append(maxdist(0,char1, char2, tree, labeling))
    diameter = max(candidates)
    return diameter

def maxdist(node, char1, char2, tree, labeling, memo={}):
    """ Helper function used by the diameter function above. """
    if leaf(node, tree): return 0
    if (node, char1, char2) in memo: return memo[(node, char1, char2)]
    leftChild, rightChild = tree[node]
    leftLabels1 = labeling[node][char1]["left set"]
    rightLabels1 = labeling[node][char1]["right set"]
    leftLabels2 = labeling[node][char2]["left set"]
    rightLabels2 = labeling[node][char2]["right set"]
    leftOptions = []
    for c1 in leftLabels1:
        for c2 in leftLabels2:
            leftOptions.append(maxdist(leftChild, c1, c2, tree, labeling))
    rightOptions = []
    for c1 in rightLabels1:
        for c2 in rightLabels2:
            rightOptions.append(maxdist(rightChild, c1, c2, tree, labeling))
    solution = max(leftOptions) + max(rightOptions) + diff(char1, char2)
    memo[(node, char1, char2)] = solution
    return solution

def pdv(tree, labeling):
    """ Takes a tree and labeling (as computed by the maxparsimony function) and returns
        the pairwise distance vector for the space of solutions. """
    optimum = min([labeling[0][c]["score"] for c in labeling[0]])
    optimumRootChars = set([c for c in labeling[0] if labeling[0][c]["score"] == optimum])
    candidates = []
    dim = (len(tree)+1)//2
    solution = Vector(dim)
    pairs = set()
    for char1 in optimumRootChars:
        for char2 in optimumRootChars:
            if (char1, char2) in pairs or (char2, char1) in pairs: continue
            solution = solution + pdvHelper(0, char1, char2, tree, labeling)
            pairs.add((char1, char2))
    return solution

def pdvHelper(node, char1, char2, tree, labeling, memo={}):
    """ Helper function used by pdv above. """
    dim = (len(tree)+1)//2
    solution = Vector(dim)
    if leaf(node, tree): 
        solution.set(0, 1)
        return solution
    if (node, char1, char2) in memo: return memo[(node, char1, char2)]
    leftChild, rightChild = tree[node]
    leftLabels1 = labeling[node][char1]["left set"]
    rightLabels1 = labeling[node][char1]["right set"]
    leftLabels2 = labeling[node][char2]["left set"]
    rightLabels2 = labeling[node][char2]["right set"]
    leftSolution = Vector(dim)
    for c1 in leftLabels1:
        for c2 in leftLabels2:
            leftSolution = leftSolution + pdvHelper(leftChild, c1, c2, tree, labeling)
    rightSolution = Vector(dim)
    pairsConsidered = set()
    for c1 in rightLabels1:
        for c2 in rightLabels2:
            rightSolution = rightSolution + pdvHelper(rightChild, c1, c2, tree, labeling)
    solution = (leftSolution * rightSolution).shift(diff(char1, char2))

    if node == 0 and char1 == char2:
        for i in range(1, solution.dim):
            solution.v[i] = solution.v[i]//2

    memo[(node, char1, char2)] = solution
    return solution


def generate(numLeaves, numCharacters):
    """ Generates a random tree, list of characters, and a dictionary mapping
        leaves to characters. Takes the number of leaves and the number of 
        characters as input.
    """
    boundary = [0]  # List of nodes that can take children. Initially, only the root.
    tree = {}       # The tree is a dictionary that associates a list of two children with
                    # each internal node.
    nodeCounter = 1
    for _ in range(numLeaves-1):    
        nextNode = random.choice(boundary)
        tree[nextNode] = [nodeCounter, nodeCounter+1]
        boundary.remove(nextNode)
        boundary.extend([nodeCounter, nodeCounter+1])
        nodeCounter += 2
    for leaf in boundary:
        tree[leaf] = []
    leafSet = boundary
    characters = range(numCharacters)
    leaves = {}
    for leaf in leafSet:
        leaves[leaf] = random.choice(characters)
    return tree, list(characters), leaves

def verify(vector):
    """ Takes a pairwise distance vector as input and verifies that the total number of 
    pairs of solutions indicated by that vector is consistent. That is, that the {vector_0 choose 2}, 
    which is the number of pairs of solutions is equal to the total number of pairs at distance 1 or more."""
    numPairs1 = vector.get(0) * (vector.get(0) - 1) // 2
    numPairs2 = sum([vector.get(i) for i in range(1, vector.dim)])
    return numPairs1 == numPairs2

def printLabeling(labeling):
    for node in labeling:
        print(node, labeling[node])
        print()

def main(tree, characters, leaves, verbose):
    """ Main function for the small parsimony problem"""
    labeling = {}
    maxparsimony(0, tree, characters, leaves, labeling)
    optimum = min([labeling[0][c]["score"] for c in characters])
    numSolutions = sum([labeling[0][c]["num solutions"] for c in characters if labeling[0][c]["score"] == optimum])
    printLabeling(labeling)
    d = diameter(tree, labeling)
    if verbose:
        print("Optimal score: ", optimum)
        print("Number of optimal solutions: ", numSolutions)
        print("Diameter: ", d)
        vector = pdv(tree, labeling)
        print("PDV: ", vector)
        assert(verify(vector))
        plt.bar(x = list(range(0, d+1)), height = vector.v[:d+1])
        plt.show()
    return numSolutions

# Generate a random tree with 101 leaves labeled at random from a set of 20 characters.
# Returns the tree, characters, and leaves.
tree, characters, leaves = generate(101, 20)

# Computing the diameter, pdv and plot results.
main(tree, characters, leaves, True)


            

    
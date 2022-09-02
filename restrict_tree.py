#Enter name of the gene program is working on
gene = "A0A075B6I0"
geneFileName = gene + ".fasta"
geneTreeFileName = gene + ".raxml.bestTree"
outputSubTreeFileName = gene + "_SubTree.txt"

#Enter how many genes needed for sub-tree
restrictNum = 50

geneSeqFile = open(geneFileName, "r")

geneName = geneSeqFile.readline()
geneName = geneName[1: geneName.find('|')] + '_' + geneName[geneName.find('|') + 1: geneName.rfind('|')] + '_' +geneName[geneName.rfind('|') + 1: geneName.find(' ')]
geneSeqFile.close()

with open(geneTreeFileName, "r") as file:
    tree = file.read().rstrip()

startSeed = tree.find(geneName)
endSeed = startSeed

oldStart, oldEnd = startSeed, endSeed

neededClosedBracket = 1
neededOpenBracket = 1
numOfGene, oldNumOfGene = 0, 0
subTree = ''

while numOfGene < restrictNum:
    for index in range(endSeed + 1, len(tree)):
        letter = tree[index]

        if letter == ')':
            neededClosedBracket -= 1
            if neededClosedBracket == 0:
                endSeed = index
                break
        elif letter == '(':
            neededClosedBracket += 1
        elif letter == ':' and tree[index - 1] != ')':
            numOfGene += 1

    for index in range(startSeed - 1, -1, -1):
        letter = tree[index]

        if letter == '(':
            neededOpenBracket -= 1
            if neededOpenBracket == 0:
                startSeed = index
                break
        elif letter == ')':
            neededOpenBracket += 1
        elif letter == ':' and tree[index - 1] != ')':
            numOfGene += 1
    
    if numOfGene <= restrictNum:
        subTree = tree[startSeed: endSeed+1]
        oldStart, oldEnd = startSeed, endSeed
    else:
        break
    neededClosedBracket, neededOpenBracket = 1, 1

outputFile = open(outputSubTreeFileName, "w")
outputFile.write(subTree)
from ncbitaxonomy import NCBITaxa

ncbi = NCBITaxa()
#database.update_taxonomy_database()
mammalianId = ncbi.get_name_translator(['Mammalia'])['Mammalia'][0]

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
oldRemoveTaxas = []
removeTaxas = []

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
            taxID = int(tree[tree[: index].rfind('_') + 1: index])
            if mammalianId in ncbi.get_lineage(taxID):
                numOfGene += 1
                if numOfGene > restrictNum:
                    break
            elif taxID not in removeTaxas:
                removeTaxas.append(taxID)

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
            taxID = int(tree[tree[: index].rfind('_') + 1: index])
            if mammalianId in ncbi.get_lineage(taxID):
                numOfGene += 1
                if numOfGene > restrictNum:
                    break
            elif taxID not in removeTaxas:
                removeTaxas.append(taxID)
    
    if numOfGene <= restrictNum:
        subTree = tree[startSeed: endSeed+1]
        oldStart, oldEnd = startSeed, endSeed
        oldRemoveTaxas = removeTaxas.copy()
        if index == (len(tree) - 1):
            break
    else:
        break
    neededClosedBracket, neededOpenBracket = 1, 1

if len(oldRemoveTaxas) != 0:
    for taxID in oldRemoveTaxas:
        if mammalianId not in ncbi.get_lineage(taxID):
            while subTree.find(str(taxID)) != -1:
                commaend = subTree.find(',', subTree.find(str(taxID)))
                bracketend = subTree.find(')', subTree.find(str(taxID)))
                commastart = subTree[: subTree.find(str(taxID))].rfind(',')
                bracketstart = subTree[: subTree.find(str(taxID))].rfind('(')

                if commaend < bracketend and commaend != -1:
                    if commastart > bracketstart:
                        subTree = subTree[:commastart] + subTree[commaend:]
                    else:
                        subTree = subTree[: bracketstart + 1] + subTree[commaend + 1:]
                else:
                    if commastart > bracketstart:
                        subTree = subTree[:commastart] + subTree[bracketend:] #may cause problem?
                    else:
                        subTree = subTree[: bracketstart - 1] + subTree[bracketend + 1:]

outputFile = open(outputSubTreeFileName, "w")
outputFile.write(subTree)
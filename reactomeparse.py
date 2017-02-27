from collections import defaultdict
import csv
import re


reactomeID = set()

# Dictionaries for quick storage/lookup whether or not a node is a child of another.
parentsDict = defaultdict(list)
childrensDict = defaultdict(list)

# Sets for identifying existing relationships, so no duplicates are output.
existingChebi = set()
existingNcbi = set()
existingRels = set()

# Function for removal of unnecessary characters from the IDs.
def cleanID(str):
    str = re.sub('[-]','',str)
    return str

# General function for parsing nodes and basic xref relationships out of
# the reactome files.
def parseReactome(nodesIn, nodesOut, edgesOut, source):
    reactomeReader = csv.reader(nodesIn, delimiter="\t")
    for line in reactomeReader:
        if source == "reactome":
            currentID = cleanID(line[0])
            if not currentID in reactomeID:
                reactomeID.add(currentID)
                nodesOut.write(currentID + "|" + line[1] + "|||reactome\n")
        elif source == "chebi":
            currentID = cleanID(line[1])
            if not (line[0],currentID) in existingChebi:
                existingChebi.add((line[0],currentID))
                edgesOut.write("CHEBI:" + line[0] + "|xref|reactome|" + currentID + "\n")
                if not currentID in reactomeID:
                    print("chebi")
                    reactomeID.add(currentID)
                    nodesOut.write(currentID + "|" + line[3] + "|||reactome\n")
        elif source == "ncbi":
            currentID = cleanID(line[1])
            if not (line[0],currentID) in existingNcbi:
                existingNcbi.add((line[0],currentID))
                edgesOut.write("nih.nlm.ncbi.gene.id:" + line[0] + "|xref|reactome|" + currentID + "\n")
                if not currentID in reactomeID:
                    print("ncbi")
                    reactomeID.add(currentID)
                    nodesOut.write(currentID + "|" + line[3] + "|||reactome\n")
        elif source == "reactrel":
            parentID = cleanID(line[0])
            childID = cleanID(line[1])
            parentsDict[parentID].append(childID)
            existingRels.add((parentID,childID))
            childrensDict[childID].append(parentID)
            edgesOut.write(childID + "|is_a|reactome|" + parentID + "\n")

# Function recursively seeks additional parents until it can find no more.
def isParentAChild(parent, children, edgesOut):
    if parent in childrensDict:
        for grandparent in childrensDict[parent]:
            if not (grandparent,parent) in existingRels:
                edgesOut.write(parent + "|is_a|reactome|" + grandparent + "\n")
                existingRels.add((grandparent,parent))
            if type(children) is list:
                for child in children:
                    if not (grandparent,child) in existingRels:
                        existingRels.add((grandparent,child))
                        edgesOut.write(child + "|is_a|reactome|" + grandparent + "\n")
                    isParentAChild(grandparent, child, edgesOut)
            elif type(children) is str:
                if not (grandparent,children) in existingRels:
                    existingRels.add((grandparent,children))
                    edgesOut.write(children + "|is_a|reactome|" + grandparent + "\n")
                isParentAChild(grandparent, children, edgesOut)
        return True
        # If top parent is found, end recursion
    else:
        return

def recursiveParentFinder(edgesOut):
    for parent, child in parentsDict.items():
        isParentAChild(parent, child, edgesOut)

def main():
    chebiIn = open("ChEBI2Reactome_All_Levels.txt","r")
    ncbiIn = open("NCBI2Reactome_All_Levels.txt","r")
    reactomePathIn = open("ReactomePathways.txt","r")
    reactomeRelsIn = open("ReactomePathwaysRelation.txt","r")

    nodesOut = open("reactomeNodesOut.csv","w")
    edgesOut = open("reactomeEdgesOut.csv","w")
    edgesXrefOut = open("reactomeEdgesOut.xref.csv","w")

    # write output file headers
    nodesOut.write("source_id:ID|name:string|synonyms:string[]|definition:string|:LABEL\n")
    edgesXrefOut.write(":START_ID|:TYPE|source:string|:END_ID\n")
    edgesOut.write(":START_ID|:TYPE|source:string|:END_ID\n")

    # parse reactome nodes from both chebi and ncbi gene files
    parseReactome(reactomePathIn, nodesOut, edgesXrefOut, "reactome")
    parseReactome(chebiIn, nodesOut, edgesXrefOut, "chebi")
    parseReactome(ncbiIn, nodesOut, edgesXrefOut, "ncbi")
    parseReactome(reactomeRelsIn, nodesOut, edgesOut, "reactrel")

    print(len(existingRels))
    recursiveParentFinder(edgesOut)
    print(len(existingRels))


    chebiIn.close()
    ncbiIn.close()
    reactomePathIn.close()
    reactomeRelsIn.close()
    nodesOut.close()
    edgesOut.close()

if __name__ == "__main__":
    main()

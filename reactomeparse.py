from collections import defaultdict
import csv
import re

reactomeID = set()
# parents = set()
# children = set()
parentsDict = defaultdict(list)
childrensDict = defaultdict(list)
existingRels = set()

def cleanID(str):
    str = re.sub('[-]','',str)
    return str

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
            edgesOut.write("CHEBI:" + line[0] + "|xref|reactome|" + currentID + "\n")
            if not currentID in reactomeID:
                print("chebi")
                reactomeID.add(currentID)
                nodesOut.write(currentID + "|" + line[3] + "|||reactome\n")
        elif source == "ncbi":
            currentID = cleanID(line[1])
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
            edgesOut.write("reactome:" + childID + "|is_a|reactome|" + parentID + "\n")

def isParentAChild(parent, children, edgesOut):
    if parent in childrensDict:
        for grandparent in childrensDict[parent]:
            # print(type(children))
            if not (grandparent,parent) in existingRels:
                # print(grandparent + " " + parent)
                edgesOut.write("reactome:" + parent + "|is_a|reactome|" + grandparent + "\n")
                existingRels.add((grandparent,parent))
            if type(children) is list:
                for child in children:
                    if not (grandparent,child) in existingRels:
                        existingRels.add((grandparent,child))
                        print(grandparent + " " + child)
                        edgesOut.write("reactome:" + child + "|is_a|reactome|" + grandparent + "\n")
                    isParentAChild(grandparent, child, edgesOut)
            elif type(children) is str:
                if not (grandparent,children) in existingRels:
                    existingRels.add((grandparent,child))
                    print(grandparent + " " + children)
                    edgesOut.write("reactome:" + children + "|is_a|reactome|" + grandparent + "\n")
                isParentAChild(grandparent, children, edgesOut)
        return True
    else:
        return

def recursiveParentFinder(edgesOut):
    for parent, child in parentsDict.items():
        if isParentAChild(parent, child, edgesOut) is True:
            pass
#
# def recursion(newparent, newchild, youngestchild):
#     if

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
    edgesOut.write(":START_ID|:TYPE|source:string|:END_ID\n")

    # parse reactome nodes from both chebi and ncbi gene files
    parseReactome(reactomePathIn, nodesOut, edgesXrefOut, "reactome")
    parseReactome(chebiIn, nodesOut, edgesXrefOut, "chebi")
    parseReactome(ncbiIn, nodesOut, edgesXrefOut, "ncbi")
    parseReactome(reactomeRelsIn, nodesOut, edgesOut, "reactrel")

    # print(existingRels)
    # print(('RMMU110373', 'RMMU110362') in existingRels)
    # for x,y in childrensDict.items():
    #     print(x)
    #     print(str(y))
    recursiveParentFinder(edgesOut)

    #     if x in childrensDict:
    #         for w in y:
    #             for z in childrensDict[x]:
    #                 edgesOut.write("reactome:" + w + "|is_a|reactome|" + z + "\n")

    # close dem shits
    chebiIn.close()
    ncbiIn.close()
    reactomePathIn.close()
    reactomeRelsIn.close()
    nodesOut.close()
    edgesOut.close()

if __name__ == "__main__":
    main()


# nih.nlm.ncbi.gene.id:
# The "Pathway hierarchy relationship" file consists of two columns of Reactome Stable identifiers (ST_ID), defining the relationship between pathways within the pathway hierarchy. The first column provides the parent pathway stable identifier,
# whereas the second column provides the child pathway stable identifier.

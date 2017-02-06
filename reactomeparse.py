from collections import defaultdict
import csv
import re

reactomeID = set()

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
                nodesOut.write(currentID + "|" + line[2] + "|||reactome\n")
        elif source == "chebi":
            currentID = cleanID(line[1])
            edgesOut.write("CHEBI:" + line[0] + "|part_of|reactome|" + currentID + "\n")
            if not currentID in reactomeID:
                print("chebi")
                reactomeID.add(currentID)
                nodesOut.write(currentID + "|" + line[3] + "|||reactome\n")
        elif source == "ncbi":
            currentID = cleanID(line[1])
            edgesOut.write("nih.nlm.ncbi.gene.id:" + line[0] + "|part_of|reactome|" + currentID + "\n")
            if not currentID in reactomeID:
                print("ncbi")
                reactomeID.add(currentID)
                nodesOut.write(currentID + "|" + line[3] + "|||reactome\n")



def main():
    chebiIn = open("ChEBI2Reactome_All_Levels.txt","r")
    ncbiIn = open("NCBI2Reactome_All_Levels.txt","r")
    reactomePathIn = open("ReactomePathways.txt","r")
    reactomeRelsIn = open("ReactomePathwaysRelation.txt","r")

    nodesOut = open("reactomeNodesOut.csv","w")
    edgesOut = open("reactomeEdgesOut.csv","w")

    # write output file headers
    nodesOut.write("source_id:ID|name:string|synonyms:string[]|definition:string|:LABEL\n")
    edgesOut.write(":START_ID|:TYPE|source:string|:END_ID\n")

    # parse reactome nodes from both chebi and ncbi gene files
    parseReactome(reactomePathIn, nodesOut, edgesOut, "reactome")
    parseReactome(chebiIn, nodesOut, edgesOut, "chebi")
    parseReactome(ncbiIn, nodesOut, edgesOut, "ncbi")

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

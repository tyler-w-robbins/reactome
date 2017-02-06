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
        currentID = cleanID(line[1])
        if not currentID in reactomeID:
            reactomeID.add(cleanID(line[1]))
            nodesOut.write(cleanID(line[1]) + "|" + line[3] + "|||reactome\n")
        if source == "chebi":
            edgesOut.write("CHEBI:" + line[0] + "|part_of|reactome|" + currentID + "\n")
        else:
            edgesOut.write("nih.nlm.ncbi.gene.id:" + line[0] + "|part_of|reactome|" + currentID + "\n")

def main():
    chebiIn = open("ChEBI2Reactome_All_Levels.txt","r")
    ncbiIn = open("NCBI2Reactome_All_Levels.txt","r")
    reactomeIn = open("ReactomePathwaysRelation.txt","r")

    nodesOut = open("reactomeNodesOut.csv","w")
    edgesOut = open("reactomeEdgesOut.csv","w")

    # write output file headers
    nodesOut.write("source_id:ID|name:string|synonyms:string[]|definition:string|:LABEL\n")
    edgesOut.write(":START_ID|:TYPE|source:string|:END_ID\n")

    # parse reactome nodes from both chebi and ncbi gene files
    parseReactome(chebiIn, nodesOut, edgesOut, "chebi")
    parseReactome(ncbiIn, nodesOut, edgesOut, "ncbi")

    # close dem shits
    chebiIn.close()
    ncbiIn.close()
    reactomeIn.close()
    nodesOut.close()
    edgesOut.close()

if __name__ == "__main__":
    main()


# nih.nlm.ncbi.gene.id:
# The "Pathway hierarchy relationship" file consists of two columns of Reactome Stable identifiers (ST_ID), defining the relationship between pathways within the pathway hierarchy. The first column provides the parent pathway stable identifier,
# whereas the second column provides the child pathway stable identifier.

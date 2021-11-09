"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file1 = open(filename, 'r')
    lines = file1.readlines()
    dna_str = ""
    for line in lines:
        dna_str += line.strip()
    file1.close()
    # print(dna_str)
    return dna_str


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    dna = dna.replace("T", "U")
    rna_lst = []
    stop_lst = ["UAA", "UAG", "UGA"]
    for i in range(startIndex, len(dna), 3):
        rna_lst.append(dna[i:i+3])
        if dna[i:i+3] in stop_lst:
            break
    # print("rna_lst=", rna_lst)
    return rna_lst


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f = open(filename)
    data_dict = json.load(f)
    codon_amino_dict = {}
    for key, value in data_dict.items():
        for i in range(len(value)):
            codon = value[i].replace("T", "U")
            codon_amino_dict[codon] = key
    f.close()
    # print(codon_amino_dict)
    return codon_amino_dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein_lst = []
    stop_lst = ["UAA", "UAG", "UGA"]
    for i in range(len(codons)):
        if codons[i] == "AUG" and len(protein_lst)==0:
            protein_lst.append("Start")
            if codons[i] in stop_lst:
                protein_lst.append("Stop")
        else:
            protein_lst.append(codonD[codons[i]])
    # print("protein_lst=", protein_lst)
    return protein_lst


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna_str = readFile(dnaFilename)
    codon_dict = makeCodonDictionary(codonFilename)
    protein_lst = []
    unused_base_count=0
    i = 0
    while i<len(dna_str):
        if(dna_str[i:i+3]=="ATG"):
            rna_lst = dnaToRna(dna_str, i)
            protiens_lst = generateProtein(rna_lst, codon_dict)
            protein_lst.append(protiens_lst)
            # i+=len(protein_lst[-1])*3
            i = i+(len(protein_lst[-1])*3)
        else:
            i = i+1
            unused_base_count+=1
    # print("unused_base_count=", unused_base_count)
    return protein_lst


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common_proteins_lst = []
    for list1 in proteinList1:
        for list2 in proteinList2:
            if list1 == list2:
                if list1 not in common_proteins_lst:
                    common_proteins_lst.append(list1)
    # print("common_proteins_lst==", common_proteins_lst)
    return common_proteins_lst
'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    return


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    return


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    return


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """
    test.testCommonProteins()
    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """

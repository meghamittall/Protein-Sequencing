"""
Protein Sequencing Project
Name:
Roll Number:
"""

from typing import List
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
    combine_proteins_lst = []
    for list in proteinList:
        for i in range(len(list)):
            combine_proteins_lst.append(list[i])
    # print("combine_proteins_lst==", combine_proteins_lst)
    return combine_proteins_lst


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    amino_dict = {}
    for amino in aaList:
        if amino not in amino_dict:
            amino_dict[amino] = 0
        amino_dict[amino] += 1
    # print("amino's dictionary==", amino_dict)
    return amino_dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    comb_p1, comb_p2 = combineProteins(proteinList1), combineProteins(proteinList2)
    amino_dct1, amino_dct2 = aminoAcidDictionary(comb_p1), aminoAcidDictionary(comb_p2)
    total1, total2= len(comb_p1), len(comb_p2)
    combined_list = comb_p2+comb_p1
    x=set(combined_list)
    l=list(x)
    final_list=[]
    for i in l:
        if (i!="Start" and i!="Stop"):
            if i in amino_dct1.keys():
                freq1=amino_dct1[i]/total1
            else:
                freq1 =0
            if i in amino_dct2.keys():
                freq2=amino_dct2[i]/total2
            else:
                freq2=0
            if abs(freq2-freq1)>=cutoff:
                temp=[]
                temp.append(i)
                temp.append(freq1)
                temp.append(freq2)
                final_list.append(temp)
    print("final_list==", len(final_list))
    return final_list


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("given list==", differences)
    temp = []
    difference = []
    for i in range(len(commonalities)):
         for j in range(len(commonalities[i])):
            if (commonalities[i][j] != "Start" and commonalities[i][j] != "Stop"):
                temp.append(commonalities[i][j])
    temp1=sorted(temp)
    print("The following proteins occured in both DNA sequences:")
    print(*temp1, sep = "\n")
    for i in range(len(differences)):
        difference.append(differences[i][0]+":")
        difference.append(str(round(differences[i][1]*100, 2))+" in Seq1")
        difference.append(str(round(differences[i][2]*100, 2))+" in Seq2")
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    print(*difference, sep = "\n")

    


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
    combine_p_lst1, combine_p_lst2 = combineProteins(proteinList1), combineProteins(proteinList2)
    dict1, dict2 =aminoAcidDictionary(combine_p_lst1), aminoAcidDictionary(combine_p_lst2) 
    labels_lst = []
    for key in dict1:
        labels_lst.append(key)
    for key in dict2:
        if key not in labels_lst:
            labels_lst.append(key)
    # print("aminoacidslist===", sorted(labels_lst))
    return sorted(labels_lst)


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    frequency_lst = []
    lst = combineProteins(proteinList)
    dct = aminoAcidDictionary(lst)
    for i in range(len(labels)):
        if labels[i] in dct.keys():
            freq = dct[labels[i]]/len(lst)
            frequency_lst.append(freq)
        else:
            frequency_lst.append(0)
    # print("frequency_lst===", frequency_lst)
    return frequency_lst


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    X_axis = np.arange(len(xLabels))
    plt.bar(X_axis, freqList1,width = -0.4, align="edge", label = label1)
    plt.bar(X_axis, freqList2, width = 0.4, align="edge", label = label2)
    plt.xticks(ticks=list(range(len(xLabels))), labels = xLabels)
    plt.title("Graph of Amino acids and their freq")
    plt.legend()
    plt.show()
    return None


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    import numpy as np
    ini_array1 = np.array(biggestDiffs)
    result = ini_array1.flatten()
    labels_color_lst = []
    # print("length of labels==", len(labels))
    for i in range(len(labels)):
        if labels[i] in result:
            labels_color_lst.append("black")
        else:
            labels_color_lst.append("white")
    # print("colors of labels===", labels_color_lst)
    # print("length of labels_color_lst==", len(labels_color_lst))
    return labels_color_lst


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
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    
    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testFindAminoAcidDifferences()
    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
    # test.testMakeAminoAcidLabels()
    test.testMakeEdgeList()

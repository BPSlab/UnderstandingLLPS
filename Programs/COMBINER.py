# Import Libraries

from pathlib import Path
import pandas
import os
from itertools import permutations
from itertools import combinations
import math
# User Defined Variables

# Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

# Define Output Directory
outputFilePath = (Path(mainProjectDirectory)/'Results'/'Peptide Design')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)

# Data Paths
proteinDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'), 'Protein')
dprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'), 'DPR')

# Motif Data Paths

motifDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/"Motif Picks.xlsx"))  # MOTIFS

#Function Reads from a motif dataframe and outputs a list of motifs

def motifListMaker(motif_dataframe):

    motifList = []

    for motifDataframeLine in range(0, motif_dataframe.shape[0]):
    
        motif = (motif_dataframe.iloc[motifDataframeLine, 0])
    
        motifList.append(motif)

    return(motifList)

#Function Reads from a 3 motif list and outputs a list of motif combinations

def motifCombinator(motif_list):

    motifCombinationList = []

    motifCombinations = combinations(motif_list, 3)

    for motifCombination in motifCombinations:

        motifCombinationList.append(list(motifCombination))

    return(motifCombinationList)

#Given 2 vectors and respective adjusters and a max Vetor outputs a score as detailed in User Guide/Motif Coexistence Score.docx

def vectorScorer(vector_1, vector_2, lengthScoreAdjuster, slopeScoreAdjuster, max_value):

    vectorLength = math.sqrt((vector_1**2)+(vector_2**2))
    vectorLengthScore = (vectorLength/max_value)*100

    try:
        vectorSlope = vector_2/vector_1
    except:
        vectorSlope = 0

    vectorAngle = math.degrees(math.atan(vectorSlope))
    vectorSlopeScore = 100-(((math.fabs(vectorAngle-45))/45)*100)

    vectorScore = (vectorLengthScore*lengthScoreAdjuster) + \
        (vectorSlopeScore*slopeScoreAdjuster)

    return(vectorScore)

#Finds max Vetor to be used in vectorScorer

def maxFinder(motif_combination_list, input_dataframe, control_variable,motif_list):

    maxValue = 0
    index = 1

    for motifCombination in motif_combination_list:

        index += 1

        motif1 = motifCombination[0]
        motif2 = motifCombination[1]
        motif3 = motifCombination[2]

        totalMotif13Interactions = 0
        totalMotif23Interactions = 0
        totalMotif21Interactions = 0

        for input_dataframeLineIndex in range(input_dataframe.shape[0]):

            sequence = str(input_dataframe.loc[input_dataframeLineIndex,control_variable])

            motif1Count = motifCounter(motif1, sequence,motif_list)
            motif2Count = motifCounter(motif2, sequence,motif_list)
            motif3Count = motifCounter(motif3, sequence,motif_list)

            if motif1Count > 0 and motif3Count > 0 and (motif1 != motif3):
                totalMotif13Interactions = totalMotif13Interactions + 1
            if motif1Count > 1 and motif3Count > 1 and (motif1 == motif3):
                totalMotif13Interactions = totalMotif13Interactions + 1

            if motif2Count > 0 and motif3Count > 0 and (motif2 != motif3):
                totalMotif23Interactions = totalMotif23Interactions + 1
            if motif2Count > 1 and motif3Count > 1 and (motif2 == motif3):
                totalMotif23Interactions = totalMotif23Interactions + 1

            if motif2Count > 0 and motif1Count > 0 and (motif2 != motif1):
                totalMotif21Interactions = totalMotif21Interactions + 1
            if motif2Count > 1 and motif1Count > 1 and (motif2 == motif1):
                totalMotif21Interactions = totalMotif21Interactions + 1

        vectorLength1 = math.sqrt(
            (totalMotif13Interactions**2)+(totalMotif21Interactions**2))
        vectorLength2 = math.sqrt(
            (totalMotif21Interactions**2)+(totalMotif23Interactions**2))
        vectorLength3 = math.sqrt(
            (totalMotif13Interactions**2)+(totalMotif23Interactions**2))

        vectorLengthList = [vectorLength1, vectorLength2, vectorLength3]

        if max(vectorLengthList) > maxValue:

            maxValue = max(vectorLengthList)
            maxMotifCombination = motifCombination
    
    print(maxMotifCombination,maxValue)
    return(maxValue)


#Counts distinct motif apearances in the sequence based on the list of motifs

def motifCounter(imput_motif,input_sequence,motif_list):

    principalMotif = imput_motif

    principalMotifCount = 0

    principalMotifSize = len(principalMotif)

    secundaryMotifList = motif_list.copy()

    secundaryMotifList.remove(principalMotif)

    for sequenceIndex in range(0, (len(input_sequence)-(principalMotifSize-1))):

        iteratingMotif = ''

        for aminoAcidIndex in range(0, (principalMotifSize)):

            aminoAcid = str(input_sequence[sequenceIndex + aminoAcidIndex])

            iteratingMotif += aminoAcid

        if iteratingMotif == principalMotif:

            isMotifUnique = True

            for secundaryMotif in secundaryMotifList:

                secundaryMotifSize = len(secundaryMotif)

                motifLengthDiference = abs(
                    principalMotifSize-secundaryMotifSize)

                if motifLengthDiference != 0:

                    for operationIndex in range(motifLengthDiference+1):

                        splitLowerBound = sequenceIndex - \
                            ((motifLengthDiference)-operationIndex)

                        splitUpperBound = sequenceIndex + principalMotifSize + operationIndex

                        splittedSequence = input_sequence[splitLowerBound:splitUpperBound]

                        if splittedSequence == secundaryMotif:

                            isMotifUnique = False

            if isMotifUnique:

                principalMotifCount += 1

    return(principalMotifCount)


#Finds overlap between two motifs and joins them by it, returning the resulting string

def motifJoiner(motif_1, motif_2):

    listElement1 = list(motif_1)
    listElement2 = list(motif_2)

    minLength = min(len(listElement1), len(listElement2))

    listOfFirst = []
    listOfOverLap = []
    listOfLast = []

    for index in range(len(listElement1)-minLength, len(listElement1)):

        if (listElement1[index:]) == (listElement2[:(len(listElement1[index:]))]):

            listOfOverLap.extend(listElement1[index:])

            listOfFirst.extend(
                listElement1[:(len(listElement1)-len(listElement1[index:]))])

            listOfLast.extend(
                listElement2[(len(listElement2)-(len(listElement2)-len(listOfOverLap))):])

            break

    listOfFinalMotif = []
    listOfFinalMotif.extend(listOfFirst)
    listOfFinalMotif.extend(listOfOverLap)
    listOfFinalMotif.extend(listOfLast)

    if len(listOfFinalMotif) == 0:

        listOfFinalMotif.extend(listElement1)
        listOfFinalMotif.extend(listElement2)

    stringOfFinalMotif = ''

    for peptide in listOfFinalMotif:

        stringOfFinalMotif += peptide

    return(stringOfFinalMotif)

#Given a list of motifs runs motifJoiner for motif permutations of the list, returns the shortest of the permutations

def tripleMotifJoiner(motif_list):

    lengthList = []
    peptideList = []

    for permutation in permutations(motif_list):

        initialPeptide = ''

        permutationList = list(permutation)

        while len(permutationList) > 0:

            initialPeptide = motifJoiner(initialPeptide, permutationList[0])

            del permutationList[0]

        peptideList.append(initialPeptide)
        lengthList.append(len(initialPeptide))

    minimumLength = max(lengthList)
    minimumPeptide = ''

    for listIndex in range(len(lengthList)):

        length = lengthList[listIndex]
        peptide = peptideList[listIndex]

        if length <= minimumLength:

            minimumPeptide = peptide
            minimumLength = length

    return(minimumPeptide, minimumLength)


"""
Given the list of motif combinations and a input dataframe, as well as the Score adjustment Values,
finds the Score as per Motif Coexistence Score.docx and the shortest Peptide formed by the motif list, outputs a Data frame.

"""
def outputMaker(motif_combination_list, input_dataframe, lengthScoreAdjuster, slopeScoreAdjuster, control_variable,motif_list):
    #Create empty dataframe with these column names
    outputDataframecolumnList = ['Motif List','Score','Shortest Peptide']
    outputDataframe = pandas.DataFrame(columns = outputDataframecolumnList)
    outputDataframeLineIndex = 0
    maxValue = maxFinder(motif_combination_list, input_dataframe, control_variable,motif_list)
    
    for motifCombination in motif_combination_list:

        motif1 = motifCombination[0]
        motif2 = motifCombination[1]
        motif3 = motifCombination[2]

        totalMotif13Interactions = 0
        totalMotif23Interactions = 0
        totalMotif21Interactions = 0

        for input_dataframeLineIndex in range(input_dataframe.shape[0]):

            sequence = str(input_dataframe.loc[input_dataframeLineIndex,
                                           control_variable])

            motif1Count = motifCounter(motif1, sequence,motif_list)
            motif2Count = motifCounter(motif2, sequence,motif_list)
            motif3Count = motifCounter(motif3, sequence,motif_list)

            if motif1Count > 0 and motif3Count > 0 and (motif1 != motif3):
                totalMotif13Interactions = totalMotif13Interactions + 1
            if motif1Count > 1 and motif3Count > 1 and (motif1 == motif3):
                totalMotif13Interactions = totalMotif13Interactions + 1

            if motif2Count > 0 and motif3Count > 0 and (motif2 != motif3):
                totalMotif23Interactions = totalMotif23Interactions + 1
            if motif2Count > 1 and motif3Count > 1 and (motif2 == motif3):
                totalMotif23Interactions = totalMotif23Interactions + 1

            if motif2Count > 0 and motif1Count > 0 and (motif2 != motif1):
                totalMotif21Interactions = totalMotif21Interactions + 1
            if motif2Count > 1 and motif1Count > 1 and (motif2 == motif1):
                totalMotif21Interactions = totalMotif21Interactions + 1

        vectorScore13 = vectorScorer(totalMotif21Interactions, totalMotif23Interactions, lengthScoreAdjuster, slopeScoreAdjuster, maxValue)
        vectorScore23 = vectorScorer(totalMotif21Interactions, totalMotif13Interactions, lengthScoreAdjuster, slopeScoreAdjuster, maxValue)
        vectorScore12 = vectorScorer(totalMotif13Interactions, totalMotif23Interactions, lengthScoreAdjuster, slopeScoreAdjuster, maxValue)
        
        finalScore = (vectorScore13*(1/3))+(vectorScore23*(1/3))+(vectorScore12*(1/3))
        
        outputDataframe.at[outputDataframeLineIndex,'Motif List'] = motifCombination
        outputDataframe.at[outputDataframeLineIndex,'Score'] = finalScore
        outputDataframe.at[outputDataframeLineIndex,'Shortest Peptide'] = tripleMotifJoiner(motifCombination)
        
        outputDataframeLineIndex += 1
        
    return(outputDataframe)


"""
Given input datframe , motif dataframe and control variable  runs outputMaker on those parameters and saves its 
output to a output_filename named excel file.

"""
def main(input_dataframe,motif_dataframe,control_variable,output_filename):    
    
    motifList = motifListMaker(motif_dataframe)
    motifCombinationList = motifCombinator(motifList)

    
    outputDataframe = outputMaker(motifCombinationList, input_dataframe, 0.5, 0.5, control_variable,motifList)
    
    outputPathFileName = (Path(outputFilePath)/output_filename)
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:

        outputDataframe.to_excel(writer, sheet_name='Peptide', index=True)

main(dprDataframe,motifDataframe,"DPR",'Motif Combinations.xlsx')        





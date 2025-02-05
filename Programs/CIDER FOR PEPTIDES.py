#Import Libraries
from pathlib import Path
import pandas
import os
from itertools import permutations
from localcider.sequenceParameters import SequenceParameters



# Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

outputFilePath = (Path(mainProjectDirectory)/'Results'/'Cider for peptides')

if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)

#Read Databases from Excel sheets     
peptideDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Peptide Design'/'Motif Combinations.xlsx'))

def ciderAnalisis(sequence):
    
    SeqObject = SequenceParameters(sequence)
    
    fractionOfDisorder = SeqObject.get_fraction_disorder_promoting()
    ncpr = SeqObject.get_NCPR(pH=None)
    isoelectricPoint = SeqObject.get_isoelectric_point()
    fcr = SeqObject.get_FCR()
    omega = SeqObject.get_Omega()
    ppiiPropensity = SeqObject.get_PPII_propensity()
    scd = SeqObject.get_SCD()
    wwHydropathy = SeqObject.get_WW_hydropathy()
    fractionExpanding = SeqObject.get_fraction_expanding()
    fractionNegative = SeqObject.get_fraction_negative()
    fractionPositive = SeqObject.get_fraction_positive()
    kappa = SeqObject.get_kappa()
    length = SeqObject.get_length()
    if len(sequence) >= 5 :
        linearSigma = SeqObject.get_linear_sigma()
    else:
        linearSigma = None
    meanHydrophaty = SeqObject.get_mean_hydropathy()
    netCharge = SeqObject.get_mean_net_charge()
    
    return([ncpr,
            isoelectricPoint,
            fcr,
            omega,
            ppiiPropensity,
            scd,
            wwHydropathy,
            fractionExpanding,
            fractionNegative,
            fractionPositive,
            kappa,
            length,
            linearSigma,
            meanHydrophaty,
            netCharge,
            fractionOfDisorder])

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

#For a Motif List calls the motifJoiner to join them with the overlap removed

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
        
    return(peptideList)

#Cleans, Filters and Resets the index of a Data Frame

def dataFrameFormater(input_dataframe):

    input_dataframe.drop(columns='Shortest Peptide', inplace=True)
    input_dataframe.drop(columns='Unnamed: 0', inplace=True)
    input_dataframe.drop(input_dataframe[input_dataframe['Score'] < 10].index, inplace = True)
    input_dataframe.reset_index(drop = True,inplace = True )
    return(input_dataframe)

'''
For each Motif list calculates the possible peptides, calls ciderAnalisis for
each of those and, together with the Score writes to outputDataframe and saves to output_filename
'''

def outputMaker(input_dataframe,output_filename):
    
    columnNames = ["Motif List",
                   "Score",
                   "Peptide",
                   
                   "NCPR",
                   "Isoelectric Point", 
                   "FCR",
                   "Omega",
                   "PPII Propensity",
                   "SCD",
                   "Hydropathy",
                   "Fraction Expanding",
                   "Fraction Negative",
                   "Fraction Positive",
                   "Kappa",
                   "Length",
                   "linear Sigma",
                   "Mean Hydropathy",
                   "Mean Net Charge",
                   "Fraction Of Disorder Pormoting Residues"]
    
    outputDataframe = pandas.DataFrame(columns = columnNames)
    
    
    formatedDataframe = dataFrameFormater(input_dataframe)
    
    outputDataframeVerticalIndex = 0
    
    for formatedDataframeLineIndex in range(formatedDataframe.shape[0]):
        
        print(formatedDataframeLineIndex+1,'/',formatedDataframe.shape[0])
        
        motifList = (formatedDataframe.loc[formatedDataframeLineIndex,'Motif List']).replace('[','').replace(']','').replace(' ','').replace("'","").split(',')
        motifListScore = formatedDataframe.loc[formatedDataframeLineIndex,'Score']

        peptideList = tripleMotifJoiner(motifList)
        
        
        for peptide in peptideList:
            
            sequenceCiderResults = ciderAnalisis(str(peptide))          
            
            outputDataframe.at[outputDataframeVerticalIndex,"Motif List"] = motifList
            outputDataframe.at[outputDataframeVerticalIndex,"Score"] = motifListScore
            outputDataframe.at[outputDataframeVerticalIndex,"Peptide"] = peptide

            outputDataframe.at[outputDataframeVerticalIndex,"NCPR"] = sequenceCiderResults[0]
            outputDataframe.at[outputDataframeVerticalIndex,"Isoelectric Point"] = sequenceCiderResults[1]
            outputDataframe.at[outputDataframeVerticalIndex,"FCR"] = sequenceCiderResults[2]
            outputDataframe.at[outputDataframeVerticalIndex,"Omega"] = sequenceCiderResults[3]
            outputDataframe.at[outputDataframeVerticalIndex,"PPII Propensity"] = sequenceCiderResults[4]
            outputDataframe.at[outputDataframeVerticalIndex,"SCD"] = sequenceCiderResults[5]
            outputDataframe.at[outputDataframeVerticalIndex,"Hydropathy"] = sequenceCiderResults[6]
            outputDataframe.at[outputDataframeVerticalIndex,"Fraction Expanding"] = sequenceCiderResults[7]
            outputDataframe.at[outputDataframeVerticalIndex,"Fraction Negative"] = sequenceCiderResults[8]
            outputDataframe.at[outputDataframeVerticalIndex,"Fraction Positive"] = sequenceCiderResults[9]
            outputDataframe.at[outputDataframeVerticalIndex,"Kappa"] = sequenceCiderResults[10]
            outputDataframe.at[outputDataframeVerticalIndex,"Length"] = sequenceCiderResults[11]
            outputDataframe.at[outputDataframeVerticalIndex,"linear Sigma"] = sequenceCiderResults[12]
            outputDataframe.at[outputDataframeVerticalIndex,"Mean Hydropathy"] = sequenceCiderResults[13]
            outputDataframe.at[outputDataframeVerticalIndex,"Mean Net Charge"] = sequenceCiderResults[14]
            outputDataframe.at[outputDataframeVerticalIndex,"Fraction Of Disorder Pormoting Residues"] = sequenceCiderResults[15]

            outputDataframeVerticalIndex += 1
            
    outputPathFileName = (Path(outputFilePath)/output_filename)
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:
    
        outputDataframe.to_excel(writer, sheet_name='DPR', index=True)


outputMaker(peptideDataframe,'Cider for peptides.xlsx')









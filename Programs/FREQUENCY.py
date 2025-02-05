#Import Libraries

from pathlib import Path
import pandas
import os

#Define Main Project and ouput Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent
outputFilePath = (Path(mainProjectDirectory)/'Results'/'Motif Frequency')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)
    
#Make Dataframes from excel in Data Paths 
proteinDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'),'Protein')
dprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'),'DPR')
nodprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'),'NODPR')
negativeDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Negative Database.xlsx'))
rnabindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'RNA binding')
dnabindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'DNA binding')
chromatinbindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Chromatin binding')
regulationDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Regulation')
hydrolaseDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Hydrolase')
structureDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Structure')
motifDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/"Motif Picks.xlsx"))


#User Defined Functions 
"""
Function that takes a motif and turns it into a full file name in the format /motif.xslx
It joins the file name with the Results file path for the full output file path.
It then initializes a dataframe variable named "resultDataframe" where the results of the 
computations detailed bellow will be stored. It initializes two list variables where one is 
the list of the Uniprot IDs of proteins that the motif appears in and the other is the list 
of the Uniprot IDs of proteins where the motif is absent. Next it cycles trought the lines 
of the DPR data frame. In each line it takes the Uniprot ID and DPR sequence, counts motif
presence. If presence is positive it is added to the list of the Uniprot IDs. Then writes to 
motifDataframe Unidprot ID,Sequence, Motif Frequency and Presence, as well sa and the number
of Uniprot IDs the motif is found/absent. Finally it writes motifDataframe to Motif excell with name motif_name.xlsx

""" 
   
def makeMotifExcel(motif_name,input_dataframe,output_path_name,control_variable):
    
    motifDataframeColumnNames = ["Uniprot ID", "Sequence", "Motif Frequency", "Motif Presence", "In Uniprot", "Not In Uniprot"]
    motifDataframe = pandas.DataFrame(columns = motifDataframeColumnNames)
    listOfUniIdsIn=[]
    listOfUniIdsNotIn=[]
    
    for dataframeLineIndex in range (0 , input_dataframe.shape[0]):
        
        countDpr = 0
        
        uniprotId = str(input_dataframe.loc[dataframeLineIndex,'UniprotID'+' '+control_variable])
        
        sequence = str(input_dataframe.loc[dataframeLineIndex,control_variable])
        
        countDprTotal = sequence.count(motif_name)
        
        if countDprTotal != 0:
            
            countDpr += 1
            
        if uniprotId not in listOfUniIdsIn and countDpr != 0:
            
            listOfUniIdsIn.append(uniprotId)
            
        if uniprotId not in listOfUniIdsNotIn:
            
            listOfUniIdsNotIn.append(uniprotId)
            
        motifDataframe.at[dataframeLineIndex,"Uniprot ID"] = uniprotId 
        motifDataframe.at[dataframeLineIndex,"Sequence"] = sequence 
        motifDataframe.at[dataframeLineIndex,"Motif Frequency"] = countDprTotal
        motifDataframe.at[dataframeLineIndex,"Motif Presence"] = countDpr
        
    motifDataframe.at[0,"In Uniprot"] = len(listOfUniIdsIn)
    motifDataframe.at[0,"Not In Uniprot"] = len(listOfUniIdsNotIn)-len(listOfUniIdsIn)
    
    
    outputFilePathComplete = (Path(outputFilePath)/output_path_name)
    outputFileName = motif_name + ".xlsx"

    outputPathFileName = (Path(outputFilePathComplete)/outputFileName)    
    if not os.path.exists(outputFilePathComplete):
        os.makedirs(outputFilePathComplete)
    #Write output_dataframe to excel file in outputPathFileName Path
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:
        motifDataframe.to_excel(writer, sheet_name='Motif Search', index=False)      
    

'''
Function takes input_dataframe, inicializes a list of motifs and 
reads lines of motifDataframe. It takes motifs from it and appends
them to the list of motifs. For each motif it calls
makeMotifExcel(motif,input_dataframe)

'''

def makeMotifList(input_dataframe,output_path_name,control_variable):
    motifList = []
    for lines in range (0 , motifDataframe.shape[0]):
        motif = motifDataframe.loc[lines,'Motifs']
        if motif not in motifList:
            motifList.append(motif)
    for motifs in motifList:
        makeMotifExcel(motifs,input_dataframe,output_path_name,control_variable)

    
# Use of Frequency function for data processing
# Change as needed to respetive Dataframe to be analised defined above

makeMotifList(dprDataframe,'DPR','DPR')
makeMotifList(rnabindingDataframe,'RNA binding','DPR')
makeMotifList(dnabindingDataframe,'DNA binding','DPR')
makeMotifList(chromatinbindingDataframe,'Chromatin binding','DPR')
makeMotifList(regulationDataframe,'Regulation','DPR')
makeMotifList(hydrolaseDataframe,'Hydrolase','DPR')
makeMotifList(structureDataframe,'Structure','DPR')
    
     

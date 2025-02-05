#Import Libraries
from pathlib import Path
import pandas 
import os


#Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent


#Define peptide size dictionary
peptideSizeDictionary = {1:'Amino Acid',
                         2:'Dipeptide',
                         3:'Tripeptide',
                         4:'Tetrapeptide',
                         5:'Pentapeptide',
                         6:'Hexapeptide'
                         }    

#Read from excel
negativeDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Negative Database.xlsx'))
dprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'),'DPR')
rnabindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'RNA binding')
dnabindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'DNA binding')
chromatinbindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Chromatin binding')
regulationDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Regulation')
hydrolaseDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Hydrolase')
structureDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database By Families.xlsx'),'Structure')

def selector(input_dataframe,peptide_size,control_variable):
    
    #Data Paths 
    outputFilePath = (Path(mainProjectDirectory)/'Results'/'Peptides')
    if not os.path.exists(outputFilePath):
        os.makedirs(outputFilePath)
        
    #Define columns
    columnNames = ["Motif","Presence","Frequency"]
    
    motifDataframe = pandas.DataFrame(columns = columnNames)
    
    
    #Get all the protein sequences together
    fullSequence = ""
    
    for input_dataframeLineIndex in range (input_dataframe.shape[0]):
        
        smallSequence = input_dataframe.loc[input_dataframeLineIndex,control_variable]
        
        if type(smallSequence) == str:
        
            fullSequence = fullSequence + '-' + smallSequence
        
    print("Full Sequence Ready")   
     
    #get patterns and analyse
    listOfMotifs = [] 
    
    for sequenceIndex in range(0,(len(fullSequence)-(peptide_size-1))):
        
        motif = ''
        
        for aminoAcidIndex in range (0,(peptide_size)):
        
            aminoAcid = str(fullSequence[sequenceIndex + aminoAcidIndex])
            
            motif += aminoAcid
            
        if (motif not in listOfMotifs) and ('-' not in  motif):
            listOfMotifs.append(motif)
            
    print("List of Motifs Ready")
           
    index = 1
    
    for motifs in listOfMotifs:
    
        print((index/(len(listOfMotifs)))*100)
        motifPresence = 0
        motifFrequency = 0
    
        
        for input_dataframeLineIndex in range(0,input_dataframe.shape[0]):
            
            sequence = input_dataframe.loc[input_dataframeLineIndex,control_variable]
            
            if type(sequence) == str:
            
                if motifs in sequence:
                    
                    motifPresence = motifPresence + 1
                    motifFrequency = motifFrequency + sequence.count(str(motifs))               
        
        motifDataframe.at[index,"Motif"] = motifs
        motifDataframe.at[index,"Presence"] = motifPresence
        motifDataframe.at[index,"Frequency"] = motifFrequency
           
        index = index + 1
    
    fileName = (str(peptideSizeDictionary[peptide_size])+'.xlsx')
    # Create a Pandas Excel writer using openpyxl as the engine.
    with pandas.ExcelWriter(os.path.join(outputFilePath, fileName), engine='openpyxl'  ) as writer:
        # Write the dataframe to the 'Data' sheet without index.
        motifDataframe.to_excel(writer, sheet_name=str(peptideSizeDictionary[peptide_size]), index=False)

#Run Selector for desired peptide lengths 
selector(structureDataframe,3,'DPR')
selector(structureDataframe,4,'DPR')
selector(structureDataframe,5,'DPR')
selector(structureDataframe,6,'DPR')


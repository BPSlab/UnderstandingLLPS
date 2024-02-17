#Import Libraries

from pathlib import Path
import pandas
import os


#User Defined Variables

#Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

outputFilePath = (Path(mainProjectDirectory)/'Results'/'Matrix')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)
    
#Data Paths 
dprStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'DPR Stats.xlsx'))
nodprStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'NO DPR Stats.xlsx'))
proteinStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'Protein Stats.xlsx')) #DATABASE

motifDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/"Tripeptide Picks.xlsx")) #MOTIFS


# Create the folder if it doesn't exist
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)

#Define column names
matriceDataFrameColumnNames = []
for motifDataframeLine in range (0,motifDataframe.shape[0]):
    motif = (motifDataframe.iloc[motifDataframeLine,0])
    matriceDataFrameColumnNames.append(motif)

'''
define matrician function that takes an input dataframe of lines and columns
named after the Motifs in motifDataframe. For each Sequence in the 
input_dataframe it counts interactions of the motifs with eachother in the sequence
and saves that value to the currespondent cell in the matrix. Finally it writes the 
dataframe to Motif Occurence Matrix.xlsx
'''
def matrician (input_dataframe):
    #Create empty dataframe with these column names
    matrice = pandas.DataFrame(columns = matriceDataFrameColumnNames,index = matriceDataFrameColumnNames)
    
    for matriceLines in range (0,matrice.shape[0]):
        motif1 = matrice.index[matriceLines]
        
        for matriceRows in range (0,matrice.shape[1]):
            motif2 = matrice.columns[matriceRows]
            
            totalMotifInteractions = 0
            for sequence in input_dataframe.loc[:,"Sequence"]:
                
                motif1Count = sequence.count(str(motif1))
                motif2Count = sequence.count(str(motif2))
                localMotifCount = (motif1Count + motif2Count)
                if motif1Count > 0 and motif2Count > 0 and (motif1 != motif2):           
                    totalMotifInteractions = totalMotifInteractions + 1
                if motif1Count > 1 and motif2Count > 1 and (motif1 == motif2):
                    totalMotifInteractions = totalMotifInteractions + 1
            matrice.at[str(motif1),str(motif2)] = totalMotifInteractions
            
    
    outputPathFileName = (Path(outputFilePath)/'Motif Occurrence Matrix.xlsx')
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:

        matrice.to_excel(writer, sheet_name='Matrix', index=True)

matrician(dprStatisticsDataframe)

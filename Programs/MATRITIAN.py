from pathlib import Path
import pandas
import os

# Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent


# Define and or Create outputFilePath

outputFilePath = (Path(mainProjectDirectory)/'Results'/'Converted Matrix')

if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)

# Read peptideDataframe and Create outputDataframe
    
peptideDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Peptide Design'/'Motif Combinations.xlsx'))

outputDataframe = pandas.DataFrame()


'''
For each line in peptideDataframe, reads Motif List and Motif Score, Then finds or creates
the corresponding cell in the Matrix and writees Motif Score to it, finaly writes dataframe to'MATRIX.xlsx.
'''

for peptideDataframeLines in range (0,peptideDataframe.shape[0]):
    
    peptideDataframeMotifListString = peptideDataframe.loc[peptideDataframeLines,'Motif List']
    peptideDataframeScore = peptideDataframe.loc[peptideDataframeLines,'Score']

    peptideDataframeMotifListString = peptideDataframeMotifListString.replace('[','').replace(']','').replace("'","").replace(' ','')
    peptideDataframeMotifList = peptideDataframeMotifListString.split(',')
    
    peptideDataframeMotifListFirstMotif = str(peptideDataframeMotifList[0])
    
    peptideDataframeMotifListSecondThrirdMotif = str(peptideDataframeMotifList[-2:])
    
    outputDataframe.at[peptideDataframeMotifListFirstMotif,peptideDataframeMotifListSecondThrirdMotif] = peptideDataframeScore

outputPathFileName = (Path(outputFilePath)/'MATRIX.xlsx')

# Create a Pandas Excel writer using XlsxWriter as the engine.
with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:

    outputDataframe.to_excel(writer, sheet_name='Matrix', index=True)

   
            
        
            

        
  

            
    
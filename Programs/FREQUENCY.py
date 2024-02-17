#Import Libraries

from pathlib import Path
import pandas
import os


#User Defined Variables

#Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

outputFilePath = (Path(mainProjectDirectory)/'Results'/'Motif Presence')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)
    
#Make Dataframes from excel in Data Paths 
dprStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'DPR Stats.xlsx'))
nodprStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'NO DPR Stats.xlsx'))
proteinStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'Protein Stats.xlsx'))
negativeStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'Negative Stats.xlsx'))
rnaBindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'RNA binding.xlsx'))
dnaStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'DNA binding.xlsx'))
regulationDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Regulation.xlsx'))
chromatinbindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Chromatin binding.xlsx'))
hydrolaseDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Hydrolase.xlsx'))
structureDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Structure.xlsx'))

motifDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/"Motif Picks.xlsx"))


#User Defined Functions 
"""
Function that takes a motif ,then turns it into a full file name in the 
format /motif.xslx after that it joins the file name with the Results file
path for the full output file path.
Then it initializes a dataframe variable named "resultDataframe" where the 
results of the computations detailed bellow will be stored, it then 
initializes two list variables one is the list of the Uniprot IDs of 
Proteins that the target sequence appears in,the other is the list of the 
Uniprot IDs of Proteins that the target sequence does not appear in, next it 
cycles trought the lines of the droplet promoting regions data frame in each line
takes uniId and dpr sequence , counts Motif apearences in that dpr and adds it
to the list of the Uniprot IDs if the motif is found in the sequence and to the 
list of the Uniprot IDs of Proteins that the target sequence appears in.
Then writes to motifDataframe Unidprot ID,Sequence, Motif Frequency and Presence
and the number of Uniprots the Motif is found in and the number of Uniprod Ids 
its not, finally it writes motifDataframe to Motif excell with name motif_name.xlsx


""" 
   
def makeMotifExcel(motif_name,target_sequences_dataframe):
    
    motifFileName ="/" + motif_name + ".xlsx"
    fullFileName = str(outputFilePath) + motifFileName
    print(fullFileName)
    
    motifDataframeColumnNames = ["Uniprot ID", "DPR", "Motif Frequency", "Motif Presence", "In Uniprot", "Not In Uniprot"]
    motifDataframe = pandas.DataFrame(columns = motifDataframeColumnNames)
    listOfUniIdsIn=[]
    listOfUniIdsNotIn=[]
    
    for lines in range (0 , target_sequences_dataframe.shape[0]):
        
        countDpr = 0
        uniId = target_sequences_dataframe.iloc[lines,0]
        dpr = target_sequences_dataframe.iloc[lines,1]
        countDprTotal = dpr.count(motif_name)
        
        if countDprTotal != 0:
            
            countDpr += 1
            
        if uniId not in listOfUniIdsIn and countDpr != 0:
            
            listOfUniIdsIn.append(uniId)
            
        if uniId not in listOfUniIdsNotIn:
            
            listOfUniIdsNotIn.append(uniId)
            
        motifDataframe.at[lines,"Uniprot ID"] = uniId 
        motifDataframe.at[lines,"DPR"] = dpr 
        motifDataframe.at[lines,"Motif Frequency"] = countDprTotal
        motifDataframe.at[lines,"Motif Presence"] = countDpr
        
    motifDataframe.at[0,"In Uniprot"] = len(listOfUniIdsIn)
    motifDataframe.at[0,"Not In Uniprot"] = len(listOfUniIdsNotIn)-len(listOfUniIdsIn)
    
    
    with pandas.ExcelWriter(fullFileName, engine='openpyxl') as writer:
        motifDataframe.to_excel(writer, sheet_name="1", index=False)
    print(motif_name)   


'''
Function takes target_sequences_dataframe inicializes a list of motifs and 
read lines of target_sequences_dataframe takes motifs from it and appends them
to the list of motifs , for each motif it calls
makeMotifExcel(motif,target_sequences_dataframe)

'''

def makeMotifList(target_sequences_dataframe):
    motifList = []
    for lines in range (0 , motifDataframe.shape[0]):
        motif = motifDataframe.iloc[lines,0]
        if motif not in motifList:
            motifList.append(motif)
    for motifs in motifList:
        makeMotifExcel(motifs,target_sequences_dataframe)
    print("END") 
    

makeMotifList(dprStatisticsDataframe)#Change to respetive Dataframe to be analised defined above
        
     

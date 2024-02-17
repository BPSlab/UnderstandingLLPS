#Import Libraries
from pathlib import Path
import pandas
import os

#Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

#Data Paths 
outputFilePath = (Path(mainProjectDirectory)/'Results'/'Statistics')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)

#Read Databases from Excel sheets 
dprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Protein Database.xlsx'),'DPR')
nodprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Protein Database.xlsx'),'NODPR')
negativeDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Negative Database.xlsx'))
familiesDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Families Data.xlsx'))

#Define columns
columnNamesList = ["Uniprot ID", "Sequence", "Total aa", "A", "%A", "G", "%G", "V", "%V", "L", "%L", "I", "%I", "P", "%P", "F", "%F", 
                "M", "%M", "W", "%W", "C", "%C", "D", "%D", "E", "%E", "R", "%R", "K", "%K", "H", "%H", "N", "%N", "Q", "%Q", "S", "%S", 
                "T", "%T", "Y", "%Y", "Hydrofobic", "% Hydrofobic", "Negative Charged", "% Negative Charged", "Positive Charged", 
                "% Positive Charged", "Polar", "% Polar", "Aromatic", "% Aromatic"]

for familiesDataFrameLines in range(0,familiesDataframe.shape[0]):
    familyName = str(familiesDataframe.iloc[familiesDataFrameLines,0])
    if familyName not in columnNamesList:
        columnNamesList.append(familyName)


#Create statistics dataframes with columns named defined above
proteinStatisticsDataframe=pandas.DataFrame(columns=columnNamesList)
dprStatisticsDataframe=pandas.DataFrame(columns=columnNamesList)
nodprStatisticsDataframe=pandas.DataFrame(columns=columnNamesList)
negativeStatisticsDataframe=pandas.DataFrame(columns=columnNamesList)


"""
Write statitian Function that takes input_dataframe,output_dataframe,
output_file_name and control_variable and analises sequences in input_dataframe
saves them in output_dataframe and outputs a excell sheet with output_dataframe 
information that being the size of the sequence, residue count and percentage
in said sequence and family and percentage of residues in that family for the sequence
"""

def statitian(input_dataframe,output_dataframe,output_file_name,control_variable):
    
    """
    Iterate trought lines in input data frame and take uniprot Id and 
    sequence from it , acording to control_variable
    """
    for dataframeLineIndex in range(0, input_dataframe.shape[0]):
        
        uniprotId = str(input_dataframe.iloc[dataframeLineIndex,0])
        
        output_dataframe.at[dataframeLineIndex, "Uniprot ID"] = uniprotId
        
        if control_variable == "PROTEIN" :
            sequence = str(input_dataframe.iloc[dataframeLineIndex,1])
            output_dataframe.at[dataframeLineIndex, "Sequence"] = sequence
        if control_variable == "SEQUENCE" :
            sequence = str(input_dataframe.iloc[dataframeLineIndex,2])
            output_dataframe.at[dataframeLineIndex, "Sequence"] = sequence
        
        """
        Write to Dataframe size of the sequence, Amino acid count and percentage of said Amino acid
        , family and percentage of Amino acids in that family 
        in said sequence
        """
        output_dataframe.at[dataframeLineIndex, "Total aa"] = len(sequence)     
        output_dataframe.at[dataframeLineIndex, "A"] = sequence.count("A")
        output_dataframe.at[dataframeLineIndex, "%A"] = (sequence.count("A") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "G"] = sequence.count("G")
        output_dataframe.at[dataframeLineIndex, "%G"] = (sequence.count("G") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "V"] = sequence.count("V")
        output_dataframe.at[dataframeLineIndex, "%V"] = (sequence.count("V") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "L"] = sequence.count("L")
        output_dataframe.at[dataframeLineIndex, "%L"] = (sequence.count("L") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "I"] = sequence.count("I")
        output_dataframe.at[dataframeLineIndex, "%I"] = (sequence.count("I") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "P"] = sequence.count("P")
        output_dataframe.at[dataframeLineIndex, "%P"] = (sequence.count("P") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "F"] = sequence.count("F")
        output_dataframe.at[dataframeLineIndex, "%F"] = (sequence.count("F") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "M"] = sequence.count("M")
        output_dataframe.at[dataframeLineIndex, "%M"] = (sequence.count("M") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "W"] = sequence.count("W")
        output_dataframe.at[dataframeLineIndex, "%W"] = (sequence.count("W") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "C"] = sequence.count("C")
        output_dataframe.at[dataframeLineIndex, "%C"] = (sequence.count("C") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "D"] = sequence.count("D")
        output_dataframe.at[dataframeLineIndex, "%D"] = (sequence.count("D") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "E"] = sequence.count("E")
        output_dataframe.at[dataframeLineIndex, "%E"] = (sequence.count("E") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "R"] = sequence.count("R")
        output_dataframe.at[dataframeLineIndex, "%R"] = (sequence.count("R") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "K"] = sequence.count("K")
        output_dataframe.at[dataframeLineIndex, "%K"] = (sequence.count("K") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "H"] = sequence.count("H")
        output_dataframe.at[dataframeLineIndex, "%H"] = (sequence.count("H") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "N"] = sequence.count("N")
        output_dataframe.at[dataframeLineIndex, "%N"] = (sequence.count("N") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "Q"] = sequence.count("Q")
        output_dataframe.at[dataframeLineIndex, "%Q"] = (sequence.count("Q") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "S"] = sequence.count("S")
        output_dataframe.at[dataframeLineIndex, "%S"] = (sequence.count("S") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "T"] = sequence.count("T")
        output_dataframe.at[dataframeLineIndex, "%T"] = (sequence.count("T") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "Y"] = sequence.count("Y")
        output_dataframe.at[dataframeLineIndex, "%Y"] = (sequence.count("Y") * 100) / len(sequence)
        output_dataframe.at[dataframeLineIndex, "Hydrofobic"] = (sequence.count("L")+ sequence.count("W")+ sequence.count("A")+ sequence.count("I")+ sequence.count("M")+ sequence.count("V")+ sequence.count("F")+ sequence.count("C"))
        output_dataframe.at[dataframeLineIndex, "% Hydrofobic"] = (((sequence.count("L") + sequence.count("W") + sequence.count("A") + sequence.count("I") + sequence.count("M") + sequence.count("V") + sequence.count("F") + sequence.count("C")) * 100)/ len(sequence))
        output_dataframe.at[dataframeLineIndex, "Negative Charged"] = sequence.count("D") + sequence.count("E")
        output_dataframe.at[dataframeLineIndex, "Negative Charged"] = (((sequence.count("D") + sequence.count("E")) * 100) / len(sequence))
        output_dataframe.at[dataframeLineIndex, "Positive Charged"] = (sequence.count("K") + sequence.count("R") + sequence.count("H"))
        output_dataframe.at[dataframeLineIndex, "Positive Charged"] = (((sequence.count("K") + sequence.count("R") + sequence.count("H")) * 100)/ len(sequence))
        output_dataframe.at[dataframeLineIndex, "Polar"] = (sequence.count("G") + sequence.count("S") + sequence.count("T") + sequence.count("N") + sequence.count("Q") + sequence.count("Y"))
        output_dataframe.at[dataframeLineIndex, "% Polar"] = (((sequence.count("G") + sequence.count("S") + sequence.count("T") + sequence.count("N") + sequence.count("Q") + sequence.count("Y")) * 100)/ len(sequence))
        output_dataframe.at[dataframeLineIndex, "Aromatic"] = sequence.count("W") + sequence.count("Y") + sequence.count("F")
        output_dataframe.at[dataframeLineIndex, "% Aromatic"] = (((sequence.count("W") + sequence.count("Y") + sequence.count("F")) * 100)/ len(sequence))
        
        for familiesIndex in range(0,familiesDataframe.shape[0]):
            if str(familiesDataframe.iloc[familiesIndex,1]) ==  uniprotId:
                for columnIndex in columnNamesList:
                    if str(familiesDataframe.iloc[familiesIndex,0]) == str(columnIndex):
                        output_dataframe.at[dataframeLineIndex,str(columnIndex)] = "X" 
        
    if control_variable == "PROTEIN" :
        
        output_dataframe.drop_duplicates(inplace=True)
    
    
    
    #Write output_dataframe to excel file in outputPathFileName Path
    outputPathFileName = (Path(outputFilePath)/output_file_name)
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:
        output_dataframe.to_excel(writer, sheet_name='STATS', index=False)
        
# Use of statitian function for data processing
statitian(dprDataframe,proteinStatisticsDataframe,"Protein Stats.xlsx","PROTEIN")
statitian(dprDataframe,dprStatisticsDataframe,"DPR Stats.xlsx","SEQUENCE")
statitian(nodprDataframe,nodprStatisticsDataframe,"NO DPR Stats.xlsx","SEQUENCE")
statitian(negativeDataframe,negativeStatisticsDataframe,"Negative Stats.xlsx","PROTEIN")

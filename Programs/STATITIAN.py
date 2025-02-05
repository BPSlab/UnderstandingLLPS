#Import Libraries
from pathlib import Path
import pandas
import os

#Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

#Define output Paths 
outputFilePath = (Path(mainProjectDirectory)/'Results'/'Statistics')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)

#Read Databases from Excel sheets 
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


"""
# Write statitian Function that takes input_dataframe,output_file_name and control_variable and analyses sequences in input_dataframe,
saves them in output_dataframe and outputs an excel sheet with output_dataframe information including the Uniprot ID, the sequence, 
and residue count and percentage 
"""

def statitian(input_dataframe,output_file_name,control_variable):
    
    columnNamesList = ["Uniprot ID", "Sequence", "Total aa", "A", "%A", "G", "%G", "V", "%V", "L", "%L", "I", "%I", "P", "%P", "F", "%F", 
                    "M", "%M", "W", "%W", "C", "%C", "D", "%D", "E", "%E", "R", "%R", "K", "%K", "H", "%H", "N", "%N", "Q", "%Q", "S", "%S", 
                    "T", "%T", "Y", "%Y", "Hydrofobic", "% Hydrofobic", "Negative Charged", "% Negative Charged", "Positive Charged", 
                    "% Positive Charged", "Polar", "% Polar", "Aromatic", "% Aromatic"]


    #Create statistics dataframes with columns named defined above
    outputDataframe=pandas.DataFrame(columns=columnNamesList)
    
    """
    Iterate trought lines in input data frame and take Uniprot ID and 
    sequence from it , acording to control_variable
    """
    for dataframeLineIndex in range(0, input_dataframe.shape[0]):
        
        uniprotId = str(input_dataframe.loc[dataframeLineIndex,'UniprotID'+' '+control_variable])
        
        outputDataframe.at[dataframeLineIndex, "Uniprot ID"] = uniprotId
        
        sequence = str(input_dataframe.loc[dataframeLineIndex,control_variable])
        
        outputDataframe.at[dataframeLineIndex, "Sequence"] = sequence

        
        """
        Write to Dataframe size of the sequence, residue count and percentage, 
        """
        outputDataframe.at[dataframeLineIndex, "Total aa"] = len(sequence)     
        outputDataframe.at[dataframeLineIndex, "A"] = sequence.count("A")
        outputDataframe.at[dataframeLineIndex, "%A"] = (sequence.count("A") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "G"] = sequence.count("G")
        outputDataframe.at[dataframeLineIndex, "%G"] = (sequence.count("G") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "V"] = sequence.count("V")
        outputDataframe.at[dataframeLineIndex, "%V"] = (sequence.count("V") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "L"] = sequence.count("L")
        outputDataframe.at[dataframeLineIndex, "%L"] = (sequence.count("L") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "I"] = sequence.count("I")
        outputDataframe.at[dataframeLineIndex, "%I"] = (sequence.count("I") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "P"] = sequence.count("P")
        outputDataframe.at[dataframeLineIndex, "%P"] = (sequence.count("P") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "F"] = sequence.count("F")
        outputDataframe.at[dataframeLineIndex, "%F"] = (sequence.count("F") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "M"] = sequence.count("M")
        outputDataframe.at[dataframeLineIndex, "%M"] = (sequence.count("M") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "W"] = sequence.count("W")
        outputDataframe.at[dataframeLineIndex, "%W"] = (sequence.count("W") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "C"] = sequence.count("C")
        outputDataframe.at[dataframeLineIndex, "%C"] = (sequence.count("C") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "D"] = sequence.count("D")
        outputDataframe.at[dataframeLineIndex, "%D"] = (sequence.count("D") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "E"] = sequence.count("E")
        outputDataframe.at[dataframeLineIndex, "%E"] = (sequence.count("E") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "R"] = sequence.count("R")
        outputDataframe.at[dataframeLineIndex, "%R"] = (sequence.count("R") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "K"] = sequence.count("K")
        outputDataframe.at[dataframeLineIndex, "%K"] = (sequence.count("K") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "H"] = sequence.count("H")
        outputDataframe.at[dataframeLineIndex, "%H"] = (sequence.count("H") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "N"] = sequence.count("N")
        outputDataframe.at[dataframeLineIndex, "%N"] = (sequence.count("N") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "Q"] = sequence.count("Q")
        outputDataframe.at[dataframeLineIndex, "%Q"] = (sequence.count("Q") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "S"] = sequence.count("S")
        outputDataframe.at[dataframeLineIndex, "%S"] = (sequence.count("S") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "T"] = sequence.count("T")
        outputDataframe.at[dataframeLineIndex, "%T"] = (sequence.count("T") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "Y"] = sequence.count("Y")
        outputDataframe.at[dataframeLineIndex, "%Y"] = (sequence.count("Y") * 100) / len(sequence)
        outputDataframe.at[dataframeLineIndex, "Hydrofobic"] = (sequence.count("L")+ sequence.count("W")+ sequence.count("A")+ sequence.count("I")+ sequence.count("M")+ sequence.count("V")+ sequence.count("F")+ sequence.count("C"))
        outputDataframe.at[dataframeLineIndex, "% Hydrofobic"] = (((sequence.count("L") + sequence.count("W") + sequence.count("A") + sequence.count("I") + sequence.count("M") + sequence.count("V") + sequence.count("F") + sequence.count("C")) * 100)/ len(sequence))
        outputDataframe.at[dataframeLineIndex, "Negative Charged"] = sequence.count("D") + sequence.count("E")
        outputDataframe.at[dataframeLineIndex, "% Negative Charged"] = (((sequence.count("D") + sequence.count("E")) * 100) / len(sequence))
        outputDataframe.at[dataframeLineIndex, "Positive Charged"] = (sequence.count("K") + sequence.count("R") + sequence.count("H"))
        outputDataframe.at[dataframeLineIndex, "% Positive Charged"] = (((sequence.count("K") + sequence.count("R") + sequence.count("H")) * 100)/ len(sequence))
        outputDataframe.at[dataframeLineIndex, "Polar"] = (sequence.count("G") + sequence.count("S") + sequence.count("T") + sequence.count("N") + sequence.count("Q") + sequence.count("Y"))
        outputDataframe.at[dataframeLineIndex, "% Polar"] = (((sequence.count("G") + sequence.count("S") + sequence.count("T") + sequence.count("N") + sequence.count("Q") + sequence.count("Y")) * 100)/ len(sequence))
        outputDataframe.at[dataframeLineIndex, "Aromatic"] = sequence.count("W") + sequence.count("Y") + sequence.count("F")
        outputDataframe.at[dataframeLineIndex, "% Aromatic"] = (((sequence.count("W") + sequence.count("Y") + sequence.count("F")) * 100)/ len(sequence))
        
    outputDataframe = outputDataframe.replace('nan',None)
    outputDataframe.dropna(inplace=True)
    #Write output_dataframe to excel file in outputPathFileName Path
    outputPathFileName = (Path(outputFilePath)/output_file_name)
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:
        outputDataframe.to_excel(writer, sheet_name='STATS', index=False)
        
# Use of statitian function for data processing
statitian(proteinDataframe,"Protein Stats.xlsx","Protein")
statitian(dprDataframe,"Protein DPR Stats.xlsx","DPR")
statitian(nodprDataframe,"Protein NODPR Stats.xlsx","NODPR")
statitian(negativeDataframe,"Negative Stats.xlsx","Protein")

statitian(rnabindingDataframe,"RNA Binding Protein Stats.xlsx","Protein")
statitian(rnabindingDataframe,"RNA Binding DPR Stats.xlsx","DPR")
statitian(rnabindingDataframe,"RNA Binding NODPR Stats.xlsx","NODPR") 

statitian(dnabindingDataframe,"DNA Binding Protein Stats.xlsx","Protein")
statitian(dnabindingDataframe,"DNA Binding DPR Stats.xlsx","DPR")
statitian(dnabindingDataframe,"DNA Binding NODPR Stats.xlsx","NODPR") 

statitian(chromatinbindingDataframe,"Chromatin Binding Protein Stats.xlsx","Protein")
statitian(chromatinbindingDataframe,"Chromatin Binding DPR Stats.xlsx","DPR")
statitian(chromatinbindingDataframe,"Chromatin Binding NODPR Stats.xlsx","NODPR") 

statitian(regulationDataframe,"Regulation Protein Stats.xlsx","Protein")
statitian(regulationDataframe,"Regulation DPR Stats.xlsx","DPR")
statitian(regulationDataframe,"Regulation NODPR Stats.xlsx","NODPR") 

statitian(hydrolaseDataframe,"Hydrolase Protein Stats.xlsx","Protein")
statitian(hydrolaseDataframe,"Hydrolase DPR Stats.xlsx","DPR")
statitian(hydrolaseDataframe,"Hydrolase NODPR Stats.xlsx","NODPR") 

statitian(structureDataframe,"Structure Protein Stats.xlsx","Protein")
statitian(structureDataframe,"Structure DPR Stats.xlsx","DPR")
statitian(structureDataframe,"Structure NODPR Stats.xlsx","NODPR") 

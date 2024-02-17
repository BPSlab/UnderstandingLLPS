#Import Libraries
from pathlib import Path
import pandas
import os


#Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

#Data Paths 
outputFilePath = (Path(mainProjectDirectory)/'Results'/'Motifs')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)
    
#Read Statistic Data frames
dprStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Statistics'/'DPR Stats.xlsx'))
rnaBindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'RNA binding.xlsx'))
dnaStatisticsDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'DNA binding.xlsx'))
regulationDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Regulation.xlsx'))
chromatinbindingDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Chromatin binding.xlsx'))
hydrolaseDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Hydrolase.xlsx'))
structureDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Structure.xlsx'))


#Define column names and make Data frame with columns using said column names
columnNames = ["Format","Motif", "Motif Presence", "Motif Frequency"]
motifDataframe = pandas.DataFrame(columns = columnNames)

"""
define the generateCombinations function manages the generation of combinations 
based on a given input string, z_value, and value_list.
It initializes an empty list to store the combinations and calls the
generateCombinationsHelper function to generate said combinations.
Once the combinations are generated, it returns the list for further processing
described bellow.
"""
def generateCombinations(string, z_value, value_list):
    combinations = []
    generateCombinationsHelper(string, z_value, value_list, 0, [], combinations)
    return combinations

"""
Define the generateCombinationsHelper function to constructing combinations.
With parameters specifying the input string, substitution values, 
current index, and ongoing combination, it navigates through the input string,
appending values as appropriate.
When encountering placeholders 'X' , it explores all possible 
substitutions constructing and appending combinations until 
the entire string is traversed. Once a combination is complete, it's added to 
the list of combinations. 
"""
def generateCombinationsHelper(string, z_value, value_list, index, current_combination, combinations):
    if index == len(string):
        combinations.append(''.join(current_combination))
        return

    if string[index] == 'Z':
        current_combination.append(z_value)
        generateCombinationsHelper(string, z_value, value_list, index + 1, current_combination, combinations)
        current_combination.pop()
    elif string[index] == 'X':
        for value in value_list:
            current_combination.append(value)
            generateCombinationsHelper(string, z_value, value_list, index + 1, current_combination, combinations)
            current_combination.pop()

def main(input_dataframe):
    #Take user input first asking the format for sequences to be generated then asking the Fixed Value 'Z'
    inputString = input("Enter the string with 'X' and 'Z' characters: ")
    zValue = input("Enter the value of 'Z': ")
    outputPathName = input("Enter the name of the OutputFolder: ")
    #Create listOfAminoacids containing all possible values for the 'X' Variable
    listOfAminoacids = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    #Run generateCombinations function to Make the List of combinations relative to the Motifs 
    combinations = generateCombinations(inputString, zValue, listOfAminoacids)
    
    verticalIndex = 0    
    #Iterate trougth Motifs and inicialize count variables for each motif
    for combo in combinations:
        masterCount = 0 
        localCount = 0
        totalCount = 0
        motifToCount = str(combo)
    
        """
        iterate trougth lines in dprStatisticsDataframe and take Sequence to
        count current iteration Motif presence(number of sequences the motif appears in)
        , and frequency (number of times the motif appears in a given sequence)
        """
        for inputDataframeLines in range (0,input_dataframe.shape[0]):
            sequence = str(input_dataframe.iloc[inputDataframeLines,1])
            masterCount = sequence.count(motifToCount)        
            if masterCount != 0 :
                localCount = localCount + 1
                totalCount = totalCount + masterCount                     
    
        #Write format motif String and Motif Presence and frequency      
        motifDataframe.at[verticalIndex,"Format"] = inputString
        motifDataframe.at[verticalIndex, "Motif"] = combo  
        motifDataframe.at[verticalIndex, "Motif Presence"] = localCount
        motifDataframe.at[verticalIndex, "Motif Frequency"] = totalCount
        verticalIndex += 1
    
    #Make file named based on sequence Format
    
    
    outputFilePath = (Path(mainProjectDirectory)/'Results'/'Motifs'/outputPathName)
    outputFileName = inputString.replace('Z',zValue) + '.xlsx'
    outputPathFileName = (Path(outputFilePath)/outputFileName)    
    if not os.path.exists(outputFilePath):
        os.makedirs(outputFilePath)
    #Write output_dataframe to excel file in outputPathFileName Path
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:
        motifDataframe.to_excel(writer, sheet_name='Motif Search', index=False)      

main(dprStatisticsDataframe)
    
     
    

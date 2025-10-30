from pathlib import Path
import pandas
import os

# Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

outputFilePath = (Path(mainProjectDirectory)/'Results'/'Peptides Cider')
if not os.path.exists(outputFilePath):
    os.makedirs(outputFilePath)
    
peptideDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Results'/'Peptides Cider'/'Peptides Cider.xlsx'))

def dataFrameFormater(input_dataframe,output_filename):

    input_dataframe.drop(input_dataframe[input_dataframe['FCR'] < 0.13].index, inplace = True)
    input_dataframe.drop(input_dataframe[input_dataframe['FCR'] > 0.28].index, inplace = True)
    
    print(input_dataframe)
    
    input_dataframe.drop(input_dataframe[input_dataframe['NCPR'] < 0].index, inplace = True)
    input_dataframe.drop(input_dataframe[input_dataframe['NCPR'] > 0.084].index, inplace = True)
    
    print(input_dataframe)
    
    input_dataframe.drop(input_dataframe[input_dataframe['Mean Hydropathy'] < 3.15].index, inplace = True)
    input_dataframe.drop(input_dataframe[input_dataframe['Mean Hydropathy'] > 3.76].index, inplace = True)
     
    print(input_dataframe)  
         
    input_dataframe.drop(input_dataframe[input_dataframe['Fraction Of Disorder Pormoting Residues'] < 0.79].index, inplace = True)
    input_dataframe.drop(input_dataframe[input_dataframe['Fraction Of Disorder Pormoting Residues'] > 0.86].index, inplace = True)
    
    print(input_dataframe)
    
    input_dataframe.drop(input_dataframe[input_dataframe['Kappa'] < 0.163].index, inplace = True)
    input_dataframe.drop(input_dataframe[input_dataframe['Kappa'] > 0.242].index, inplace = True)

    print(input_dataframe)
 

    input_dataframe.reset_index(drop = True,inplace = True )
    
    outputPathFileName = (Path(outputFilePath)/output_filename)
    
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:
    
        input_dataframe.to_excel(writer, sheet_name='DPR', index=True)

    return(input_dataframe)

df = dataFrameFormater(peptideDataframe,'Cider for peptides filtered.xlsx')
    
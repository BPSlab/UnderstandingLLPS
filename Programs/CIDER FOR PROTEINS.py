#Import Libraries
from pathlib import Path
import pandas 
import os
from localcider.sequenceParameters import SequenceParameters

#Define Main Project Diectory
global mainProjectDirectory
mainProjectDirectory = Path.cwd().parent

#Read from excel
dprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'),'DPR')
nodprDataframe = pandas.read_excel((Path(mainProjectDirectory)/'Datasets'/'Database.xlsx'),'NODPR')

"""
Turn a Sequence String into a Sequence Object with Sequence Parameters,
then run the get function for the following localCider parameters on the Sequence Object:
       
    "NCPR"
    "Isoelectric Point"
    "FCR"
    "Omega"
    "PPII Propensity"
    "SCD"
    "Hydropathy"
    "Fraction Expanding"
    "Fraction Negative"
    "Fraction Positive"
    "Kappa"
    "Length"
    "linear Sigma"
    "Mean Hydropathy"
    "Mean Net Charge"
    "Fraction Of Disorder Pormoting Residues
    
Returns List of Parameters.  
"""

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

   #print(SeqObject.get_phasePlotRegion())
   #print(SeqObject.show_phaseDiagramPlot())


"""
Reads from input_dataframe the Type of sequence defined in control_variable.
For each sequence call the ciderAnalisis UDF and write results to output Dataframe, 
then writes Data frame to output_file_name Named .excel File.

"""   
def cider(input_dataframe,output_file_name,control_variable):
    
    #Define column names and make Data frame with columns using said column names
    columnNames = ["Sequence",
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
    
    sequenceDataframe = pandas.DataFrame(columns = columnNames)

    for inputDataframeLines in range (0,input_dataframe.shape[0]):
        
        sequence = str(input_dataframe.loc[inputDataframeLines,control_variable])
        
        sequenceCiderResults = ciderAnalisis(sequence)
           
        sequenceDataframe.at[inputDataframeLines,"Sequence"] = sequence
        sequenceDataframe.at[inputDataframeLines,"NCPR"] = sequenceCiderResults[0]
        sequenceDataframe.at[inputDataframeLines,"Isoelectric Point"] = sequenceCiderResults[1]
        sequenceDataframe.at[inputDataframeLines,"FCR"] = sequenceCiderResults[2]
        sequenceDataframe.at[inputDataframeLines,"Omega"] = sequenceCiderResults[3]
        sequenceDataframe.at[inputDataframeLines,"PPII Propensity"] = sequenceCiderResults[4]
        sequenceDataframe.at[inputDataframeLines,"SCD"] = sequenceCiderResults[5]
        sequenceDataframe.at[inputDataframeLines,"Hydropathy"] = sequenceCiderResults[6]
        sequenceDataframe.at[inputDataframeLines,"Fraction Expanding"] = sequenceCiderResults[7]
        sequenceDataframe.at[inputDataframeLines,"Fraction Negative"] = sequenceCiderResults[8]
        sequenceDataframe.at[inputDataframeLines,"Fraction Positive"] = sequenceCiderResults[9]
        sequenceDataframe.at[inputDataframeLines,"Kappa"] = sequenceCiderResults[10]
        sequenceDataframe.at[inputDataframeLines,"Length"] = sequenceCiderResults[11]
        sequenceDataframe.at[inputDataframeLines,"linear Sigma"] = sequenceCiderResults[12]
        sequenceDataframe.at[inputDataframeLines,"Mean Hydropathy"] = sequenceCiderResults[13]
        sequenceDataframe.at[inputDataframeLines,"Mean Net Charge"] = sequenceCiderResults[14]
        sequenceDataframe.at[inputDataframeLines,"Fraction Of Disorder Pormoting Residues"] = sequenceCiderResults[15]


    outputFilePath = (Path(mainProjectDirectory)/'Results'/'Cider for proteins')
    if not os.path.exists(outputFilePath):
        os.makedirs(outputFilePath)
    
    outputPathFileName = (Path(outputFilePath)/output_file_name) 
        
    #Write output_dataframe to excel file in outputPathFileName Path
    with pandas.ExcelWriter(outputPathFileName, engine='openpyxl') as writer:
        sequenceDataframe.to_excel(writer, sheet_name='Cider Analysis', index=False)      

#Run cider for the given inputs
cider(dprDataframe,'Cider DPR.xlsx','DPR')
cider(nodprDataframe,'Cider NODPR.xlsx','NODPR')




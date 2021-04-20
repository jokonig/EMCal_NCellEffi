Main code is MakeNCellEffiPureGammas.cxx which can be started with the bash script startTagging.sh

The code will process data and MC input (Tree based) from the train output from PWGGA/GammaConv/AliAnalysisTaskClusterQA.cxx

The output from MakeNCellEffiPureGammas.cxx is written into a root file (in rootFiles/NameOfRootFile.root).

The NCell efficiency can be calculated and plotted using the macro PlottingNCellEffi_Methods.cxx.
In this code different methods for the tagging procedure can be used (wide, narrow etc.). The Standard method is "Wide".

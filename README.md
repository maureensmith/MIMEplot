# MIMEplot
R routines for plotting the output tables of MIMEAnTo/MIMEAn2

Running MIMEplot
-------------------

Run the R script 
```
Rscript plotMIMEAnToResults.R <MIMEAnTo_directory> <referenceFile> [<suffix>]
```

withe following parameters

| parameter       | type          | description  |
| :---  |:---:| :----------------|
| MIMEAnTo_directory         | (string)      |   directory where the subdirectories of the MIMEAnTo/MIMEAn2 results are lying |
| referenceFile         | (string)      |   reference sequence in fasta format |
| suffix          | (string)      |   (optional) suffix of the table names, for example when using different options |

The MIMEAnToDirectory contains one or more subdirectories, where the results tables of [MIMEAnTo](https://github.com/maureensmith/mimeanto) or [MIMEAn2](https://github.com/maureensmith/mimean2) are located. 

The script generates individual plots for each of the subfolders: The smoothed maxKd and the "long plot", meaning a boxplot showing the Kds for each mutation at each position of the regarded sequence.
Another plot is generated in the MIMEAnTo_directory, overlaying all smoothed maxKds for each of the sub analyses (colour-coded).

If a suffix was given for the generation of the result tables, it has to be submitted in the plot call, too.
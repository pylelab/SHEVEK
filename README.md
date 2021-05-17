# SHEVEK #
A computational method for predicting intramolecular and intermolecular biopolymer interactions

## Install ##
The original SHEVEK was originally developed using the (now obsolete) Microsoft Visual C++ 6.0. It is not partly rewritten to make it compatible with gcc on Linux.
```bash
make
```

## Run the SHEVEK program interactively ##

```bash
./shevek
```
```
Enter alignment file name (include extension): example.pir
```

A.Input Alignment File Format

The input file format can be in FASTA/PIR/NBRF formats. The ONLY two reserved characters (i.e. characters noted especially by the program) are ‘>’ and ‘-‘. ‘>’ denotes a new sequence, which is followed by the sequence name. ‘-‘ is the only gap character recognized by SHEVEK.

NOTE: the sequence name (denoted by the character ‘>’) MUST be followed by a “carriage return” –i.e. it must be on its own separate line. There does NOT, however, need to be an extra line between the sequence name and the start of actual sequence characters. 

NOTE: all other non-reserved characters are acceptable as representatives of biopolymer units. 

Alignments can be generated using any number of software programs. Suggest CLUSTAL and BIOEDIT.

```
Enter the number of the REFERENCE sequence: 2
```

B.Reference Sequence

The reference sequence is the sequence of interest. i.e. the mouse biologist would type 2, since the second sequence in the alignment EXAMPLE.TXT is the mouse protein. This will generate the file POSITION.TXT which will contain a table that relates the alignment position of sequence 2 to its residue number and type. (i.e. residue 1 in sequence 2 may not be alignment position 1, due to the presence of gap characters.)

```
Set INTEGER PRECISION of P-values to be calculated to [default:8]: 8
```

C.Integer Precision

THIS IS AN ADVANCED FEATURE. See paper. The default value of 8 should be used. This number MUST be an integer. It specifies the precision to which p scores will be calculated. 

 The following questions allow for SHEVEK to split an
 alignment among different processors. Use default values
 to analyze the entire alignment.

```
Enter the column position START value [default:0]: 0
Enter the column position STOP value [Default:-1]: -1
```

D.Alignment Analysis

THIS IS AN ADVANCED FEATURE. See paper. Default values of 0 and –1 should be used. ‘-1’ specifies the end of the alignment. For flexibility, SHEVEK can be set to analyze only a section of an alignment. Note that if the start is 0 and the stop is 4, the program will compare position 0 against 1,2,3;4; position 1 against 2,3,4; position 2 against 3,4; position 3 against 4.

```
Scoring position #125's possible interactions. . . please wait.
   Scoring interaction w/position #126. Computing.
   Scoring interaction w/position #127. Computing.
```

You will see lines like these  until the entire alignment has been processed -- at that point, you will see:

```
********************************************************************
********************************************************************
            Scoring and Standardization Complete
********************************************************************
********************************************************************

You can view your "alldata" file using any plotting program/EXCEL.
Shevek will now calculate suggested lower thresholds.
********************************************************************
```

E. Output Data Files

At this point, the SHEVEK program will generate two files, beginning with ‘alldata’ and ‘nulldata’. These files will be followed by two numbers, which designate the start and stop of the alignment. In this case, since EXAMPLE.TXT contains sequences whose alignment length is 153, the output file names are ‘alldata0x153.txt’ and ‘nulldata0x153.txt’.

The Alldata FILE CONTAIN THE RELEVANT OUTPUT CALCULATIONS. The nulldata file contains pair-wise combinations of positions whose scores are irrelevant because: a) they contain only gaps or many gaps; b) they contain only one character type; c) they result in tables that are sensitive; d) they result in imbalanced tables; d) their probability is 1.

```
Enter alldata file name for distribution analysis: alldata0x153.txt
```

F. Making Predictions

Much of the screen output is an ADVANCED FEATURE, and can be ignored, except in the case of WARNINGS. The most relevant text is the suggested thresholds:

```
*************************************************
SUGGESTED Lower -Log(P) threshold: 3.7
SUGGESTED Lower V threshold: 0.400
SUGGESTED protein DF threshold: 9.0
*************************************************
```

The program will then ask you to re-type the above filename.:

Enter file name for progressive threshing (alldataXX.txt):  alldata0x153.txt

THIS IS AN ADVANCED FEATURE. WHY does the computer ask this repeat question? Isn’t that just stupid? Yes and No. Yes, it does not necessarily need to, but doing so allows for flexibility both in terms of programming and data anlysis. 

```
***********************************************
PARAMETER ENTRY (use numbers ABOVE as defaults)
***********************************************

Set -log(P) threshold to less than or equal to [range:0-inf]: 3.7

Set V threshold to less than or equal to [range:0-1]: .4

Set df threshold to less than [default for RNA/DNA:100]: 9

Set Sensitivity threshold to less than [default:1.0]: 1
```

THIS IS AN ADVANCED FEATURE. ANSWER these questions using the SUGGESTED parameters found above. In the case of sensitivity, use the default value of 1. Again, why bother to ask the user to type numbers the computer has just output? Advanced flexibility.

NOTE: if analyzing only  RNA or DNA sequences, the dfthreshold should be set to the default of 100.

Next follows a large amount of scrolling screen text. THIS IS AN ADVANCED FEATURE. This text can be ignored. It describes the “chains of associations” that occur at each theshold, starting with the entered p threshold and ending with the maximum p that occurs in the dataset.  Specifically, the third column of negative numbers specifies pair-wise combinations of positions that are linked (i.e. part of the same chain.) The fourth column is filled with –1 for all pair-wise combinations of positions that “lost” – i.e. are not the best for the chain to which they belong. “tcounter” denotes the robustness of a prediction (see PROGRESSIVE THRESHING in paper.)

```
********************************************************************
            PREDICTIONS COMPLETE: SEE 'predictions.txt'
********************************************************************
NOTE: Predictions with '666' in the last column denote non-robust
predictions; these should be disregarded. Predictions with negative
values in column 3 denote predictions whose associations could not
be untangled; these should also be disregarded.

NOTE: predictions identify alignment positions. In order to identify
how these alignment positions relate to the reference sequence, see
the file 'position.txt.' Alternatively, use distcalc.exe program.
********************************************************************
```

SHEVEK has created a file entitled PREDICTIONS.TXT that contains the predictions of the SHEVEK program for the alignment EXAMPLE.TXT.

```
Apply Misalignment Algorithm to the output predictions? (Y/N) y
```

G. Misalignment Algorithm.

The misalignment algorithm is a simple method of identifying sequences with unusual character combinations. It is an additional, accessory  feature of SHEVEK.

```
********************************************************************
              Misalignment Identifcation Process
********************************************************************
SHEVEK will now attempt to identify misaligned sequences. Shevek
will output a file entitled, 'misalign.txt', which will list
sequences identified as possibly misaligned by categorical
statistical analysis.

This information can be used to re-align the sequence alignment.
After such re-alignment, rerun SHEVEK with the new alignment for
increased prediction accuracy.

Stop such iteration when the number/identity of possibly misaligned
sequences no longer changes -- or further re-alignment is not deemed
possible.
********************************************************************

Enter predictions file name (predictions.txt): predictions.txt
```

Follows a large amount of text. SHEVEK creates files entitled “MISALIGN.TXT”, which contains detailed instructions on possible regions of misalignment. SHEVEK also creates a series of files that being with “ARtables", with the numbers designating the two alignment column positions to which the table corresponds. The numbers in these tables specify which specific character combinations are favorable or unfavorable. See paper.

```
PROGRAM EXECUTION COMPLETED

Type 'exit' to exit program: exit
```

## Reference ## 
Phillip S. Pang, Eckhard Jankowsky, Leven M. Wadley, and Anna Marie Pyle.
["Prediction of functional tertiary interactions and intermolecular interfaces from primary sequence data."](https://doi.org/10.1002/jez.b.21024)
Journal of Experimental Zoology Part B: Molecular and Developmental Evolution 304, no. 1 (2005): 50-63.

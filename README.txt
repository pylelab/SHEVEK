**************************************************************
README README README README README README README
**************************************************************
I. Introduction:

	THIS SOURCE CODE IS NOT TO BE DUPLICATED WITHOUT PERMISSION

	This Disk contains the source code (written in the C/C++ programming language)
	for US Patent Application of Pang et al. for  "COMPUTATIONAL METHOD FOR PREDICTING
 			INTRAMOLECULAR AND INTERMOLECULAR BIOPOLYMER INTERACTIONS"
			https://patents.google.com/patent/US20050026145

	By Phillip S. Pang, Eckhard Jankowsky and Anna Marie Pyle
	For more information contact: phillip.pang@stanfordalumni.org

	Filing By: Baker Botts LLP, 30 Rockefeller Plaza, New York, NY, 10112.
 

II. Files:

	The Microsoft Visual C++ 6.0 Programming environment was used:

	The workspace file is: Shevek.dsw
	The project file is: Shevek.dsp

	Source files include:
		
		ApplyThresh.c
		ChainElim.c
		definitions.h
		FASTEXP2.C
		FASTEXP2.H
		FASTEXP2.INC
		GETTIMER.C
		MainExactPScore.c
		MainMisalign.c
		MainScreening.c
		MainShevek.c
		MisalignNRPang.c
		nrutilp.h
		NumRecPang.c
		NumRecUtilities.c
		PRCERR.C
		RandNumGen.c
		RCONT2p.H
		Rcount2Pang.C
		shevek.h
		SUPPORT.C
		SUPPORT.H	
	
	Text files include:

		degenmsg
		eliminatemsg
		misalignmsg
		predictmsg
		scoremsg
		threshmsg
		titlepage

 
III. General Program Structure:


	FUNCTION			PURPOSE			FOUND IN FILE
	
	main( )				[driving function]	MainShevek.c
	 |
	 |
	 -------- openfile( )     	[reads alignment]	MainShevek.c
	 |
	 -------- score_manager( )	[scoring]		MainShevek.c
	 |
	 -------- screener( )		[screening]		MainScreening.c
	 |
	 -------- misalign_identifier( )[misalignment]		MainMisalign.c

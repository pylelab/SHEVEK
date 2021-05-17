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

	Source files include:
		
		applythresh.c
		chainelim.c
		definitions.h
		fastexp2.c
		fastexp2.h
		fastexp2.inc
		gettimer.c
		mainexactpscore.c
		mainmisalign.c
		mainscreening.c
		mainshevek.c
		misalignnrpang.c
		nrutilp.h
		numrecpang.c
		numrecutilities.c
		prcerr.c
		randnumgen.c
		rcont2p.h
		rcount2pang.c
		shevek.h
		support.c
		support.h	
	
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
	
	main( )				[driving function]	mainshevek.c
	 |
	 |
	 -------- openfile( )     	[reads alignment]	mainshevek.c
	 |
	 -------- score_manager( )	[scoring]		mainshevek.c
	 |
	 -------- screener( )		[screening]		mainscreening.c
	 |
	 -------- misalign_identifier( )[misalignment]		mainmisalign.c

@ECHO OFF

SET ROOTNAME=Kaimal
MOVE  %ROOTNAME%.sum   tstfiles
MOVE  %ROOTNAME%.u     tstfiles
MOVE  %ROOTNAME%.v     tstfiles
MOVE  %ROOTNAME%.w     tstfiles

SET ROOTNAME=vonKarm
MOVE  %ROOTNAME%.dat   tstfiles
MOVE  %ROOTNAME%.hh    tstfiles
MOVE  %ROOTNAME%.sum   tstfiles

SET ROOTNAME=Kaimal_15
MOVE  %ROOTNAME%.sum   tstfiles
MOVE  %ROOTNAME%.u     tstfiles
MOVE  %ROOTNAME%.v     tstfiles
MOVE  %ROOTNAME%.w     tstfiles

SET ROOTNAME=vonKarm_15
MOVE  %ROOTNAME%.dat   tstfiles
MOVE  %ROOTNAME%.hh    tstfiles
MOVE  %ROOTNAME%.sum   tstfiles

SET ROOTNAME=smooth
MOVE  %ROOTNAME%.bin   tstfiles
MOVE  %ROOTNAME%.wnd   tstfiles
MOVE  %ROOTNAME%.sum   tstfiles
MOVE  %ROOTNAME%.cts   tstfiles

SET ROOTNAME=HYDRO_TIDAL
MOVE  %ROOTNAME%.sum   tstfiles
MOVE  %ROOTNAME%.dat   tstfiles
MOVE  %ROOTNAME%.hh    tstfiles
MOVE  %ROOTNAME%.wnd   tstfiles


SET ROOTNAME=KHtest
MOVE  %ROOTNAME%.wnd   tstfiles
MOVE  %ROOTNAME%.sum   tstfiles
MOVE  %ROOTNAME%.cts   tstfiles

SET ROOTNAME=GPLLJ
MOVE  %ROOTNAME%.bts   tstfiles
MOVE  %ROOTNAME%.sum   tstfiles
MOVE  %ROOTNAME%.cts   tstfiles
MOVE  %ROOTNAME%.dat   tstfiles


REM IF NOT "%1"=="LARGE"   GOTO DONE
IF /I NOT "%1"=="LARGE"   GOTO DONE 

set ROOTNAME=GPLLJ_Large
move  %ROOTNAME%.sum   tstfiles
move  %ROOTNAME%.u     tstfiles
move  %ROOTNAME%.v     tstfiles
move  %ROOTNAME%.w     tstfiles
move  %ROOTNAME%.cts   tstfiles
move  %ROOTNAME%.wnd   tstfiles

set ROOTNAME=IECKAI_Large
move  %ROOTNAME%.sum   tstfiles
move  %ROOTNAME%.u     tstfiles
move  %ROOTNAME%.v     tstfiles
move  %ROOTNAME%.w     tstfiles
move  %ROOTNAME%.wnd   tstfiles

:DONE
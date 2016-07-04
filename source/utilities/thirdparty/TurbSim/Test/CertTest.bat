@ECHO OFF
@ECHO.

REM ******************************************************
REM **  Test batch file for TurbSim v1.21, 1-Feb-2007   **
REM ** updated 6/11/2009 to allow MATLAB comparison of  **
REM ** CTS files.  TurbSim v1.41h-bjj                   **
REM ******************************************************

REM ================================================================================================================================
REM  Set up environment variables.  You will probably have to change these.
REM  If you want to use MATLAB to see the results, call this batch program using "CertTest MATLAB"
REM ================================================================================================================================

rem SET PROGRAM=CALL ..\TurbSim_Debug.exe
rem SET PROGRAM=CALL ..\TurbSim64.exe
SET PROGRAM=CALL ..\TurbSim.exe

SET Compare=FC
SET DateTime=DateTime.exe
SET Editor=%SystemRoot%\system32\notepad.exe
SET MATLAB=Call matlab

REM ================================================================================================================================
REM  TurbSim test sequence definition:
REM    Test 1:  Validation test for Kaimal spectrum.
REM    Test 2:  Validation test for von Karman spectrum.
REM    Test 3:  Validation test for Kaimal spectrum with user-specified TI.
REM    Test 4:  Validation test for von Karman spectrum with user-specified TI.
REM    Test 5:  Validation test for Smooth spectrum with coherent events.
REM    Test 6:  Validation test for NWTCUP spectrum with K-H test
REM    Test 7:  Validation test for GP_LLJ spectrum
REM    Test 8:  Validation test for TIDAL spectrum
REM ================================================================================================================================


IF EXIST CertTest.out  DEL CertTest.out
IF EXIST CertTest.m    DEL CertTest.m

REM ================================================================================================================================
ECHO.                                                                             >> CertTest.out
ECHO           ***************************************                            >> CertTest.out
ECHO           **  TurbSim Acceptance Test Results  **                            >> CertTest.out
ECHO           ***************************************                            >> CertTest.out
ECHO.                                                                             >> CertTest.out
ECHO ############################################################################ >> CertTest.out
ECHO # Inspect this file for any differences between your results and the saved # >> CertTest.out
ECHO # results.  Any differing lines and the two lines surrounding them will be # >> CertTest.out
ECHO # listed.  The only differences should be the time stamps at the start of  # >> CertTest.out
ECHO # each file.                                                               # >> CertTest.out
ECHO #                                                                          # >> CertTest.out
ECHO # If you are running on something other than a PC, you may see differences # >> CertTest.out
ECHO # in the last significant digit of many of the numbers.                    # >> CertTest.out
ECHO ############################################################################ >> CertTest.out
ECHO.                                                                             >> CertTest.out
ECHO Date and time this acceptance test was run:                                  >> CertTest.out
%DateTime%                                                                        >> CertTest.out
ECHO.                                                                             >> CertTest.out


ECHO set(0,'DefaultFigureWindowStyle','docked');                                  >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=1
SET ROOTNAME=Kaimal

SET FN1=%ROOTNAME%.sum
SET FN2=%ROOTNAME%.u
SET FN3=%ROOTNAME%.v
SET FN4=%ROOTNAME%.w

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for Kaimal spectrum.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR
IF NOT EXIST %FN4%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for Kaimal spectrum                                          >> CertTest.out
ECHO ---------------------------------------------                                               >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out
%Compare% %FN4% TstFiles\%FN4%                                                                   >> CertTest.out

ECHO CompareFFtxtFiles( '%FN2%', 'TstFiles\%FN2%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN3%', 'TstFiles\%FN3%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN4%', 'TstFiles\%FN4%' )                                              >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=2
SET ROOTNAME=vonKarm

SET FN1=%ROOTNAME%.sum
SET FN2=%ROOTNAME%.dat
SET FN3=%ROOTNAME%.hh

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for von Karman spectrum.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for von Karman spectrum                                      >> CertTest.out
ECHO -------------------------------------------------                                           >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out

ECHO CompareHHfiles( '%FN2%', 'TstFiles\%FN2%' )                                                 >> CertTest.m
ECHO CompareHHfiles( '%FN3%', 'TstFiles\%FN3%' )                                                 >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=3
SET ROOTNAME=Kaimal_15

SET FN1=%ROOTNAME%.sum
SET FN2=%ROOTNAME%.u
SET FN3=%ROOTNAME%.v
SET FN4=%ROOTNAME%.w

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for Kaimal spectrum with user-specified TI.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR
IF NOT EXIST %FN4%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for Kaimal spectrum with user-specified TI                   >> CertTest.out
ECHO --------------------------------------------------------------------                        >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out
%Compare% %FN4% TstFiles\%FN4%                                                                   >> CertTest.out

ECHO CompareFFtxtFiles( '%FN2%', 'TstFiles\%FN2%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN3%', 'TstFiles\%FN3%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN4%', 'TstFiles\%FN4%' )                                              >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=4
SET ROOTNAME=vonKarm_15

SET FN1=%ROOTNAME%.sum
SET FN2=%ROOTNAME%.dat
SET FN3=%ROOTNAME%.hh

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for von Karman spectrum with user-specified TI.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for von Karman spectrum with user-specified TI               >> CertTest.out
ECHO ------------------------------------------------------------------------                    >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out

ECHO CompareHHfiles( '%FN2%', 'TstFiles\%FN2%' )                                                 >> CertTest.m
ECHO CompareHHfiles( '%FN3%', 'TstFiles\%FN3%' )                                                 >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=5
SET ROOTNAME=smooth

SET FN1=%ROOTNAME%.sum
SET FN2=%ROOTNAME%.cts
SET FN3=%ROOTNAME%.bin

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for Smooth spectrum with coherent events.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for Smooth spectrum with coherent events                     >> CertTest.out
ECHO ------------------------------------------------------------------                          >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out

ECHO CompareCTSfiles('%FN2%', 'TstFiles\%FN2%' )                                                 >> CertTest.m
ECHO CompareHHfiles( '%FN3%', 'TstFiles\%FN3%' )                                                 >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=6
SET ROOTNAME=KHtest

SET FN1=%ROOTNAME%.sum
SET FN2=%ROOTNAME%.cts

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for NWTCUP spectrum with KH-test.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for NWTCUP spectrum with KH-test                             >> CertTest.out
ECHO ----------------------------------------------------------                                  >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out

ECHO CompareCTSfiles('%FN2%', 'TstFiles\%FN2%' )                                                 >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=7
SET ROOTNAME=GPLLJ

SET FN1=%ROOTNAME%.sum
SET FN2=%ROOTNAME%.bts
SET FN3=%ROOTNAME%.cts
SET FN4=%ROOTNAME%.dat

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for GP_LLJ spectrum.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR
IF NOT EXIST %FN4%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for GP_LLJ spectrum                                          >> CertTest.out
ECHO ----------------------------------------------------------                                  >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out
%Compare% %FN4% TstFiles\%FN4%                                                                   >> CertTest.out

ECHO CompareBTSfiles( '%FN2%', 'TstFiles\%FN2%' )                                                >> CertTest.m
ECHO CompareCTSfiles( '%FN3%', 'TstFiles\%FN3%' )                                                >> CertTest.m
ECHO CompareHHfiles(  '%FN4%', 'TstFiles\%FN4%' )                                                >> CertTest.m

REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=8
SET ROOTNAME=HYDRO_TIDAL

@SET FN1=%ROOTNAME%.sum
@SET FN2=%ROOTNAME%.dat
@SET FN3=%ROOTNAME%.hh
@SET FN4=%ROOTNAME%.wnd

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for Hydro-TurbSim TIDAL spectral model.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR
IF NOT EXIST %FN4%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for Hydro-TurbSim TIDAL spectral model                       >> CertTest.out
ECHO -------------------------------------------------------------------                         >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out

ECHO CompareHHfiles(  '%FN2%', 'TstFiles\%FN2%' )                                                >> CertTest.m
ECHO CompareHHfiles(  '%FN3%', 'TstFiles\%FN3%' )                                                >> CertTest.m
ECHO CompareWNDfiles( '%FN4%', 'TstFiles\%FN4%' )                                                >> CertTest.m


IF /I NOT "%2"=="LARGE"  GOTO CompareMatlabResults


REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=9
SET ROOTNAME=GPLLJ_Large

@SET FN1=%ROOTNAME%.sum
@SET FN2=%ROOTNAME%.u
@SET FN3=%ROOTNAME%.v
@SET FN4=%ROOTNAME%.w
@SET FN5=%ROOTNAME%.cts
@SET FN6=%ROOTNAME%.wnd

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for GP_LLJ spectrum with larger grid.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR
IF NOT EXIST %FN4%  GOTO ERROR
IF NOT EXIST %FN5%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for GP_LLJ spectrum with larger grid                         >> CertTest.out
ECHO -------------------------------------------------------------------                         >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out
%Compare% %FN4% TstFiles\%FN4%                                                                   >> CertTest.out
%Compare% %FN5% TstFiles\%FN5%                                                                   >> CertTest.out

ECHO CompareFFtxtFiles( '%FN2%', 'TstFiles\%FN2%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN3%', 'TstFiles\%FN3%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN4%', 'TstFiles\%FN4%' )                                              >> CertTest.m
ECHO CompareCTSfiles(   '%FN5%', 'TstFiles\%FN5%' )                                              >> CertTest.m
ECHO CompareWNDfiles(   '%FN6%', 'TstFiles\%FN6%' )                                              >> CertTest.m


REM ================================================================================================================================
REM ================================================================================================================================
SET TEST=10
SET ROOTNAME=IECKAI_Large

@SET FN1=%ROOTNAME%.sum
@SET FN2=%ROOTNAME%.u
@SET FN3=%ROOTNAME%.v
@SET FN4=%ROOTNAME%.w
@SET FN5=%ROOTNAME%.wnd

ECHO -------------------------------------------------------------------------
ECHO TurbSim Test #%TEST%:  Validation test for IECKAI spectrum with larger grid.
ECHO -------------------------------------------------------------------------
%PROGRAM% %ROOTNAME%.inp

IF ERRORLEVEL 1  GOTO ERROR

IF NOT EXIST %FN1%  GOTO ERROR
IF NOT EXIST %FN2%  GOTO ERROR
IF NOT EXIST %FN3%  GOTO ERROR
IF NOT EXIST %FN4%  GOTO ERROR

ECHO.                                                                                            >> CertTest.out
ECHO ########################################################################################### >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO Test #%TEST%:  Validation test for Kaimal spectrum with larger grid                         >> CertTest.out
ECHO -------------------------------------------------------------------                         >> CertTest.out
ECHO.                                                                                            >> CertTest.out
ECHO ------------------------------------------------------------------------------------------- >> CertTest.out
%Compare% %FN1% TstFiles\%FN1%                                                                   >> CertTest.out
%Compare% %FN2% TstFiles\%FN2%                                                                   >> CertTest.out
%Compare% %FN3% TstFiles\%FN3%                                                                   >> CertTest.out
%Compare% %FN4% TstFiles\%FN4%                                                                   >> CertTest.out

ECHO CompareFFtxtFiles( '%FN2%', 'TstFiles\%FN2%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN3%', 'TstFiles\%FN3%' )                                              >> CertTest.m
ECHO CompareFFtxtFiles( '%FN4%', 'TstFiles\%FN4%' )                                              >> CertTest.m
ECHO CompareWNDfiles(   '%FN5%', 'TstFiles\%FN5%' )                                              >> CertTest.m


REM ================================================================================================================================
:CompareMatlabResults
REM ================================================================================================================================
REM  We're finished with the tests so now let's look at the comparisons.
REM ================================================================================================================================

ECHO -------------------------------------------------------------------------
ECHO Comparing CertTest results.
ECHO -------------------------------------------------------------------------

REM ================================================================================================================================
REM IF NOT "%1"=="MATLAB" ( GOTO CompareResults )
IF /I NOT "%1"=="MATLAB" ( GOTO CompareResults )

REM ================================================================================================================================
REM  Let's look at MATLAB. The CertTest.m file writes to the Matlab Command Window & the information generated in the command 
REM  window is also written in the logfile CertTest.log.
REM ================================================================================================================================
REM SET EnterKey=^<Enter^>
SET EnterKey="Enter"

ECHO fprintf( '\n' )                                                                            >> CertTest.m
ECHO fprintf('Type "quit" to close MATLAB and return to the CertTest if necessary.\n');         >> CertTest.m 
REM !BONNIE: THIS DOESN'T WORK WITH MATLAB R2006B!  ECHO dum = input('Press %EnterKey% to close MATLAB and return to the CertTest.');               >> CertTest.m 
REM !BONNIE: THIS DOESN'T WORK WITH MATLAB R2006B!  ECHO quit                                                                                       >> CertTest.m 

ECHO.
ECHO Opening MatLab CertTest.m....
REM !BONNIE: THIS DOESN'T WORK WITH MATLAB R2006B!  ECHO   Press %EnterKey% in the MATLAB command window or close MATLAB to continue CertTest. 
ECHO   You may need to close MATLAB to continue the CertTest.
%MATLAB% /r CertTest /nosplash /logfile CertTest.log


REM ================================================================================================================================
:CompareResults
REM ================================================================================================================================
REM  Let's look at the CertTest.out file in the Editor.
REM ================================================================================================================================
ECHO.
ECHO Opening CertTest.out....
ECHO   Close the CertTest.out editor to finish this CertTest.

%Editor% CertTest.out

GOTO END


REM ================================================================================================================================
:ERROR
REM ================================================================================================================================
ECHO ** An error has occured in Test #%TEST% **


REM ================================================================================================================================
:END
REM ================================================================================================================================
TYPE Bell.txt
ECHO Processing complete.


SET Compare=
SET DateTime=
SET Editor=
SET PROGRAM=
SET MATLAB=
SET Test=
SET ROOTNAME=
SET EnterKey=

SET FN1=
SET FN2=
SET FN3=
SET FN4=
SET FN5=
SET FN6=

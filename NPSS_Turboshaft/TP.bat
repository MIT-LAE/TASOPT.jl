@echo off
@SET _calldir=%CD%
@SET EMI_TOP="C:\Users\Prashanth\NPSS\NPSS.EMI.14.0.1_LE"
@SET EMI_TOP=%EMI_TOP:"=%
@Set NPSS_VERSION_STRING=NPSS_3.2_nt_Release_VC14_64
@set NPSS_CONFIG=nt
@set NPSS_TOP=C:\Users\Prashanth\NPSS\NPSS.nt.ver32_VC14_64_LE
@set NPSS_DEV_TOP=C:\Users\Prashanth\NPSS\NPSS.nt.ver32_VC14_64_LE
@set DCLOD_PATH=%NPSS_TOP%\DLMComponents\nt
@set PATH=%PATH%;%NPSS_TOP%\bin;%NPSS_TOP%\scripts
@SET PATH=%EMI_TOP%\scripts\nt;%EMI_TOP%\scripts\AutoDoc;%PATH%
@cd %_calldir%
@SET _calldir=


REM Run file
runnpss TShaft.run -I ./src -I ./model -I ./data

PAUSE

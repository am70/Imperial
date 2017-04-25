@echo off
REM automatically generated
ECHO generated on host: WPIA-DIDE142
ECHO generated on date: 2017-04-25
ECHO didehpc version: 0.1.0
ECHO context version: 0.1.0
ECHO running on: %COMPUTERNAME%
set CONTEXT_WORKDRIVE=q:
set CONTEXT_WORKDIR=ImperialMalaria\larvalModel\packages\odinPackage
set CONTEXT_ROOT=q:\ImperialMalaria\larvalModel\packages\odinPackage\contexts
set CONTEXT_ID=ece2816d8acb6791bd876684b5ec1d82
set CONTEXT_PROPAGATE_ERROR=TRUE
set CONTEXT_BOOTSTRAP=TRUE
call setr64_3_3_2.bat
ECHO mapping Q: -^> \\fi--san02\homes\alm210
net use Q: \\fi--san02\homes\alm210 /y
ECHO mapping T: -^> \\fi--didef2\tmp
net use T: \\fi--didef2\tmp /y
ECHO mapping q: -^> \\fi--san03\homes\alm210
net use q: \\fi--san03\homes\alm210 /y
ECHO Using Rtools at T:\Rtools\Rtools33
set PATH=T:\Rtools\Rtools33\bin;T:\Rtools\Rtools33\gcc-4.6.3\bin;%PATH%
set BINPREF=T:/Rtools/Rtools33/mingw_64/bin/
ECHO This is a parallel job: will use %CPP_NUMCPUS%
set CONTEXT_CORES=%CCP_NUMCPUS%
set REDIS_HOST=12.0.0.1
set REDIS_URL=redis://12.0.0.1:6379
%CONTEXT_WORKDRIVE%
cd \%CONTEXT_WORKDIR%
ECHO working directory: %CD%
ECHO this is a single task
set CONTEXT_TASK_ID=41e8f4bc816144776ea8737a330c8739
set CONTEXT_LOGFILE=q:\ImperialMalaria\larvalModel\packages\odinPackage\contexts\logs\%CONTEXT_TASK_ID%
ECHO logfile: %CONTEXT_LOGFILE%
@REM The quoting here is necessary for paths with spaces.
ECHO on
Rscript "q:\ImperialMalaria\larvalModel\packages\odinPackage\contexts\bin\task_run" "%CONTEXT_ROOT%" %CONTEXT_TASK_ID% > "%CONTEXT_LOGFILE%" 2>&1
@ECHO off
if %ERRORLEVEL% neq 0 (
  ECHO Error running task
  EXIT /b %ERRORLEVEL%
)
@ECHO Quitting

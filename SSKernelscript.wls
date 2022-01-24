#!/usr/bin/env wolframscript
(* ::Package:: *)

workdir=If[$InputFileName=="",NotebookDirectory[],DirectoryName@$InputFileName];
Print["Working directory is ",workdir]

SetDirectory[workdir<>"\\Packages"];
Get["Configuration`"]
Get["ModelProcessing`"]
Get["Centering`"]
Get["RedundancyRemoval`"]
Get["FixedFluxes`"]
Get["RayFinding`"]
Get["FacetTesting`"]
Get["FBFsearch`"]
Get["Capping`"]
ResetDirectory[];
PackagesDirectory=workdir<>"\\Packages";

datafile=First@FileNames["*.mat"|"*.m", workdir<>"TEST"];
commandline=True;  verbose=False; Species="";
(* Print["There are "<>ToString@Length@$ScriptCommandLine<>" script arguments."]; *)
If[Length[$ScriptCommandLine]>=2, 
	datafile=DataDirectory<>$ScriptCommandLine[[2]](*;
	Print[">",datafile,"<"]*)];
If[Length[$ScriptCommandLine]==3, commands=ToExpression[$ScriptCommandLine[[3]]];
Print["A total of "<>ToString@Length[commands]<>
" parameter values are adjusted by the command line."];];
If[!StringQ[datafile]||StringLength[datafile]==0,
	Print["The file path is empty! If you use SSKernelbatch.bat, drag an input file to it."];
	Input["Press x and <Enter> to exit"];Exit[]]
If[FileExistsQ[datafile], Print["I found the data file, allright!"],
Print["No file was found at the path "<>datafile];Input["Press x and <Enter> to exit"];Exit[]];

(* STAGE 1 *)
success=True;
If[StringLength[Species]<1,  Species =FileNameTake[datafile, {-2}]]; 
result= CheckAbort[LoadModel[datafile, verbose], success = False];
If[Length[feasiblepoints] < 1, success = False;
     Print["LP solution of the imported FBA model has failed. Retry loading after
 adjusting the LP tolerance, or reversing max/min spec in the datafile.  "],
  	 {S, bounds, objectselector, maxmin} = result];
Heading = {"STAGE 1: Loaded data for " <> Species, "from " <> FileNameTake[datafile] <>
    	 " with LP test in " <> Quiet@ToString[NumberForm[Testime, 3]] <> " seconds."};
If[success,Print[Heading],Print["Load and test of the datafile has failed."];
			Input["Press x and <Enter> to exit"];Exit[]];    	 

(* STAGE 2*)
SStype = "FacetCone";
Reductiontable = {{"Stage", " Constraints ", " Variables ", " Ray Yield "}};
{RSStime, result} = Timing[
	result = CheckAbort[Check[
   		ReducedSolutionSpace[S, bounds, objectselector, maxmin, verbose],
   		 success = False],  success = False];
	If[StringQ[result], success=False; Print[result]; PrependTo[feasiblepoints,{}]]];
 Heading = {"STAGE2: Reduced Solution Space for " <> Species, 
   	"  Model " <> FileBaseName[datafile] <> " of type " <> ToString@SStype <> 
   	"\n calculated on " <> DateString["DateShort"] <> " in " <> 
      Quiet@ToString[NumberForm[RSStime, 3]] <> " seconds."};
If[success,Print[Heading],Print["RSS reduction of the solutionspace has failed."];
		SaveResults[{}, datafile];
		Input["Press x and <Enter> to exit"];Exit[]];    	 

(* STAGE 3 *)
{kerneltime, result} =  Timing[CheckAbort[Check[
    	KernelFinder[samplesize, treesize, fixtol, verbose],
    	 success = False], success = False]];
 Heading = {"STAGE 3: Capped Solution Space for " <> Species, 
     " Model " <> FileBaseName[datafile] <> " of uncapped type " <> ToString@SStype <>
       "\n calculated on " <> DateString["DateShort"] <> " in " <> 
     Quiet@ToString[NumberForm[kerneltime, 3]] <> " seconds."};
If[success,
	SSKinRSS = KernelSpace; 
    SSKinsphere = InscribedSphere;
	Print[Heading],
	
	Print["Calculation of a closed SSK has failed."];
		SaveResults[{}, datafile];
		Input["Press x and <Enter> to exit"];Exit[]];    	 

(* STAGE 4 *)
If[SStype != "Point",
    {Shapetime, result} = Timing[CheckAbort[Check[KernelShape[maxaspect, maxthin, verbose],
    	 success = False], 	 success = False]]; 
    Heading = {"STAGE 4: Shape Analysis for " <> Species, 
      " Model " <> FileBaseName[datafile] <> "  with uncapped kernel of type " <> 
      ToString@SStype <> 
       "\n calculated on " <> DateString["DateShort"] <> " in " <> 
       Quiet@ToString[NumberForm[Shapetime, 3]] <> " seconds."};,
    Shapetime = 0; 
    flatplot = Graphics@Text["Point SSK, has no shape."]];
If[success, Print[Heading],Print["Shape analysis of the closed SSK has failed."];
		SaveResults[{}, datafile];
		Input["Press x and <Enter> to exit"];Exit[]];    	 

(* STAGE 5 *)
{Chordtime, result} = Timing[CheckAbort[Check[
    	{chordpic, cenpic} = KernelDisplay[ExportResults, verbose],
    	 success = False],  success = False]];
 Heading = {"STAGE 5: Full SS Kernel Analysis for " <> Species, 
     " Model " <> FileBaseName[datafile] <> "  with uncapped kernel of type " <> 
      ToString@SStype <> "\n calculated on " <> DateString["DateShort"] <> " in " <> 
      Quiet@ToString[NumberForm[Testime + RSStime + kerneltime + Shapetime + Chordtime, 3]] <> 
      " seconds."}; 
If[success, Print[Heading],Print["Display of chords and centering has failed."];
		SaveResults[{}, datafile];
		Input["Press x and <Enter> to exit"];Exit[]];    
		
(* STAGE 6 *)
(* Only export first 8 components of exportresults; component 9
  	  is the RSS structure, not really externally useful. *)
Heading = {"Full SS Kernel Analysis for " <> Species, 
     " Model " <> FileBaseName[datafile] <> "  with uncapped kernel of type " <> 
      ToString@SStype <> "\n calculated on " <> DateString["DateShort"] <> " in " <> 
      Quiet@ToString[NumberForm[Testime + RSStime + kerneltime + Shapetime + Chordtime, 3]] <> 
      " seconds."}; 
CheckAbort[success=Check[
   	SaveResults[If[Length@ExportResults==0,{},ExportResults[[1;;8]]], datafile],
   	 False],  success = False];

If[success, Print["Output files have been stored successfully."],
	Print["Export of output files has failed."]; 
	Input["Press x and <Enter> to exit"];Exit[]];  

Input["Press x and <Enter> to exit"]; 




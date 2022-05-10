(* Wolfram Language Package *)

(* COPYRIGHT
					© Copyright 2022 Wynand Verwoerd

This file is part of SSKernel.

The SSKernel program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE . See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with this program . If not, see http://www.gnu.org/licenses/
*)

BeginPackage["DataHandler`",{"Configuration`"}]
(* Exported symbols added here with SymbolName::usage *)  
MReader::usage = "Reads input data {modelname, S, bounds, objectvector, FBAvector, ReactionNames, maxmin}
 from a Mathematica .m file" 
MatReader::usage="Reads input data from a MatLab .mat file, using COBRA spec as in  https://arxiv.org/pdf/1710.04038.pdf  "


Begin["`Private`"] (* Begin Private Context *) 


SetAttributes[MReader, HoldRest];
MReader[file_, S_, Svals_, bounds_, Object_, FBAvector_, 
	 ReactionNames_, MetabNames_, maxmin_, printresult_:False] := 
(* This function reads input data from a Mathematica .m file that contains the numeric 
arrays as a list of the structure 
		{modelname, S, Svals(Optional), bounds, objectvector, FBAvector, ReactionNames, MetabNames, maxmin}. 
Here modelname is a string, which may include a formal ID and a desciption of the model. It has to appear first.
The remaining items may appear in arbitrary order, and any additional list items are ignored.
maxmin is not an array, but a string containing either "max" or "min"
to indicate whether the objective should be minimized or maximized. 
It is optional, the default is taken as maximisation. 

From within Mathematica, having created the individual arrays somehow, produce the file  by executing: 
	Export[exportfile, {ModelName, S, Svals, bounds, objectselector, FBAvector,  ReactionNames, MetaboliteNames,
  Switch[maxmin, -1, "max", 1, "min", _, "max"]}, "List"];

*)
 Block[{stringimport, eximport, import, indims, candidates, foundvals, Sat, valsat, 
 	boundsat, cons, vars, vectorat,	namesat, FBAvecat, selectorat,SvalsTest},
 	(* Test whether a vector satisfies all constraint equations as determined by Svals *)
  SvalsTest[svals_, vector_] := Block[{equals, smallers, biggers},
  	equals = Position[svals[[All, 2]], 0];
  	smallers = Position[svals[[All, 2]], -1];
  	biggers = Position[svals[[All, 2]], 1];
  	And @@ Join[
    MapThread[Equal, {Extract[vector, equals], Extract[svals[[All, 1]], equals]}], 
    MapThread[LessEqual, {Extract[vector, smallers], Extract[svals[[All, 1]], smallers]}], 
    MapThread[GreaterEqual, {Extract[vector, biggers], Extract[svals[[All, 1]], biggers]}]]
  ];
  Print["Data imported from file " <> file];
(* Using the "List" format takes care of importing numbers written e.g. as 1.23E-10 *)
  eximport = Import[file, "List"];
  ModelName=First@eximport;
  eximport= ToExpression /@ Rest[eximport];
  (* Separate off any non-array, stemming from a string, and convert it back *)
  import=Select[eximport, ArrayQ];
  stringimport=ToString /@ Complement[eximport,import];
  stringimport = If[stringimport == {}, "max", stringimport[[1]]];
  maxmin=Which[StringContainsQ[stringimport, "max", IgnoreCase -> True], -1, 
 	StringContainsQ[stringimport, "min", IgnoreCase -> True], 1, True, -1];
  (* The imported arrays may be in arbitrary order and have 
	superfluous data, such as the FBA optimum value, interspersed. 
  	So identify the required ones by their dimensions and data type. *)
  
  indims = Dimensions /@ import;
  If[printresult, Print["The arrays as imported have the following dimensions \n"<>  ToString[indims]]]; 
  Sat = FirstPosition[indims, {_, _?(# > 2 &)}, {0}][[1]];
  If[Sat > 0, S = import[[Sat]]];
  If[Sat==0 || !MatrixQ[S, NumberQ], 
  	Print[Style["Abort - No 2D numeric array with more than 2 rows could be identified 
   as a stoichiometry matrix in the input data.",Red]]];
  {cons, vars} =indims[[Sat]];
  (* Find arrays with the right dimensions to contain values and bounds *)
  valsat = Position[indims, {cons, 2}, {1}];
  boundsat = Position[indims, {vars, 2}, {1}];
  (* Multiple cases of each may be found, in particular if vars = cons. 
	Isolate the values array by checking that it has pairs with the second element -1, 0, or 1 *)
  candidates = If[valsat != {}, Extract[import, valsat], {}];
  foundvals = FirstPosition[candidates, 
  	{Repeated[{_?NumberQ, _?(Sign[#] == # &)}, {cons}]}, {}];
  Svals = If[foundvals != {}, Extract[candidates, foundvals], {}];
  If[Svals == {}, Svals = ConstantArray[{0, 0}, cons];
  	If[printresult, Print["Default: constraints taken as = 0"]]];
  (* Remove values array, then check remaining bounds candidates for 
	the first one having low-high pairs *)
  boundsat = Complement[boundsat, {Extract[valsat, foundvals]}];
  candidates = If[boundsat != {}, Extract[import, boundsat], {}];
  foundvals = FirstPosition[candidates,
  	 _?(Count[#, {lo_, hi_} /; hi >= lo] == vars &), {}];
  bounds = If[foundvals != {}, Extract[candidates, foundvals], {}];  
  If[bounds == {},
      Print[Style["Abort - A 2D numeric array with "<>ToString[vars]<>
      " rows, each a pair of low and high numbers, to serve as bounds, "<>
      "could not be be identified in the input data.", Red]]; Abort[]];
  (* Next, vectors with reaction names, objective selector and possibly FBA solution *)
  vectorat = Flatten@Position[indims, {vars}]; 
  namesat=Flatten@Position[import, _?(VectorQ[#, StringQ]&)];  
  ReactionNames = Cases[import[[namesat]], _?(Length[#]==vars&)];
  ReactionNames = If[ReactionNames == {}, 
  	Print["No suitable list of reaction names or ID's were detected in the import data!"];
  	ConstantArray["NoName",vars], First@ReactionNames];
  MetabNames = Cases[import[[namesat]], _?(Length[#]==cons&)];
  MetabNames = If[MetabNames == {}, 
  	Print["No suitable list of metabolite names or ID's were detected in the import data!"];
  	ConstantArray["NoName",vars], First@MetabNames];
  vectorat=Complement[vectorat,namesat];  
(*  At this point, any numeric vectors in flux space are isolated.
	If any of those satisfies the constraint set, the first of those 
	is taken as the FBA flux; if none, it is set to a zero vector.  *) 
  FBAvecat = Cases[vectorat, _?(SvalsTest[Svals,S.import[[#]]] &)]; 
  FBAvector = If[FBAvecat == {}, ConstantArray[0., vars], import[[First@FBAvecat]]];
  selectorat = Cases[vectorat, _?(0 < Count[import[[#]], Except[0 | 0.]] < 5 &)];
  Object = Switch[Length[selectorat],
    0, Print@
     "Input data invalid - no objective vector that combines fewer 
than 5 fluxes could be identified"; ConstantArray[0., vars],
    1, import[[First@selectorat]],
    _, Print@
     "Multiple apparent objective vectors found in input data, only 
the first one is taken."; import[[First@selectorat]]];
  Print["Dimensions of imported stoichiometry, values, bounds, objective and FBA flux arrays: " <> 
    ToString[Dimensions /@ {S, Svals, bounds, Object, FBAvector}]]; 
  Print["The objective vector selects contributions from the flux set " <> 
  	ToString[Flatten@Position[Object, _?(# > 0 &)]]];
  ]

SetAttributes[MatReader, HoldRest];
MatReader[file_, S_, Svals_, bounds_, Object_,  FBAvector_,
	 ReactionNames_, MetabNames_, maxmin_, printresult_:False] := 
(* This function reads input data from a MatLab .mat file that specifies a metabolic model
according to the COBRA specification (See https://arxiv.org/pdf/1710.04038.pdf). 
The data item "osenseSTR" is optional and if absent, the objective is maximized. *)
 Block[{labels, modelname, description, missing, ruleform, reactionIDs, reactionnames, 
 	reactionECnos, metIDs, metnames, lowers, uppers, ineq, vals, cons, vars},
  maxmin="max";
  Print["Data imported from file " <> file];
  If[printresult,Print@Import[file, "Comments"]];
  (* The following used to work with Mathematica versions up to 12.1*)
(*
  labels = Import[file, "LabeledData"][[1, 2, All, 1]];
  ruleform = First[Import[file, "Labels"] /. Import[file, "LabeledData"]]; 
*)
  (* But since 12.2, it imports as a set of associations. Convert to rules by Normal. *)
  ruleform = Normal[Import[file, "LabeledData"]][[1, 2, 1, 1, All]];
  labels = ruleform[[All, 1]];
  If[printresult,Print["The labeled data items in the input file are: "<>ToString[labels]]];
  missing = Complement[{"S", "lb", "ub", "b", "c"}, labels];
  If[missing!={},
  	Print[Style["Imported data file is incomplete; the following essential data was not found:"<>
  	ToString[missing /. {"S" -> "Stoichiometry matrix", "lb" -> "Lower bounds", 
  		"ub" -> "Lower bounds", "b" -> "S-constraint values", "c" -> "Objective"}],Red]];
  Abort[]];
  {modelname, description, S, lowers, uppers, Object, ineq, vals, maxmin, 
  	reactionIDs, reactionnames, reactionECnos, metIDs, metnames} = 
  	{"modelID", "description", "S", "lb", "ub", "c", "csense", "b","osenseStr",
  		"rxns","rxnNames", "rxnECNumbers","mets", "metNames"} /. ruleform;
  	(* Remove superfluous braces that result from the rules *)
  {modelname, description, ineq, maxmin} = {modelname, description, ineq, maxmin} /. {string_} :> string;
  {reactionIDs, reactionnames, reactionECnos} = 
  		Replace[{reactionIDs, reactionnames, reactionECnos}, t_ /; ArrayQ[t] :> Flatten[t],{1}];
  {metIDs, metnames} = Replace[{metIDs, metnames}, t_ /; ArrayQ[t] :> Flatten[t], {1}];
  modelname= If[modelname=="modelID", FileNameTake[file], modelname];
  description=If[description == "description", "", " : "<>description];
  ModelName = Species<>" - " <> modelname <> description;
  maxmin=Which[StringContainsQ[maxmin, "max", IgnoreCase -> True], -1, 
 	StringContainsQ[maxmin, "min", IgnoreCase -> True], 1, True, -1];
  {cons, vars} =Dimensions@S; 
  
  (* The following gives precedence to reaction IDs/metabolite IDs; 
  		simply rearrange ordering if e.g. full names are preferred *)
  ReactionNames = Which[Length@reactionIDs > 0, Flatten@reactionIDs, Length@reactionECnos > 0, Flatten@reactionECnos,
  	Length@reactionnames > 0, Flatten@reactionnames, True, ConstantArray["NoName",vars]  ];
  MetabNames = Which[ Length@metIDs > 0, Flatten@metIDs, Length@metnames > 0, Flatten@metnames,
  	 True, ConstantArray["NoName",cons]  ];

  bounds = Transpose[Flatten /@ {lowers, uppers}];
  Object = Transpose@Object; 
  vals = First@Transpose@vals;
  (* If the csense item is missing, assume all constraints to be equals 
  		Also, provide for the case that ineq comes through as a single string 
  		of characters instead of a list of characters *)
(* Print[{"before processing, ineq =",ineq}]; *)
  ineq=If[ineq=="csense", ConstantArray[0,Length[vals]],
  				If[StringQ[ineq], StringPartition[ineq, 1]]];
  ineq = ineq /. {"E" -> 0, "L" -> -1, "G" -> 1};
(* Print[{"after processing, ineq =",ineq}]; *)
  Svals = Transpose[{vals, ineq}];
  If[Svals != ConstantArray[{0, 0}, Length[S]], 
   Print["WARNING: The loaded model contains inequalities and/or 
nonzero RHS vector. The SSK calculation currently does not provide for that! "]];
  If[printresult,Print["Dimensions of imported stoichiometry, bounds and objective arrays: " <> 
  	ToString[Dimensions /@ {S, bounds, Object}]]];
  objectcount=First@Dimensions[Object];
  If[objectcount==1,Object=First@Object;
  If[printresult,Print["The objective vector selects contributions from the flux set " <> 
  	ToString[Flatten@Position[Object, _?(Abs[#] > 0 &)]]]],
  Print["The objective specification contains "<>ToString[objectcount]<>" objective vectors."]]; 
  If[printresult,Print["This objective is to be "<>If[maxmin==1,"minimized.","maximized."]]];

  FBAvector = ConstantArray[0., Length[Object]];
  ]
  
 


End[] (* End Private Context *)

EndPackage[]
(* Targeted comparison with unmodified public StableWein computational cells.
   Native Wolfram Language kernel required. No unitary/NN filtering is applied.
   The pinned notebook is fetched, not redistributed in this research bundle. *)
base=DirectoryName[$InputFileName]; commit="a540c6dbe099998f53d5af098e1c4a325c54ca84";
text=Import["https://raw.githubusercontent.com/jurciukonis/StableWein/"<>commit<>"/StableWein.nb","Text"];
If[Hash[text,"SHA256","HexString"]!="acf62cd0c9c4de76f2cef5e54b127071c085cf50f2cf5dcecac1c2aa91e8d098",Print["Notebook content hash mismatch"];Exit[2]];
nb=ToExpression[text,InputForm,HoldComplete]; boxes=Cases[nb,Cell[BoxData[b_],"Input",___]:>b,Infinity];
wpCompTargetPac="WVM"; wbLoadBinariesPac=False;
held=ToExpression[#,StandardForm,HoldComplete]&/@boxes[[{33,34,35,36,37,38,39,40,42,43,44,46}]];
If[!FreeQ[held,$Failed],Print["Selected source cell failed to parse"];Exit[2]];Scan[ReleaseHold,Flatten[held]];
grid=wbWeinGridGen[{17,40,0}]; exact[x_]:=If[StringQ[x],ToExpression[x,InputForm],x];
mapInput[data_]:=Module[{a,b,c,dr,di,d,o={2,3,1}},{a,b,c,dr,di}=(exact/@data[#]&)/@{"a","b","c","d_real","d_imag"};d=dr+I di;Join[a,b[[o]],c[[o]],Abs[d[[o]]],{If[Times@@d==0,0,Arg[Times@@d]]}]];
names={"A"->"regression_A","B"->"regression_B","Bminus"->"regression_Bminus","C"->"regression_C","charged_flat"->"charged_only_flat","charged_positive"->"charged_only_positive","charged_negative"->"charged_only_negative","weighted_flat"->"weighted_charged_flat","real_projector_flat"->"real_projector_flat","singular_flat"->"singular_flat","all_zero"->"all_zero"};
inputs=Association@Table[First[r]->mapInput[Import[FileNameJoin[{base,Last[r]<>"_input.json"}],"RawJSON"]],{r,names}];
rows=Table[With[{p=N[inputs[name]]},<|"case"->name,"mapped_input"->p,"mode0_sufficient"->wbWeinSufCond[p],"mode1"->wbWeinNecCond1[p],"mode2"->(wbWeinNecCond1[p]&&wbWeinNecCond2[p]),"mode3"->(wbWeinNecCond1[p]&&wbWeinNecCond2[p]&&wbWeinNecCond3[p,Sequence@@grid])|>],{name,Keys[inputs]}];
mode4={};Do[SeedRandom[20261021+batch];Do[timed=AbsoluteTiming[TimeConstrained[wbMinimFuncPotWein[N[inputs[name],70],{1,2,5,10}],16,"TIMEOUT"]];AppendTo[mode4,<|"case"->name,"seconds"->First[timed],"result"->Last[timed],"seed"->20261021+batch|>],{name,{{"Bminus","C","B"},{"A","charged_positive","charged_negative"}}[[batch]]}],{batch,2}];
result=<|"kernel_version"->$Version,"commit"->commit,"grid_parameters"->{17,40,0},"grid_sizes"->(Length/@grid),"results"->rows,"mode4"->mode4|>;
Export[FileNameJoin[{base,"stablewein_rerun_results.json"}],result,"RawJSON"];Print[result];

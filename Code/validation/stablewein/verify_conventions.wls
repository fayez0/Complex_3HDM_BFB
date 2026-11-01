(* Exact polynomial check of the potential normalization and oriented edges. *)
nb=ToExpression[Import["https://raw.githubusercontent.com/jurciukonis/StableWein/a540c6dbe099998f53d5af098e1c4a325c54ca84/StableWein.nb","Text"],InputForm,HoldComplete];
bs=Cases[nb,Cell[BoxData[b_],"Input",___]:>b,Infinity];Scan[ReleaseHold,Flatten[ToExpression[#,StandardForm,HoldComplete]&/@bs[[{37,38,39}]]]];
pars={a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,ep};vars={1,0,cr,ci,dr,di,er,ei,fr,fi};xx={{0,cr+I ci,er+I ei},{1,dr+I di,fr+I fi}};kk=ComplexExpand[ConjugateTranspose[xx].xx];edges={{1,2},{2,3},{3,1}};
results=Table[d={d3,d1,d2};d[[{2,3,1}[[j]]]]=d[[{2,3,1}[[j]]]]Exp[I ep];pot=Sum[{a1,a2,a3}[[i]]kk[[i,i]]^2/2,{i,3}]+Sum[With[{i=edges[[e,1]],k=edges[[e,2]]},{b3,b1,b2}[[e]]kk[[i,i]]kk[[k,k]]+{c3,c1,c2}[[e]]ComplexExpand[Abs[kk[[i,k]]]^2]+ComplexExpand[Re[d[[e]]kk[[i,k]]^2]]],{e,3}];Expand[ComplexExpand[{wbPotFuncWein,wbPotFuncWein2,wbPotFuncWein3}[[j]][pars,vars]-pot]]===0,{j,3}];
Print[results];If[!And@@results,Exit[1]];

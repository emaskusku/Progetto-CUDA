n=500;
prec=10;

m=Table[Random[Integer,{0,100}],{i,n},{j,n}];

str=OpenWrite["matrix.txt"];

Write[str,n];

Do[Write[str,m[[i,j]]],{i,n},{j,n}];

Close[str];

inc=Table[Subscript[x,i],{i,n}];

term=Table[Random[Integer,{0,100}],{i,n}];

str=OpenWrite["term.txt"];

Do[Write[str,term[[i]]],{i,n}];

Close[str];

sis=m.inc;

eq=Table[sis[[i]]==term[[i]],{i,n}];

sol=Solve[eq,inc];

sol=Flatten[Table[Subscript[x,i] /. sol,{i,n}]];

sol=N[sol,prec];

str=OpenWrite["solutions.txt"];

Do[Write[str,N[sol[[i]]]],{i,n}];

Close[str];

m={{8,1,0,1},{4,5,5,7},{9,11,7,5},{7,4,6,8}};


Triangolizza[matrice_List]:=Module[{matr, N},
		matr=matrice;
		N=Length[Transpose[matr]];
		matr=CG[matr,1,N];
		MatrixForm[matr]
		];


CG[matrice_List, passo_, dim_]:=Module[{matr,pos,N,col,coeff,i,j},
		matr=matrice;
		If[passo===dim,matr,
			col=Delete[matr[[All,passo]],Map[List,Range[passo-1]]];
			col=Abs[col];
			pos=Position[Abs[matr[[All,passo]]],Max[col]];
			temp=matr[[passo]];
			matr[[passo]]=matr[[pos[[1,1]]]];
			matr[[pos[[1,1]]]]=temp;

			For[i=passo, i<dim, i++,
				coeff=-matr[[i+1,passo]]/matr[[passo,passo]];
				For[j=passo, j<=dim, j++,
					matr[[i+1,j]]=matr[[i+1,j]]+coeff*matr[[passo,j]];
				];
			];

			//Print[MatrixForm[matr]];
			matr=CG[matr,passo+1,dim];
			];
			matr
		];

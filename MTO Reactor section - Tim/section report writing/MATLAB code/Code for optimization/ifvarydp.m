%% bh from 1-10m

n=0;
dpoptim=[];
dpmatrix=[];
for dp=108*10^(-6):50*10^(-6):308*10^(-6)
    
    n=n+1;
    dpmatrix(n,1)=dp;
    p=3000000;
    t=500+273.15;
    f=3.73;
    cw=60000;
    
    outfunc=newzsm0func(p,t,f,cw);
    dpoptim(n,1:26)=outfunc;
    
end
%% generating excel

[a,b]=size(dpoptim);


dpoptimcell=num2cell(dpoptim);
dpmatrixcell=num2cell(dpmatrix);
exceldpmolfra=zeros(a+1,b+1);
exceldpmolfrac=num2cell(exceldpmolfra);
exceldpmolfrac(1,:)= {'catalystweight (kg)','ch3oh','dme','h2o','c2h4','c3h6','c4h8','c5h10','c6h12','cx','c','ch4','rxtorexitTmassflow(kt/a)','plantcap','rxtordiameter','catalystweight','inletmethanolmassflow(kt/a)','WHSV','meohconversion','bedheight (m)','u0initial','c1select','c2select','c3select','c4select','c5select','c6select'};
exceldpmolfrac(2:a+1,1)=dpmatrixcell;
exceldpmolfrac(2:a+1,2:b+1)=dpoptimcell;

xlswrite('ifvaryparticlediameter',exceldpmolfrac)
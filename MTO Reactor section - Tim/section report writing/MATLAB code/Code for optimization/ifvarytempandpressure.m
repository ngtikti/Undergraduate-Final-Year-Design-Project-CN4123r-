%% temperature from 450-650oC
n=0;
toptim=[];
tmatrix=[]
for t=(400+273.15):20:(540+273.15)
    
    n=n+1;
    tmatrix(n,1)=t;
    p=3000000;
    f=3.73;
    cw=27500;
    dp=180*10^(-6);
    
    outfunc=newzsm0func(p,t,f,cw);
    toptim(n,1:26)=outfunc;
    
end

    
%% pressure from 10atm to 40atm

n=0;
poptim=[];
pmatrix=[];

for p=100000:100000:3000000
    
    n=n+1
    pmatrix(n,1)=p;
    t=550+273.15;
    f=3.73;
    cw=27500;
    dp=180*10^(-6);
    
    outfunc=newzsm0func(p,t,f,cw);
    poptim(n,1:26)=outfunc;
    
end
%% generating excel sheet

[a,b]=size(toptim);
[c,d]=size(poptim);

toptimcell=num2cell(toptim);
tmatrixcell=num2cell(tmatrix);
exceltmolfra=zeros(a+1,b+1);
exceltmolfrac=num2cell(exceltmolfra);
exceltmolfrac(1,:)= {'Temperature(K)','ch3oh','dme','h2o','c2h4','c3h6','c4h8','c5h10','c6h12','cx','c','ch4','rxtorexitTmassflow(kt/a)','plantcap','rxtordiameter','catalystweight','inletmethanolmassflow(kt/a)','WHSV','meohconversion','bedheight (m)','u0initial','c1select','c2select','c3select','c4select','c5select','c6select'};
exceltmolfrac(2:a+1,1)=tmatrixcell;
exceltmolfrac(2:a+1,2:b+1)=toptimcell;

poptimcell=num2cell(poptim);
pmatrixcell=num2cell(pmatrix);
excelpmolfra=zeros(c+1,d+1);
excelpmolfrac=num2cell(excelpmolfra);
excelpmolfrac(1,:)= {'pressure(Pa)','ch3oh','dme','h2o','c2h4','c3h6','c4h8','c5h10','c6h12','cx','c','ch4','rxtorexitTmassflow(kt/a)','plantcap','rxtordiameter','catalystweight','inletmethanolmassflow(kt/a)','WHSV','meohconversion','bedheight (m)','u0initial','c1select','c2select','c3select','c4select','c5select','c6select'};
excelpmolfrac(2:c+1,1)=pmatrixcell;
excelpmolfrac(2:c+1,2:d+1)=poptimcell;

xlswrite('ifvarytemp',exceltmolfrac)
xlswrite('ifvarypressure',excelpmolfrac)

    

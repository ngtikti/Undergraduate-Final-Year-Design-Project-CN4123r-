%% bh from 1-10m

n=0;
cwoptim=[];
cwmatrix=[];
for cw=1000:200:4000
    
    n=n+1;
    cwmatrix(n,1)=cw;
    p=3000000;
    t=550+273.15;
    f=3.73;
    dp=180*10^(-6);
    
    outfunc=newzsm0func(p,t,f,cw);
    cwoptim(n,1:26)=outfunc;
    
end
%% generating excel

[a,b]=size(cwoptim);


bhoptimcell=num2cell(cwoptim);
bhmatrixcell=num2cell(cwmatrix);
excelbhmolfra=zeros(a+1,b+1);
excelbhmolfrac=num2cell(excelbhmolfra);
excelbhmolfrac(1,:)= {'catalystweight (kg)','ch3oh','dme','h2o','c2h4','c3h6','c4h8','c5h10','c6h12','cx','c','ch4','rxtorexitTmassflow(kt/a)','plantcap','rxtordiameter','catalystweight','inletmethanolmassflow(kt/a)','WHSV','meohconversion','bedheight (m)','u0initial','c1select','c2select','c3select','c4select','c5select','c6select'};
excelbhmolfrac(2:a+1,1)=bhmatrixcell;
excelbhmolfrac(2:a+1,2:b+1)=bhoptimcell;

xlswrite('ifvarycatweight',excelbhmolfrac)
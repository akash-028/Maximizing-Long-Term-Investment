clear ;
clc ;
%%
prompt = {'Amount that you want to Invest : ','Time Period : ','Number of Bonds : '} ;
dlgtitle = 'Input' ;
dims = [1 35];
definput = {'1000','5','4'} ;
answer = inputdlg(prompt,dlgtitle,dims,definput) ;
answer = str2double(answer) ;
Initial = answer(1) ;   % Initial Amount to be Invested
T = answer(2) ;         % Time Period for Invest this Amount
N = answer(3) ;         % Available Bonds for Investing
total_bonds = N + T ;   % Total Bonds containing T bonds with 0 Interest Rates

bond = zeros(total_bonds,1) ;   % Set Bond details
for i=1:total_bonds
    bond(i) = i ;
end

%%
% For Make a choice about Details Of bonds.
% 1. If user want to give manual details.  
% 2. If number of bonds are larger than user can randomly generate bond
%    details.
question_ans = questdlg({'1. User input bond details';...
    '2. Generate Random bond details'},'Make a choice Window','1','2','1') ;
if strcmp(question_ans,'1')
   choice = 1 ;
else
   choice = 2 ;
end
if choice==1
    % For User input
    Maturity = zeros(total_bonds,1) ;           % Maturity Period Of the Bond
    Interest_Rates = zeros(total_bonds,1) ;     % Interest Rate for the Bond
    Purchase_Years = zeros(total_bonds,1) ;     % In which year bond can be purchase
    for i=1:total_bonds
        if i<=T                     % For Extra T bonds
           Maturity(i) = 1 ; 
           Interest_Rates(i) = 0 ;
           Purchase_Years(i) = i ;
        else                        % For N given Bonds
            prompt = {'Maturity Period  : ','Interest Rate : ','Purchase Year : '} ;
            str = strcat('Bond ',48+i-T);
            dlgtitle = str ;
            dims = [1 35];
            definput = {'1','0','1'} ;
            answer = inputdlg(prompt,dlgtitle,dims,definput) ;
            answer = str2double(answer) ;
            Maturity(i) = answer(1) ;
            Interest_Rates(i) = answer(2) ;
            Purchase_Years(i) = answer(3) ;
            % If any bond exceeds given time period.
            while Purchase_Years(i)+Maturity(i)>T+1
                str = "This Bond Exceeds Given Time Period." ;
                question_ans = questdlg(str,'Error','Ok','Ok') ;
                prompt = {'Maturity Period  : ','Interest Rate : ','Purchase Year : '} ;
                str = strcat('Bond',48+i-T) ;
                dlgtitle = str ;
                dims = [1 35] ;
                definput = {'1','0','1'} ;
                answer = inputdlg(prompt,dlgtitle,dims,definput) ;
                answer = str2double(answer) ;
                Maturity(i) = answer(1) ;
                Interest_Rates(i) = answer(2) ;
                Purchase_Years(i) = answer(3) ;
            end
        end
    end
else
    % For random generation
    % Generate random maturity durations
    Maturity = randi([1 T-1],total_bonds,1) ; 
    % Bond 1 has a maturity period of 1 year
    Maturity(1:T) = 1 ; 
    % Generate random yearly interest rate for each bond
    Interest_Rates = randi([4 8],total_bonds,1) ; 
    % Bond 1 has an interest rate of 0 (not invested)
    Interest_Rates(1:T) = 0 ; 
    % Generate random purchase years for each option
    Purchase_Years = zeros(total_bonds,1) ;
    % Bond 1 is available for purchase every year
    Purchase_Years(1:T)=1:T ;  
    for i=1:N
        % Generate a random year for the bond to mature 
        % before the end of the T year period
        Purchase_Years(i+T) = randi([1 T-Maturity(i+T)+1]) ;
    end
end

% Return after one year of interest
rt = 1 + Interest_Rates/100 ; 
% Compute the return at the end of the 
% maturity period for each bond
r = rt.^Maturity ; 

% Compute the years where each bond 
% reaches maturity at the end of the year
Maturity_Years = Purchase_Years + Maturity - 1 ;

%%
% Plot the Bond Detils like : Bond Number, Interest Rate
%                             Maturity, Purchase Year
plot_Investments(N,Purchase_Years,Maturity,Interest_Rates) ;

%%
% Create a matrices A,b and c for solving this linear problem
% using Simplex Method.
% System is : Ax=b which is constraints of this LPP.
% And matrix c is generated function for maximize this LPP.
A = zeros(T,total_bonds) ;
b = zeros(T,1) ;
b(1,1) = Initial ;
c = zeros(1,total_bonds) ;
for i=1:T
    for j=1:total_bonds
       A(i,j) = (Purchase_Years(j,1) == i) ; 
       A(i,j) = A(i,j)-((Maturity_Years(j,1) == i-1)*r(j,1)) ;
    end
end
for i=1:total_bonds
   c(1,i) = -(Maturity_Years(i,1) == T)*r(i,1) ;
end

%%
% Phase one
[m,n]                = size(A);
basicvar             = n+1:n+m;
phase                = zeros(m+1,n+m+1);
phase(1:m,1:n)       = A;
phase(end,n+1:end-1) = 1;
phase(1:m,end)       = b(:);
phase(1:m,n+1:n+m)   = eye(m);
 
for i = 1:m %now make all entries in bottom row zero
    phase(end,:) = phase(end,:)-phase(i,:);
end

[phase,basicvar] = simplex(phase(1:m,1:n+m),...
    phase(1:m,end),phase(end,1:n+m),basicvar);

%%
% Phase two
A = phase(1:m,1:n);
b = phase(1:m,end);
 
[phase,basicvar]   = simplex(A,b,c,basicvar);
[nRow,nCol] = size(phase);
b = phase(1:nRow-1,nCol);

%%
% Maximized Return for given Initial Amount
sol = zeros(total_bonds,1) ;
for i=1:length(basicvar)
    if b(i)>0
       sol(basicvar(i),1) = b(i) ;  
    end
end

ans1=sol.*(-1.*(c')) ; 
ans2=sum(ans1) ;
message = sprintf...
    ('After %d years, the return for the initial $%g is $%g',T,Initial,ans2) ;
uiwait(msgbox(message,'Return'));

%%
% Plot the bond details that which bonds can be buy
plot_Investments(N,Purchase_Years,Maturity,Interest_Rates,sol)
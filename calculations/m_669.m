clearvars, clc

% wuhan spike

n = 669; % length of genome


% number of genomes after k mutations

G = zeros(4,1); 
for k=1:4
    G(k) = nchoosek(n,k)*3^k;
end

disp('QA1');

N1 = G(1) % 1st step

%{
N1 =

        2007

%}

disp('QA2');

b = 5/223; % non-neutral sites from bloom escape calculator
a = 1-b;

ntrl = a*N1
not_ntrl = b*N1

%{
ntrl =

    45


not_ntrl =

        1962
%}

disp('QA3');
N = nchoosek(n,0) + G(1) + G(2) + G(3) % 0,1,2,3-step

%{
N =

   1.3434e+09
%}

disp('QA4');

N4 = G(4) % 4-step

%{
 N4 =

   6.7000e+11
%}

disp('QA5');

%{
 

ntrl =

   6.5498e+11


not_ntrl =

   1.5022e+10
%}

ntrl = N4*a
not_ntrl = N4*b

disp('QA6')

N = n; % initial genome length

for k=1:4
    N = N + G(k)
    sprintf('k=%d, N=%d, ntrl=%.2f, not_ntrl=%.2f',k,N,N*a,N*b)
    
end

%{
N =

        2676


ans =

    'k=1, N=2676, ntrl=2616.00, not_ntrl=60.00'


N =

     2013690


ans =

    'k=2, N=2013690, ntrl=1968540.00, not_ntrl=45150.00'


N =

   1.3434e+09


ans =

    'k=3, N=1343360028, ntrl=1313239848.00, not_ntrl=30120180.00'


N =

   6.7135e+11


ans =

    'k=4, N=671345855859, ntrl=656293258194.00, not_ntrl=15052597665.00'
%} 
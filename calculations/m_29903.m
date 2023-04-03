clearvars, clc

% wuhan  (full)

n = 29903; % length of genome


% number of genomes after k mutations

G = zeros(4,1); 
for k=1:4
    G(k) = nchoosek(n,k)*3^k;
end

disp('QA1');

N1 = G(1) % 1st step

%{
N1 =

       89709

%}

disp('QA2');

a = 20/64; 
b = 1-a;

ntrl = a*N1
not_ntrl = b*N1

%{
ntrl =

   2.8034e+04


not_ntrl =

   6.1675e+04
%}

disp('QA3');
N = nchoosek(n,0) + G(1) + G(2) + G(3) % 0,1,2,3-step

%{
                    1
                89709
           4023717777
      120313185250077

N =

   1.2032e+14

%}

disp('QA4');

N4 = G(4) % 4th step

%{

 N4 =

   2.6980e+18
%}

disp('QA5');

ntrl = N4*a
not_ntrl = N4*b

%{
 ntrl =

   8.4313e+17


not_ntrl =

   1.8549e+18
%}


disp('QA6')

N = n; % initial genome length

for k=1:4
    N = N + G(k)
    sprintf('k=%d, N=%d, ntrl=%.2f, not_ntrl=%.2f',k,N,N*a,N*b)
    
end


%{


    'k=1, N=119612, ntrl=37378.75, not_ntrl=82233.25'


    'k=2, N=4023837389, ntrl=1257449184.06, not_ntrl=2766388204.94'


    'k=3, N=120317209087466, ntrl=37599127839833.12, not_ntrl=82718081247632.88'


    'k=4, N=2698143496442064384, ntrl=843169842638145152.00, not_ntrl=1854973653803919360.00'

%}

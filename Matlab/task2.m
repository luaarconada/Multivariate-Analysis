%Load the clean data (cleaned in R in task1)
data_aux=readtable('clean_data.csv');

%Eliminate the first row, identity variable
data_aux = data_aux(:,2:26);

%Set the qualitative variables as categorical
data_aux.rbc=categorical(data_aux.rbc);
data_aux.pc=categorical(data_aux.pc);
data_aux.pcc=categorical(data_aux.pcc);
data_aux.ba=categorical(data_aux.ba);
data_aux.htn=categorical(data_aux.htn);
data_aux.dm=categorical(data_aux.dm);
data_aux.cad=categorical(data_aux.cad);
data_aux.appet=categorical(data_aux.appet);
data_aux.pe=categorical(data_aux.pe);
data_aux.ane=categorical(data_aux.ane);
data_aux.class=categorical(data_aux.class);
data_aux.bp=categorical(data_aux.bp);
data_aux.sg=categorical(data_aux.sg);
data_aux.al=categorical(data_aux.al);
data_aux.su=categorical(data_aux.su);

%Number of continuous variables
p1=10;
%Number of binary variables
p2=11;
%Number of categorical variables
p3=4;

%Put data as arrays creating in between different arrays acording to the
%type of variable
X1=table2array(data_aux(:,1:p1));
X2=table2array(data_aux(:,p1+1:p1+p2));
X3=table2array(data_aux(:,p1+p2+1:end));
X=[X1,double(X2),double(X3)];

%Important dimensions
[n,p]=size(X);

%Take care of the binary data
% Map categories to binary representation
binaryData = zeros(size(X2));
% Map "normal" to 0, "abnormal" to 1 for variables 1 and 2
binaryData(:, [1, 2]) = double(X2(:, [1, 2]) == "abnormal");
% Map "notpresent" to 0, "present" to 1 for variables 3 and 4
binaryData(:, [3, 4]) = double(X2(:, [3, 4]) == "present");
% Map "yes" to 0, "no" to 1 for variables 5, 6, 7, 9, and 10
binaryData(:, [5, 6, 7, 9, 10]) = double(X2(:, [5, 6, 7, 9, 10]) == "no");
% Map "poor" to 0, "good" to 1 for variable 8
binaryData(:, 8) = double(X2(:, 8) == "good");
% Map "notckd" to 0, "ckd" to 1 for variable 11
binaryData(:, 11) = double(X2(:, 11) == "ckd")

%%%%%%%%%%%%%%%%%%%%%%%% Distances %%%%%%%%%%%%%%%%%%%%%%%%
%Gower's distance
Sgower=gower2(X,p1,p2,p3);
Dgower=sqrt(2*(ones(n)-Sgower));
D2gower=2*(ones(n)-Sgower);


%%%%%%%% Joint metrics %%%%%%%%
%Mahalanobis' distance matrix for NUMERICAL variables
Dmaha=pdist(X1, 'mahal');
D2maha=squareform(pdist(X1, 'mahal').^2)

%Jaccard similarity and distance matrices for BINARY variables
[n,p]=size(binaryData);
Sjaccard=jaccard(binaryData);
D2jaccard=2*(ones(n)-Sjaccard);
Djaccard=sqrt(D2jaccard);

%Sokal-Michener similarity and distance matrices for binary variables
Ssokal=sokal(binaryData);
D2sokal=2*(ones(n)-Ssokal);
Dsokal=sqrt(D2sokal)
%I think we'll use Jaccard

%Compute similarity coefficients SC1 and SC4 for CATEGORICA variables
[n,p]=size(X3);
SC1=coincidencias(X3);
D2SC1=2*(ones(n)-SC1);
DSC1=sqrt(D2SC1);

SC4=coincidencias4(X3);
D2SC4=2*(ones(n)-SC4);
DSC4=sqrt(D2SC4);
%I think we'll use DSC1

%We won't, fuck

%Compute joint metric using relms2 and Mahalanobis' distance for
%quantitative variables and Hamming distance for qualitative variables
%We already have the Mahalanobis' distance matrix D_maha
%We compute the Hamming' distance matrix for categorical AND binary
%variables
Shamming=coincidencias(X(:,p1+1:end));
D2hamming=2*(ones(n)-Shamming);
%First term of joint metric (Pythagorean sum)
vg1=sum(sum(D2maha))/(2*n^2);
vg2=sum(sum(D2hamming))/(2*n^2);
D2=D2maha/vg1+D2hamming/vg2;
%Joint metric by RelMS
D2joint=relms2(D2maha,D2hamming);
%% 
%%%%%%%%%%%%%%%%%%%%%%%% MDS %%%%%%%%%%%%%%%%%%%%%%%%
%[Y,vaps,percent,cum]=coorp(D2gower)
%size(Y)

%corr_table=correlaciones2(X,Y(:,1:3),p1,3);

%influence(X,p1,p2);

%identif_cuantis(X(:,1:p1),Y);
%Xcuali=[double(X2),double(X3)];
%identif_cualis(Xcuali,Y);

%sensitividad_Gower(X,p1,p2)

%RelMDS mapping (Related Metric Scaling);
%[Y_joint,vaps_joint,percent_joint,acum_joint]=coorp(D2joint)

%corr_table_joint=correlaciones2(X,Y_joint(:,1:3),p1,3);
%% 
%%%%%%%%%%%%%%%%%%%%%%%% Clustering %%%%%%%%%%%%%%%%%%%%%%%%


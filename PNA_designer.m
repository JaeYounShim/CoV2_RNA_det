%2023.01, PNA designer
%Tmp:Predicted melting temp in DNA/PNA 
%TmnnDNA :Melting temp from calculated from DNA/DNA duplex (by SantaLucia, 1996)
%l: length[bp]
%fpyr=fractional pyrimidine contents(U,T,C), cf) purine contents(A,G)
%co=20.79;c1=0.83;c2=-26.13;c3=0.44;
%These calculation is based on the paper by Freier, Susan M. in 1986 (RNA) and Ursula Giesen in 1998 (DNA).
%This code can be adpated for RNA sequences.
% Tms=co+c1*TmnnDNA+C2*fpyr+c3*l;
%l=15; % the length of probes [mer], typically PNA probe is designed no longer than 40 bp (15-20 mer is usual) 


%% TmnnDNA [K]
%{ TmnnDNA rule 
% 1. Tm between 20-65C, 
% 2. no more than three consecutive G sequences 
% 3. terminal GC base pairs to minimize helix fraying 
% (which could invalidate the two-state approximation of double helix and random coil)
% 4. oligonucleotide is dissolved in 1M NaCl 
% 5. TmnnDNA=H/s+Rln(Ct/4), If self-complementary molecules ct/4 -> ct. R is gas constatnt [1.9787 cal/(Kmol)] 
%}
%initiation effects are ignored.
%TmnnDNA=H/(S+R*ln(ct/4); H is enthalpy change and S is entropy change.
clear
H(:,1)=[-8.4,-8.4,-6.5,-6.3,-7.4,-7.4,-8.6,-8.6,-6.1,-6.1,-7.7,-7.7,-10.1,-11.1,-6.7,-6.7]'; 
% S(:,1)=[-23.6,-23.6,-18.8,-18.5,-19.3,-19.3,-23.0,-23.0,-16.1,-16.1,-20.3,-20.3,-25.5,-28.4,-15.6,-15.6]';

S(:,1)=[-23.6,-23.6,-18.8,-18.5,-19.3,-19.3,-23.0,-23.0,-16.1,-16.1,-20.3,-20.3,-25.5,-28.4,-15.6,-15.6]';
Seq(:,1)=["AA","TT","AT","TA","CA","TG","GT","AC","CT","AG","GA","TC","CG","GC","GG","CC"]'; %from 5' to 3'
R=1.987; % [cal/(K*mol)] 

%%  probe designer 
probe_result={};probe_TmnnDNA=[];
Temp=45;range=5;%degree of celcius
gapbps=0; %size of gap between probes
z=1; %initial position 
l=1; %initial length 
zz=1; %save column sequence


%contstraints 
   %1. No longer than 40bps 
   %2. escape high purine content over 50%
   %3. No Purine stretch over 5bp (A&G)
   %4. No consecutive G residues. 
   %5. ideal legnth is between 12 and 18 mer. 


% RNA -> DNA 
clear Total_seq
Total_seq=splitGraphemes("AUGUCUGAUAAUGGACCCCAAAAUCAGCGAAAUGCACCCCGCAUUACGUUUGGUGGACCCUCAGAUUCAACUGGCAGUAACCAGAAUGGAGAACGCAGUGGGGCGCGAUCAAAACAACGUCGGCCCCAAGGUUUACCCAAUAAUACUGCGUCUUGGUUCACCGCUCUCACUCAACAUGGCAAGGAAGACCUUAAAUUCCCUCGAGGACAAGGCGUUCCAAUUAACACCAAUAGCAGUCCAGAUGACCAAAUUGGCUACUACCGAAGAGCUACCAGACGAAUUCGUGGUGGUGACGGUAAAAUGAAAGAUCUCAGUCCAAGAUGGUAUUUCUACUACCUAGGAACUGGGCCAGAAGCUGGACUUCCCUAUGGUGCUAACAAAGACGGCAUCAUAUGGGUUGCAACUGAGGGAGCCUUGAAUACACCAAAAGAUCACAUUGGCACCCGCAAUCCUGCUAACAAUGCUGCAAUCGUGCUACAACUUCCUCAAGGAACAACAUUGCCAAAAGGCUUCUACGCAGAAGGGAGCAGAGGCGGCAGUCAAGCCUCUUCUCGUUCCUCAUCACGUAGUCGCAACAGUUCAAGAAAUUCAACUCCAGGCAGCAGUAGGGGAACUUCUCCUGCUAGAAUGGCUGGCAAUGGCGGUGAUGCUGCUCUUGCUUUGCUGCUGCUUGACAGAUUGAACCAGCUUGAGAGCAAAAUGUCUGGUAAAGGCCAACAACAACAAGGCCAAACUGUCACUAAGAAAUCUGCUGCUGAGGCUUCUAAGAAGCCUCGGCAAAAACGUACUGCCACUAAAGCAUACAAUGUAACACAAGCUUUCGGCAGACGUGGUCCAGAACAAACCCAAGGAAAUUUUGGGGACCAGGAACUAAUCAGACAAGGAACUGAUUACAAACAUUGGCCGCAAAUUGCACAAUUUGCCCCCAGCGCUUCAGCGUUCUUCGGAAUGUCGCGCAUUGGCAUGGAAGUCACACCUUCGGGAACGUGGUUGACCUACACAGGUGCCAUCAAAUUGGAUGACAAAGAUCCAAAUUUCAAAGAUCAAGUCAUUUUGCUGAAUAAGCAUAUUGACGCAUACAAAACAUUCCCACCAACAGAGCCUAAAAAGGACAAAAAGAAGAAGGCUGAUGAAACUCAAGCCUUACCGCAGAGACAGAAGAAACAGCAAACUGUGACUCUUCUUCCUGCUGCAGAUUUGGAUGAUUUCUCCAAACAAUUGCAACAAUCCAUGAGCAGUGCUGACUCAACUCAGGCCUAA"); %total sequences of target region which we want to extract sequences from.  
    for i=1:size(Total_seq,1)
        if Total_seq(i,1)=="U"
            Total_seq(i,1)="T";
        end
    end


while z+l<size(Total_seq,1)
clear Targ Ct
Ct=4*10^(-6); % mole/L
Targ=Total_seq(z:z+l,1);%minimum probe length is 2bps 

% conver to complementary sequences 
clear Targ_com 
for i=1:1:size(Targ,1)
    if Targ(size(Targ,1)-i+1,1)=="G"
        Targ_com(i,1)="C";
    elseif Targ(size(Targ,1)-i+1,1)=="C"
        Targ_com(i,1)="G";
    elseif Targ(size(Targ,1)-i+1,1)=="A"
        Targ_com(i,1)="T";
    elseif Targ(size(Targ,1)-i+1,1)=="T"
        Targ_com(i,1)="A"; 
    end
end

% Calculation of Enthalpy and Entropy       
%based on seequences 
clear TmnnDNA_S TmnnDNA_Seq TmnnDNA_H TmnnDNA
for i=1:1:size(Targ_com,1)-1
    for ii=1:size(Seq,1)
        if strcat(Targ_com(i,1),Targ_com(i+1,1))==Seq(ii,1)
            TmnnDNA_Seq(i,1)=strcat(Targ_com(i,1),Targ_com(i+1,1));
            TmnnDNA_S(i,1)=S(ii,1);
            TmnnDNA_H(i,1)=H(ii,1);
        end   
    end
end

%based on initiation  
if any(Targ_com(:,1)=="G") || any(Targ_com(:,1)=="C")
    TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
    TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=-5.9;
elseif ~any(Targ_com(:,1)=="G") && ~any(Targ_com(:,1)=="C")
    TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
    TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=-9.0;
% elseif
%     TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
%     TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=0;
end

%based on symmetric correction: self-complementary sequences         
if all(Targ_com==Targ)
   TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
   TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=-1.4;
%    Ct=Ct*4;
end

%based on 5'-T
if Targ_com(1,1)=="T" 
    TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0.4;
    TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=0;
end
%based on A-3'
if Targ_com(size(Targ_com,1),1)=="A"
    TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0.4;
    TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=0;
end


% TmnnDNA extraction
TmnnDNA=sum(TmnnDNA_H(:,1))*1000/(sum(TmnnDNA_S(:,1))+R*log(Ct)/log(exp(1)))-273.15; % initiation component missing

%Pyrimidine percentage (Pyrimidine: T&C)
fpry=0;
for i=1:size(Targ_com,1)
    if Targ_com(i,1)=="T"||Targ_com(i,1)=="C"
        fpry=fpry+1;
    end
end
fpry_per=fpry/size(Targ_com,1);

%Tmns calculation
co=20.79;c1=0.83;c2=-26.13;c3=0.44;
Tms=co+c1*TmnnDNA+c2*fpry_per+c3*size(Targ_com,1);

% Tms between Temperature range or not 
if Temp-range<Tms && Tms<Temp+range && 1-fpry_per<=0.5
    probe_result{zz,1}=z;probe_result{zz,2}=z+l;probe_result{zz,3}=l+1;probe_result{zz,4}=strjoin(Targ_com(1:size(Targ_com,1),1),"");
    probe_result{zz,5}=Tms;
    probe_TmnnDNA(zz,1)=sum(TmnnDNA_H);probe_TmnnDNA(zz,2)=sum(TmnnDNA_S);probe_TmnnDNA(zz,3)=TmnnDNA;
    z
    zz=zz+1;
    z=z+l+1+gapbps;    
else
   l=l+1;
   if l>20
       z=z+1;
       l=1;
   end
end

end

%% constraints 
%컴포넌트 중에 GC비율이 5개 이상이다. 
%self complementary 넣기 (7bp 이상)
%no more consecutive 3G residues
%terminate with G or C. 
%oligonucleotides were dissolved in 1M NaCl

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


%% Random seq 
clear test 
test_temp=[];
for i=1:100
clear operators
operators = ["A","T","C","G"];
% Get 14 random indices from that list.
indexes = randi(length(operators), 1, 20);
test(i,1) = strjoin((operators(indexes)),"");
end

%% Extra) Input direct PNA seq.: measurement of melting temperature (Tms) 
% Do not Run repeatedly.
for k=1:size(test,1)
clear Targ
Targ=splitGraphemes(test(k,1)); %sequences of fragement to bring out probes, RNA sequence base
% Targ=splitGraphemes("ACATATG"); %sequences of fragement to bring out probes, RNA sequence base
Ct=4*10^(-6); % oligonucleotide concentration [mol/L]

% RNA -> DNA     
    for i=1:size(Targ,1)
        if Targ(i,1)=="U"
            Targ(i,1)="T";
        end
    end


% Calculation of Enthalpy and Entropy 

    %based on seequences 
    clear TmnnDNA_S TmnnDNA_Seq TmnnDNA_H TmnnDNA
    for i=1:1:size(Targ,1)-1
        for ii=1:size(Seq,1)
            if strcat(Targ(i,1),Targ(i+1,1))==Seq(ii,1)
                TmnnDNA_Seq(i,1)=strcat(Targ(i,1),Targ(i+1,1));
                TmnnDNA_S(i,1)=S(ii,1);
                TmnnDNA_H(i,1)=H(ii,1);
            end   
        end
    end  
    %based on initiation  
    if any(Targ(:,1)=="G") || any(Targ(:,1)=="C")
        TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
        TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=-5.9;
    elseif ~any(Targ(:,1)=="G") && ~any(Targ(:,1)=="C")
        TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
        TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=-9.0;
%     elseif
%         TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
%         TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=0;
    end

    %based on symmetric correction: self-complementary sequences 
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

       if all(Targ_com==Targ)
           TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0;
           TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=-1.4;
           Ct=Ct;
%        else
%            Ct=Ct/4;
       end
    
   %based on 5'-T, A-3'
       if Targ(1,1)=="T" 
            TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0.4;
            TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=0;
       end
   %based on A-3'  
       if  Targ(size(Targ,1),1)=="A"
            TmnnDNA_H(size(TmnnDNA_H,1)+1,1)=0.4;
            TmnnDNA_S(size(TmnnDNA_S,1)+1,1)=0;
       end



% TmnnDNA extraction
TmnnDNA=sum(TmnnDNA_H(:,1))*1000/(sum(TmnnDNA_S(:,1))+R*log(Ct)/log(exp(1)))-273.15; % initiation component missing

%Pyrimidine percentage (Pyrimidine: T&C)
fpry=0;
for i=1:size(Targ,1)
    if Targ(i,1)=="T"||Targ(i,1)=="C"
        fpry=fpry+1;
    end
end
fpry_per=fpry/size(Targ,1);

%Tmns calculation
co=20.79;c1=0.83;c2=-26.13;c3=0.44;
% co=-20.7167;c1=0.53627;c2=-45.133;c3=9.34722;

Tms=co+(c1*TmnnDNA)+(c2*fpry_per)+(c3*size(Targ,1));
[sum(sum(TmnnDNA_H)),sum(sum(TmnnDNA_S)), TmnnDNA, Tms];

test_temp(k,1)=Tms;
end
mean(test_temp)
std(test_temp)
function p = getPval(table1,table2,table12,table22,x1,x2)

%Full-blown LOGRANK procedure
%Merge the first columns of Table1 and Table2 (time intervals)
%and pick-up unique values
A=unique([table1(:,1);table2(:,1)]);
table=zeros(length(A),9); %matrix preallocation
%Out in the first column the time intervals
table(:,1)=A; 
%Put in the columns 2 and 3 and in the proper rows the deaths and alive
%taken from table1 columns 2 and 3
[~, ia ib]=intersect(table1(:,1),A);
table(ib,2:3)=table1(ia,2:3);
%Put in the columns 4 and 5 and in the proper rows the deaths and alive
%taken from table2 columns 2 and 3
[~, ia ib]=intersect(table2(:,1),A);
table(ib,4:5)=table2(ia,2:3);
%remove the rows where there arent't deaths in both treatments
table((table(:,2)==0 & table(:,4)==0),:)=[];
clear A c ia ib table1 table2
%fill the "pigeon-holes"
c=find(table(:,3)==0); %find the "pigeon-holes" of treatment 1
for I=1:length(c)
    if c(I)~=1
        %find the first interval time before the hole where there is almost 1
        %death
        J=find(table(1:c(I)-1,3)>0,1,'last');
        table(c(I),3)=table(J,3)-table(J,2);
        if ~isempty(table12)
        %find eventually censored data
            K=find((table12(:,1)<table(c(I),1) & table12(:,1)>=table(J,1)),1,'last');
            %Put in the hole how many subject were alive before the interval time
            %of the hole
            if ~isempty(K)
                table(c(I),3)=table(c(I),3)-sum(table12(K,2));
            end
        end
    else
        table(1,3)=length(x1);
    end
end
%Do the same for tratment 2
c=find(table(:,5)==0);
for I=1:length(c)
    if c(I)~=1
        J=find(table(1:c(I)-1,5)>0,1,'last');
        table(c(I),5)=table(J,5)-table(J,4);
        if ~isempty(table22)
            K=find((table22(:,1)<table(c(I),1) & table22(:,1)>=table(J,1)),1,'last');
            if ~isempty(K)
                table(c(I),5)=table(c(I),5)-sum(table22(K,2));
            end
        end
    else
        table(1,5)=length(x2);
    end
end
clear c I J K table12 table22

%Fill the table and compute the statistic variable
%Compute the total deaths and alive before the i-th time interval
table(:,6:7)=[sum(table(:,[2 4]),2) sum(table(:,[3 5]),2)];
%Compute the difference between observed deaths for treatment 1 and
%expected deaths in the hyphthesis that the treatments are similar
table(:,8)=table(:,2)-table(:,3).*table(:,6)./table(:,7);
%Log-rank statistic is the sum of column 8 values
J=sum(table(:,8)); UL=abs(J);
%Compute the contribute to the standard error
table(:,9)=prod(table(:,[3 5 6]),2).*(table(:,7)-table(:,6)) ...
    ./(table(:,7).^2.*(table(:,7)-ones(size(table,1),1)));
%find if there is some NaN (i.e. 0/0)
loc=isnan(table(:,9));
if any(loc)
    table(loc,9)=0;
end
V=sum(table(:,9)); SUL=sqrt(V); %Compute the totale standard error
K=J/V; HR=exp(K); HRci=[exp(K-1.96/SUL) exp(K+1.96/SUL)];
z=abs((UL-0.5)/SUL); %normalized UL with Yates'es correction
p=2*(1-0.5*erfc(-z/realsqrt(2))); %p-value

end
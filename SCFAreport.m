function SCFAreport(AS)
% [Fct,Fcor]=SCFAreport(AS);
Fct=AS.Fct;
o=ordreFct(Fct);
nf=numel(o);
F=Fct(:,o);
fprintf('\nFct:')
for k=1:size(F,1)
    fprintf('\n%d',k)
    for j=1:nf
        fprintf('\t%.3f',F(k,j));
    end
end
F=AS.CorFct(o,o);
printCor(F,'fCorr');
[r,rc]=triU(AS.R);
ar=atanh(r);
m=triU(AS.reprodR);
am=atanh(m);
z=abs(ar-am).*sqrt(AS.N-3);
[z,oo]=sort(z,1,"descend");
rc=rc(oo,:);
r=r(oo);
m=m(oo);
v=numel(AS.pertinent);
zcrit=sqrt(chi2inv(.95.^(1./[v,v*(v-1)/2]),1));
fprintf('\nr\tc\tRobs\tRmodl\t|z_diff|      z_crit(/var,gobal) = (%.4f, %.4f)\n',zcrit)
for j=1:max(5,sum(z>2))
    jr=rc(j,1);
    jc=rc(j,2);
    fprintf('%d\t%d\t%+.3f\t%+.3f\t%.4f\n',jr,jc,r(j),m(j),z(j));
end
% F=AS.R-AS.reprodR;
% sig=printCor(F,'CorrResid',AS.N);
dendr(AS,1);
end



function sig=printCor(F,titre,N)
if nargin>1
 fprintf('\n%s:\n',titre)
end
nf=size(F,2);
if nargin>2
    rcrit=sqrt(chi2inv(.95^(1/(nf-1)),1)/(N-1));
else
    rcrit=1;
end
for k=1:nf
    fprintf('\t%d',k);
end
sig=[];
for j=1:nf
    fprintf('\n%d',j);
    for k=1:nf
        fprintf('\t%.3f',F(j,k));
        if j~=k && abs(F(j,k))>rcrit
            fprintf('*');
            if k>j
                sig=[sig; [j k F(j,k)]];
            end
        end
    end
end
fprintf('\n');
for k=1:size(sig,1)
    fprintf('\n%2d %2d %.3f',sig(k,:))
end 
fprintf('\n');
end
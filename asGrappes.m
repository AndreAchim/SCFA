function AS=asGrappes(AS,alpha)
% AS=grappes(AS,alpha);
% AS.Z a été produit par linkage.m (linkage complet)
% alpha (défaut: .05), mais (1-alpha) adpté pour le nombre de paires dont
% le critère est le maximum
if nargin<2
    alpha=.05;
end
p=1-alpha;
nv=numel(AS.pertinent);
Gr{nv}=[];
for k=1:nv
    Gr{k}=AS.pertinent(k);
end
for k=1:size(AS.Z,1)
    a=Gr{AS.Z(k,1)};
    b=Gr{AS.Z(k,2)};
    nx=numel(a)*numel(b);
    crit=chi2inv(p.^(1/nx),nv-2);  % vaut 0 si nx==0
    if AS.Z(k,3)<crit
        Gr{end+1}=sort([a b]);
        Gr{AS.Z(k,1)}=[];
        Gr{AS.Z(k,2)}=[];
        AS.X2crit=crit;   % garder le dernier critère satisfait
    else
        Gr{end+1}=[];
    end
end
AS.reste=[];
nb=[];
m=numel(Gr);
for k=m:-1:1
    n=numel(Gr{k});
    if n<2
        if n==1
            AS.reste=[Gr{k} AS.reste];
        end
        Gr(k)=[];
    else
        nb=[n nb];
    end
end
[~,oo]=sort(nb,'descend');
AS.Gr=Gr(oo);
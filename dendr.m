function dendr(AS,fct)
% dendr(AS,fct);
% si fct est spécifié (e.g., 'sqrt'), cela devient l'échelle verticale
% si fct n'est pas de type char, il devient 'sqrt'
% au lieu des X2
% produit le dendrogramme des regroupements à partir de paires dans SCEFA
figure
% dendrogram(AS.Z,'ColorThreshold',chi2inv(.99,numel(AS.pertinent)));
Z=AS.Z;
crit=AS.X2crit; 
AxeY='X^2';
if nargin>1
    if ~ischar(fct), fct='sqrt'; end
    Z(:,3)=feval(fct,Z(:,3));
    crit=feval(fct,crit);
    AxeY=[fct ' (X^2)'];
end
dendrogram(Z);
lb=xticklabels;
v=str2num(lb);
nv=AS.pertinent(v);
if isfield(AS,'GrCoplan') && ~isempty(AS.GrCoplan)
    for j=1:numel(AS.GrCoplan)
        cop=AS.GrCoplan{j};
        for k=1:numel(cop)
            nv(nv==cop(k))=-nv(nv==cop(k));
        end
    end
end
noms{numel(AS.pertinent)}='';
for k=1:size(v)
    noms{k}=num2str(nv(k));
end
xticklabels(noms);
hold on
plot([0 numel(v)+1],crit*[1 1],'--k');
ylabel(AxeY);
titre=sprintf('N=%d',AS.N);
if isfield(AS,'titre')
    titre=[AS.titre '  ' titre];
end
title(titre);

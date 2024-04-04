function AS=asTuples(AS,k)
% AS=asTuplesAFD(AS,k);
% AS.GS(v,v) est la GrammSchmidt (en triangle supérieur) de R(v,v)
% Si k est un scalaire, optimise les variables de la rangée k de AS.Cpaire pour en minimiser le signal
% Autrement, optimise les valeurs de k(1:end-1) pour minimiser le signal dans k(end)
% Le champ AS.tmp retourne alors le critère, les crrélations avec les
% autres variables et les poids optimisés

% ajoute AS.Crit
% Po, si présent, contient les poids optimaux des numel(melange) variables
% corr(1,v) contient les corrélations de SPo avec chacune des autres variables
options=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolFun',1e-6,'TolX',1e-6);
if ~isfield(AS,'pertinent')
    AS.pertinent=1:AS.nv;
end
% cible=1:numel(AS.pertinent);
cible=AS.pertinent;
if numel(k)==1
    melange=AS.Cpaires(k,:);
else
    melange=k;
end
% si les deux plus grandes corrélations ne sont pas essentiellement égales,
% essayer avec nouvelles initialisations
crit=9e9;
mul=[1 .5 2];
for e=1:numel(mul)
    P=sign(AS.GS(:,melange(1:end-1))'*AS.GS(:,melange(end)));  % initialiser à +ou - 1 selon signe de corrélations
    P=mul(e)*P./sqrt(numel(P));
    [P,cr]=fminsearch(@(P) asCrit(P,AS.GS,melange,cible),P,options);
    [~,cor]=asCrit(P,AS.GS,melange,cible);
    cor=cor';
    c=sort(abs(cor));
    if cr<crit
        crit=cr;
        if numel(k)==1
            % AS.Crit(k)=cr;
            AS.Crit(k)=cor(:)'*cor(:); % bien qu'on ait minimisé la valeur absolue maximale des corrélations
            AS.Corr(:,k)=cor;
            AS.Ppaires(k,:)=P;
        else
        AS.tmp.Crit=cr;
        AS.tmp.Corr=cor;
        AS.tmp.Poids=P';
        end
    end
    if c(end)-c(end-1)<.0001
        break;
    end
end
if c(end)-c(end-1)>.0001
    keyboard  % l'optimisation n'a pas trouvé un vrai creux
end
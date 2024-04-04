function AS=asMultiSatur(AS,np)
% AS=asMultiSatur(AS,np);
% AS contient le champ .GS et le vecteur .reste
% Cette commande cherche à annuler la colonne AS.reste(1) de AS.GS
% à partir de np prédicteurs dont les rangs sont donnés dans AS.Var
% Si une combinaison de np prédicteurs annule suffisamment AS.reste(1)
% les saturations associées sont inscrites dans out.Fct et AS.reste(1) est
% enlevé
% autrement, AS.reste(1) est mis en négatif pour indiquer que pas résssi
% avec np prédicteurs et garder l'information pour essayer avec un np plus grand 
if isempty(AS.reste)
    return
end
nv=AS.nv;
C=nchoosek(AS.Var,np); % tous les croisements de Var
nc=size(C,1);
Crit=zeros(nc,1);
Corr=zeros(nv,nc);
P=zeros(nc,np);
G=AS.GS;
prob=zeros(nc,1); % des X2 avec dl différents, on comparera les probabilités
for j=1:nc
    melange=[C(j,:),AS.reste(1)];
    cible=setdiff(AS.pertinent,melange);
    AS=asTuples(AS,melange);
    Crit(j)=AS.tmp.Crit;
    P(j,:)=AS.tmp.Poids;
    Corr(:,j)=AS.tmp.Corr;
    prob(j)=1-chi2cdf(Crit(j)*(AS.N-1),numel(cible));
end
Crit=Crit*(AS.N-1);
z2=Corr.^2*(AS.N-1);
% ici, évaluer les annulations de signal
% et mettre en négatif les rangs de variables à retirer de reste
f=find(prob==max(prob),1);
if isempty(f) || prob(f)<.05 || sum(z2(:,f))>chi2inv(.95,numel(cible)) % solution d'annulation pas trouvée
% if isempty(f) || prob(f)<.05 || max(z2(:,f)>3.8415) % solution d'annulation pas trouvée
    AS.reste(1)=-AS.reste(1);
    return
end
% si les prédicteurs sont membres de groupes déjà formés, calculer tous les croisements
% et utiliser la moyenne des saturations des croisements de prédicteurs
gd=AS.GrDe(C(f,:));  % groupes des variables qui annulent le signal de reste(1)
if all(gd>0)
    crois=croise(AS.Gr{gd(1)}',AS.Gr{gd(2)}');
    for c=3:numel(gd)  % si plus que 2 prédicteurs
        crois=croise(crois,AS.Gr{gd(c)}');
    end
    nt=size(crois,1);
    CritTuple=zeros(nt,1);
    saturTuple=zeros(nt,np);
    for f=1:nt
        melange=[crois(f,:) AS.reste(1)];
        cible=setdiff(AS.pertinent,melange);
        % [CritTuple(f),Po,Co]=tuplesAFD(G,melange,setdiff(cible,melange));
        AS=asTuples(AS,melange);
        CritTuple(f)=AS.tmp.Crit;
        Po=AS.tmp.Poids;
        Co=AS.tmp.Corr;
        if abs(log10(Po))>10
            % f=[];
            break;
        end
        CritTuple(f)=Co'*Co*(AS.N-1);
        saturTuple(f,:)=Po.*sum(AS.Fct(melange(1:np),:),2)';
    end
    % if isempty(f)
    %    break;
    % end
    % saturations du meilleur critère ou la moyenne des saturations ou moy pondérée par 1/crit? 
    % AS.Fct(AS.reste(1),gd)=saturTuple(CritTuple==min(CritTuple),:);
    po=(1./CritTuple.^2)';   % mis au carré pour insister davantage sur les bas critères
    po=po/sum(po);
    AS.Fct(AS.reste(1),gd)=po*saturTuple;
    % AS.Fct(AS.reste(1),gd)=mean(saturTuple);
    AS.reste(1)=[];
else
    % keyboard
    AS.reste=[AS.reste(2:end) -AS.reste(1)];
    return
end

% keyboard
%     [cr,P,cor]=tuplesAFD([G IC],[v+f reste(k)],pertinent);
%     saturTuple(end,:)=Po'.*SaturIC(f)';
%     cor(cibleX2==0)=0;  % est-ce encore pertinent?
%     cor(reste)=0;
%     cor=cor.^2*(N-1);
%     cr=sum(cor);
%     X2=chi2inv(.95,numel(cibleX2));
%     if cr>X2      % si les indicatrices combinées n'annulent pas bien. utiliser
%         Fct(reste(k),f)=mean(saturTuple(1:end-1,:));  % la moyenne des indicatrices seules
%     else
%         Fct(reste(k),f)=saturTuple(end,:);
%     end
%     reste(k)=-reste(k);
% end
% reste(reste<0)=[];
end

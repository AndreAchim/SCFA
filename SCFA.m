function AS=SCFA(R,N)
% AS=SCFA(R,N);    ou AS=SCFA(dat);
% R peut être une matrice de covariances
% AS est une structure avec divers champs dont
% .R, .N, .GS, .Fct, .CorFct, .Gr, .Ppaires, .Cpaires
AS.dat=R; 
if nargin==1
    N=size(R,1);
    R=cov(R);
end
AS.et=sqrt(diag(R));  % ceci est destiné à pouvoir exprimer la solution factorielle en termes des variables d'origine
iet=1./AS.et;
AS.R=R.*(iet*iet');
AS.N=N;
AS.nv=numel(iet);
AS.GS=chol(R);
AS=asOrphelines(AS);
% Identifier les variables indicatrices, s'il y en a
AS=asPairesIndicatrices(AS);
% On pourrait avoir ici une vérification de variables oporphelines, avec
% un poids de signe opposé à corrélation, pour paire dont l'autre membre a une nette corrélation
% mais les signes opposés apparaissent peut-être uniquement pour paires de
% deux variables orphelines
% Regrouper les variables à partir des paires non-significatives
AS=asDistances(AS);
AS.Z=linkage(AS.Dist,'complete');
AS=asGrappes(AS);  % défini aussi AS.reste
AS=asCoplanaire(AS); % si trois groupes n'occupent qu'un plan, en enlever un
AS.ng=numel(AS.Gr);
AS.Fct=zeros(AS.nv,AS.ng);
% Estimer saturations des facteurs, puis leurs corrélations
AS=asSaturations(AS);   % prépare aussi AS.Var: une variable par groupe et AS.GrDe
AS=asCorrFct(AS);
% Annuler les variables restantes à partir des groupes identifiés
AS=asVariablesMulti(AS); 
AS.reprodR=AS.Fct*AS.CorFct*AS.Fct';
% étape 2 encore à implanter
%% 2) à partir de groupes, s'il y en a, et de variables restantes
%% en (2), s'il y a des groupes, on n'utilisera qu'une variable de chacun
if ~isempty(AS.reste)
    mesg='Les facteurs identifiés n''expliquent pas les variables: ';
    for j=AS.reste(:)'
        mesg=[mesg sprintf(' %d',j)];
    end
    warning(mesg);
end
%     % while ~isempty(AS.reste)
%     pertinent=AS.pertinent;
%     if ng>0
%         for k=1:ng
%             pertinent=setdiff(pertinent,AS.Gr{k}(2:end));
%         end
%     end
%     AS=multiSatur(AS,np,pertinent);
% end

% réestimer les saturations des facteurs doublets
% for k=1:ng
%     paire=find(AS.Fct(:,k)~=0);
%     % if max(abs(Fct(paire,k)))>.95
%         if numel(paire)==2
%             % [~,P]=tuplesAFD(G,paire,setdiff(1:v,paire));
%             [~,P]=tuplesAFD(AS.GS,paire);
%             [IC(:,k),Fct(paire,k),SaturIC(k)]=fusionne(AS.GS,paire,P,paire(:)');
%             % Fct(paire,k)=satur(G(:,paire),P);
%         endMM
%     % end
% end
% CorFct=IC'*IC;
% %% 2:   Si pas ou peu d'indicatrices unifactorielles
%

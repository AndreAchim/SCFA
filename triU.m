function [V,rc]=triU(R)
% [V,rc]=triU(R);
% Pour une matrice R(v,v), retourne V(v**(v-1)/2,1) des valeurs en haut de
% la diagonale
% rc(v**(v-1)/2,2) donne rangée,colonne pour chaque entrée
v=size(R,1);
V=zeros(v*(v-1)/2,1);
rc=zeros(v*(v-1)/2,2);
j=0;
for k=1:v-1
    V(j+1:j+v-k)=R(k,k+1:v);
    rc(j+1:j+v-k,1)=k;
    rc(j+1:j+v-k,2)=k+1:v;
    j=j+v-k;
end

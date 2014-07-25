% ******************* Adaptive Matrix Completion **********************
% Input:
%   X         :  (d x n) original matrix 
%   m         :  size of the sampling of column
% Output:
%   Xopt      :  recovery matrix
%
% **********************************************************************
% by Jamie Zemin Zhang
% 07/23/2014

clear all
close all
clc

addpath('test_figures/')                   ;
addpath('toolbox')                         ;

%% Problem set up
X             =  double(imread('coins.png')) ;    % original matrix
r             =  20                          ;

[U,S,V]       =  svd(X)                      ;
S(r+1:end,r+1:end) = 0                       ;
X             =  U*S*V'                      ;

[d,n]         =  size(X);
p             =  0.2                         ;
m             =  round(d*p)                  ;
U             =  X(:,1)/norm(X(:,1))         ;
X_result      =  zeros(size(X))              ;   
X_result(:,1) =  X(:,1)                      ;
epsilon       =  1e-1                        ;
Omega         =  randi(d,[1,m])              ;

flag          =  zeros(1,n)                  ;    % sample all or not
flag(1)       =  1                           ;

%% Process
fprintf('%3s\t%16s\t%7s\n',...
            'col','error','flag')       ;
        
for col = 2:n
    Uo      =  U(Omega,:)              ;
    Puo     =  Uo*inv(Uo'*Uo)*Uo'      ;
    xt      =  X(:,col)  ;
    xo      =  xt(Omega) ;
    
    val     =  norm(xo-Puo*xo)^2       ;
    if val > epsilon
        v   =  xt    ;
        for i = 1:size(U,2) ;
            v = v - gs_proj(U(:,i),v)  ;
        end
        U     =  [ U , v ]             ;
        Omega =  randi(d,[1,m])        ;
        X_result(:,col) = X(:,col)     ;
        flag(col)  =  1;
    else
        X_result(:,col) = U*inv(Uo'*Uo)*Uo'*xo ;
        flag(col)  =  0;
    end
    fprintf('%3d\t%10.4f\t%10.4f\n',...
                col,val,flag(col))  ;
            
end

figure;
subplot(121);imagesc(X);colormap(gray);title('Original');
subplot(122);imagesc(X_result);colormap(gray);title('Rec');

sampling_rate = (sum(flag)*d + (n-sum(flag))*m)/(d*n)
    
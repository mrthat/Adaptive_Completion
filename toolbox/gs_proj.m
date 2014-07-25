function proj = gs_proj( u , v )
% **** projection function for Gram-Schimit ****
% u and v are both vectors
% by Jamie Zemin Zhang
% 07/23/0214
% **********************************************

proj = (u'*v)/(u'*u)*u ;

end


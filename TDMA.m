%==========================================================================
%
% TDMA:  Solves the tridiagonal linear system Ax = b for x using the
% tridiagonal matrix algorithm (i.e. the Thomas algorithm).
%
% Copyright © 2022 Llorente-Lázaro, Víctor Javier
% Last Update: 2022-04-24
% Website: https://tamaskis.github.io
% Contact: victor.javier.llorente@gmail.com
%
%--------------------------------------------------------------------------
%
% ------------
% Description:
% ------------
%
%  [ bm(1) cm(1)                                              ] [   vecsol(1)    ]   [   dm(1)    ]
%  [ am(2) bm(2) cm(2)                                        ] [   vecsol(2)    ]   [   dm(2)    ]
%  [       am(3) bm(3) cm(3)                                  ] [                ]   [            ]
%  [             ...   ...   ...                              ] [      ...       ] = [    ...     ]
%  [                   ...   ...        ...                   ] [                ]   [            ]
%  [                         am(nequ-1) bm(nequ-1) cm(nequ-1) ] [ vecsol(nequ-1) ]   [ dm(nequ-1) ]
%  [                                    am(nequ)   bm(nequ)   ] [  vecsol(nequ)  ]   [  dm(nequ)  ]
%
% Conditions: 
% 1. Strictly Diagonally Dominant (SDD) matrix or Symmetric Positive Definite (SPD) matrix
% 2. bm(1) ~= 0
% The condition 1 is a sufficient condition but not a necessary condition
%
% ----------------
% MATLAB function: 
% ----------------
%   vecsol = TDMA( am, bm, cm, dm )
%
% ------
% INPUT:
% ------
%   Double, (nequ)x1 array   :: bm       - diagonal of A
%   Double, (nequ-1)x1 array :: am, cm   - lower and upper diagonal of A (Note: am(1) = cm(nequ) = 0)
%   Double, (nequ)x1 array   :: dm       - source b
%
% -------
% OUTPUT:
% -------
%   Double, (nequ)x1 array   :: vecsol   - solution x
%
% ----------------
% LOCAL VARIABLES:
% ----------------
%   Integer                  :: nequ     - number of equations of the system
%   Double                   :: parsmall - small parameter to avoid division by 0
%   Double, (nequ)x1 array   :: vecaux   - auxiliar vector
%
% -----------
% REFERENCES:
% -----------
%   1. https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
%   2. http://www.thevisualroom.com/tri_diagonal_matrix.html
%   3. https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
%
%==========================================================================
function vecsol = TDMA( am, bm, cm, dm )

    %% determines number of equations
    nequ = length( dm );
    
    %% Check dimensions
    if ( length(am) ~= nequ ) || ( length(bm) ~= nequ ) || ( length(cm) ~= nequ )
        warning('Dimension of arrays do not match')
        return
    end
    
    %% allocates
    vecaux = zeros( nequ, 1 );   
    vecsol = vecaux;
    
    %% small number
    parsmall = 1e-15;   
    
    %% extracts first element
    vecaux(1) = cm(1) / bm(1);
    vecsol(1) = dm(1) / bm(1);

    %% forward elimination
    for i = 2:nequ
        denom = bm(i) - am(i) * vecaux(i - 1);
        if abs( denom ) < parsmall
            denom = parsmall * sign( denom );
        end
        vecaux(i) = cm(i) / denom;
        vecsol(i) = ( dm(i) - am(i) * vecsol(i - 1) ) / denom;
    end
    
    %% backward substitution
    for i = nequ-1:-1:1
        vecsol(i) = vecsol(i) - vecaux(i) * vecsol(i + 1);
    end
    
end
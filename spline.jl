using LinearAlgebra

function spline(S,X)
#      DIMENSION X[N],XS[N],S[N]
#      PARAMETER (NMAX=5001)
#      DIMENSION A(NMAX),B(NMAX),C(NMAX)
#
#-------------------------------------------------------
#     Calculates spline coefficients for X(S).          |
#     Natural end conditions are used (zero 3rd         |
#      derivative over first, last intervals).          |
#                                                       |
#     To evaluate the spline at some value of S,        |
#     use SEVAL and/or DEVAL.                           |
#                                                       |
#     S        independent variable array (input)       |
#     X        dependent variable array   (input)       |
#     XS       dX/dS array                (calculated)  |
#     N        number of points           (input)       |
#                                                       |
#-------------------------------------------------------
#     

N = length(S)
B = zeros(N-1)
A = zeros(N)
C = zeros(N-1)
D = zeros(N)

for i=2:N-1
    DSM = S[i] - S[i-1]
    DSP = S[i+1] - S[i]
    B[i] = DSP
    A[i] = 2.0*(DSM+DSP)
    C[i] = DSM
    D[i] = 3.0*((X[i+1]-X[i])*DSM/DSP + (X[i]-X[i-1])*DSP/DSM)
end

#---- set zero 3rd derivative end conditions
  A[1] = 1.0
  C[1] = 1.0
  D[1] = 2.0*(X[2]-X[1]) / (S[2]-S[1])

  B[N-1] = 1.0
  A[N] = 1.0
  D[N] = 2.0*(X[N]-X[N-1]) / (S[N]-S[N-1])

  if(N==2)
#----- if only two points are present, specify zero 2nd derivative instead
#-     (straight line interpolation will result)
   B[N] = 1.0
   A[N] = 2.0
   D[N] = 3.0*(X[N]-X[N-1]) / (S[N]-S[N-1])
  end

#---- solve for derivative array XS
SYS = Tridiagonal(B,A,C)
XS = SYS\D    

  return XS
end # SPLINE


function SEVAL(SS,X,XS,S)
#--------------------------------------------------
#     Calculates X(SS)                             |
#     XS array must have been calculated by SPLINE |
#--------------------------------------------------
i_low = 1
i = N

while(i-i_low>1)
i_mid = Int((i+i_low)/2)	
  if(SS < S[i_mid])
   i = i_mid
  else
   i_low = i_mid
  end
end

DS = S[i] - S[i-1]
T = (SS - S[i-1]) / DS

CX1 = DS*XS[i-1] - X[i] + X[i-1]
CX2 = DS*XS[i]   - X[i] + X[i-1]

XX = T*X[i] + (1.0-T)*X[i-1] + (T-T*T)*((1.0-T)*CX1 - T*CX2)

return XX 
end # SEVAL

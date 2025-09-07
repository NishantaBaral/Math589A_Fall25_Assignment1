import math
import quartic_solver 


def two_ladders_width(L1, L2, h):
    '''Given two ladders of lengths L1 and L2 leaning against opposite walls
    of a narrow alley, and crossing at height h, return the width of the alley.'''

    a = (L1/h)**2
    b = (L2/h)**2
    
    #The dimensionless form of the mono-quadrtic equation: X = Xi
    #X^4+AX^3+BX^2+CX+D=0

    A = -2.0*(a+b-2.0)
    B = b**2 + 4*a*b - 6*b + a**2 - 6*a
    C = 2.0*(-a*b**2 + b**2 - a**2*b + a**2 + 4*a*b)
    D = a**2*b**2 - 2*a*b**2 + b**2 - 2*a**2*b - 2*a*b + a**2

    roots = quartic_solver.solve_quartic(1, A, B, C, D)


    tol = 1e-12

    width = []
    for i in roots:
        z = complex(i)
        if abs(z.imag) <= tol and z.real > 0.0:
            w = h * math.sqrt(z.real)     # Xi = (w/h)^2  -> w = h*sqrt(Xi)
            #comapre to see if w actuallt satisfies the original equation in case there're more than one real roots
            left_side = 1.0/(math.sqrt(L1**2 - w**2)) + 1.0/(math.sqrt(L2**2 - w**2))
            right_side = 1.0/h
            if abs(left_side - right_side) <= tol:
                width.append(w)
            
    if width == []:
        raise ValueError("dimensions aren't physically possible")

    return width
    
    


if __name__ == "__main__":
    print(two_ladders_width(40.0, 30.0, 10.0)) 
import math, cmath
from cubic_solver import solve_cubic as solve_cubic
from quadratic_solver import solve_quadratic as solve_quadratic


def cleaning_function(vals, tol=1e-12, ndigits=12):
    return tuple(
        round(z.real, ndigits) if abs(z.imag) < tol
        else complex(round(z.real, ndigits), round(z.imag, ndigits))
        for z in map(complex, vals)
    )

def solve_quartic(a,b,c,d,e):
    if a == 0:
        return solve_cubic(b,c,d,e)  # Degenerates to cubic
    
    if a==0 and b ==0:
        return solve_quadratic(c,d,e)  # Degenerates to quadratic
    if a==0 and b==0 and c==0:
        return False  # Degenerates to linear or constant equation


    tol = 1e-12

    #quartic equation ax^4+bx^3+cx^2+dx+e
    #mono quartic equation x^4+a_3x^3+a_2x^2+a_1x+a_0
    a0 = e/a
    a1 = d/a
    a2 = c/a
    a3 = b/a

    #coefficients of resolvent cubic b_3x^3+b_2x^2+b_1x+b_0

    b3 = 1
    b2 = -a2
    b1 = a1*a3 - 4*a0
    b0 = 4*a0*a2 - a1*a1 - a0*a3**2

    roots_of_resolvent_cubic = solve_cubic(b3,b2,b1,b0)

    real_root = [r for r in roots_of_resolvent_cubic if abs(r.imag) < tol]
    #I've picked the first real root but we can pick any real root. And we're guaranteed at least one
    #real root as a solution to the cubic. Hence, no error handling necessary.
    beta = real_root[0]

    sigma1, sigma2 = solve_quadratic(1.0, a3, (a2 - beta))   # s^2 + a3s + (a2-beta)=0
    eta1,   eta2   = solve_quadratic(1.0, -beta, a0)         # p^2 - beta*p + a0=0

    # Cleaning floating point errors
    sigma1, sigma2 = cleaning_function([sigma1,sigma2])
    eta1, eta2 = cleaning_function([eta1,eta2])


    # Residuals for the *actual x-coefficient*  target is σ1*η2 + σ2*η1  ≈ -a1
    res_keep  = abs((sigma1*eta2 + sigma2*eta1) + a1)
    # swap (eta2,eta1):  target is σ1*η1 + σ2*η2  ≈ -a1
    res_swap  = abs((sigma1*eta1 + sigma2*eta2) + a1)

    # Pick the assignment that better matches -a1
    if res_swap < res_keep:
        eta1, eta2 = eta2, eta1

    # Final two quadratics
    root1, root2 = solve_quadratic(1.0, -sigma1, eta1)
    root3, root4 = solve_quadratic(1.0, -sigma2, eta2)

    return cleaning_function([root1, root2, root3, root4])

def main():
    tests = [
        (1, 0, 0, 0, -1),  # roots of x^4 - 1 = 0 (2 real and 2 complex roots)
        (1, 0, 1, 0, -1),  # roots of x^4 + x^2 - 1 = 0 (2 real and 2 complex roots)
        (0, 1, -3, 2, 0),  # a=0 => cubic: x^3 - 3x^2 + 2x = 0  (roots 0,1,2)
    ]
    for a, b, c, d, e in tests:
        roots = solve_quartic(a, b, c, d, e)
        print(f"solve_quartic({a}, {b}, {c}, {d}, {e}) -> {roots}")

if __name__ == "__main__":
    main()

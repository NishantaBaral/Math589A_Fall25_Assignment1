import math, cmath
from quadratic_solver import solve_quadratic as solve_quadratic

def cleaning_function(vals, tol=1e-12, nd=12):
    xs = vals if isinstance(vals, (list, tuple)) else [vals]
    out = []
    for v in xs:
        z = complex(v)
        re = round(z.real, nd)
        im = 0.0 if abs(z.imag) < tol else round(z.imag, nd)
        out.append(complex(re, im))
    return tuple(out)

def half_power(z):
   return cmath.exp(0.5 * cmath.log(z))


def solve_cubic(a,b,c,d,tol=1e-12):

  if a == 0:
    return cleaning_function(solve_quadratic(b,c,d))


  #coefficient of the depressed cubic x^3+px+q
  p = (3*a*c-b*b)/(3*a*a)
  q = (2*b**3 - 9*a*b*c + 27*a*a*d)/(27*a**3)
  if abs(p.imag) <= tol:p = p.real
  if abs(q.imag) <= tol: q = q.real

  #handling edge case where p=0. depressed cubic is x^3+q=0
  if abs(p) == 0:

    r = abs(-q+0j)
    theta = cmath.phase(-q)

    alpha1 = r**(1/3) * (cmath.cos(theta/3) + 1j*math.sin(theta/3))
    alpha2 = r**(1/3) * (cmath.cos(theta/3 + 2*math.pi/3) + 1j*cmath.sin(theta/3 + 2*math.pi/3))
    alpha3 = r**(1/3) * (cmath.cos(theta/3 + 4*math.pi/3) + 1j*cmath.sin(theta/3 + 4*math.pi/3))

    root1 = alpha1 - b/(3*a)
    root2 = alpha2 - b/(3*a)
    root3 = alpha3 - b/(3*a)

    return cleaning_function([root1,root2,root3])
  
  #handling the case where the cofficients are complex numbers
  if any(abs(complex(v).imag) > tol for v in (p, q)):
     k = complex(half_power((-4*p)/3))
     z = complex(((3*q)/p)*(half_power(-3/(4*p))))

     alpha1 = 1/3*cmath.acos(z)
     alpha2 = 1/3*(cmath.acos(z) + 2*math.pi)
     alpha3 = 1/3*(cmath.acos(z) + 4*math.pi)

     root1 = k*cmath.cos(alpha1) - b/(3*a)
     root2 = k*cmath.cos(alpha2) - b/(3*a)
     root3 = k*cmath.cos(alpha3) - b/(3*a)

     return cleaning_function([root1,root2,root3])
     

  if p<0:
    #changing depressed cubic into 4cos^3(theta)-3cos(theta) = z; where x is scaled by kcos(theta)
    z = complex(((3*q)/p)*(half_power(-3/(4*p))))
    k = complex(half_power((-4*p)/3))


    if abs(z) < 1:
      alpha1 = 1/3*cmath.acos(z)
      alpha2 = 1/3*(cmath.acos(z) + 2*math.pi)
      alpha3 = 1/3*(cmath.acos(z) + 4*math.pi)

      root1 = k*cmath.cos(alpha1) - b/(3*a)
      root2 = k*cmath.cos(alpha2) - b/(3*a)
      root3 = k*cmath.cos(alpha3) - b/(3*a)

      return cleaning_function([root1,root2,root3])

    else:
      alpha1 = (1/3) * cmath.acosh(z+0j)
      alpha2 = (1/3) * (cmath.acosh(z+0j) + 2j*math.pi)
      alpha3 = (1/3) * (cmath.acosh(z+0j) + 4j*math.pi)

      root1 = k * cmath.cosh(alpha1) - b/(3*a)
      root2 = k * cmath.cosh(alpha2) - b/(3*a)
      root3 = k * cmath.cosh(alpha3) - b/(3*a)

      return cleaning_function([root1, root2, root3])   

  if p > 0:
    #changing depressed cubic into 4cos^3(theta)-3cos(theta) = z; where x is scaled by kcos(theta)
    z = complex(((-3*q)/p)*(half_power(3/(4*p))))
    k = complex(half_power((4*p)/3))
    if abs(z) < 1:
      alpha1 = 1/3*cmath.asin(z)
      alpha2 = 1/3*(cmath.asin(z) + 2*math.pi)
      alpha3 = 1/3*(cmath.asin(z) + 4*math.pi)

      root1 = k * cmath.sin(alpha1) - b/(3*a)
      root2 = k * cmath.sin(alpha2) - b/(3*a)
      root3 = k * cmath.sin(alpha3) - b/(3*a)

      return cleaning_function([root1, root2, root3])
      
    else:
      alpha1 = (1/3) * cmath.asinh(z+0j)
      alpha2 = (1/3) * (cmath.asinh(z+0j) + 2j*math.pi)
      alpha3 = (1/3) * (cmath.asinh(z+0j) + 4j*math.pi)

      root1 = k * cmath.sinh(alpha1) - b/(3*a)
      root2 = k * cmath.sinh(alpha2) - b/(3*a)
      root3 = k * cmath.sinh(alpha3) - b/(3*a)

      return cleaning_function([root1, root2, root3])

  
def main():
    tests = [
        (1, 0, 0, -1),     # roots of x^3 - 1 = 0 (1 and two complex cube roots)
        (1, -6, 11, -6),   # roots [1.0, 2.0, 3.0]
        (1,2,3,4),
        (1+1j,-6+10j,11,-6+2j),
    ]
    

    for a, b, c, d in tests:
        roots = solve_cubic(a, b, c, d)
        print(f"solve_cubic({a}, {b}, {c}, {d}) -> {roots}")

if __name__ == "__main__":
    main()

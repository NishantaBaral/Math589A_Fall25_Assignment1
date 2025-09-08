import math,cmath

def cleaning_function(vals, tol=1e-12, nd=12):
    xs = vals if isinstance(vals, (list, tuple)) else [vals]
    out = []
    for v in xs:
        z = complex(v)
        re = round(z.real, nd)
        im = 0.0 if abs(z.imag) < tol else round(z.imag, nd)
        out.append(complex(re) if im == 0.0 else complex(re, im))
    return tuple(out)

def solve_quadratic(a,b,c):

  if a == 0:
    return cleaning_function([-c/b])

  discriminant = b**2-4*a*c
  d = (discriminant - 2*a**2)/(2*a**2)


  if abs(d) < 1:
    alpha1 = 0.5*cmath.acos(d)
    alpha2 = 0.5*(cmath.acos(d) + 2*math.pi)

    root1 = cmath.cos(alpha1) - b/(2*a)
    root2 = cmath.cos(alpha2) - b/(2*a)

    return cleaning_function([root1,root2])

  else:
      alpha1 = 0.5*cmath.acosh(d+0j)
      alpha2 = 0.5*(cmath.acosh(d+0j) + 2j*math.pi)

      root1 = cmath.cosh(alpha1) - b/(2*a)
      root2 = cmath.cosh(alpha2) - b/(2*a)

      return cleaning_function([root1,root2])


def main():
    tests = [
        (1, 0, 1),  # roots of x^2 + 1 = 0 are [i, -i]
        (1,4,-8),    #  roots of x^2 +4x+8=0 are [1.46, -5.46]
        (1, -8, 15),   # roots [3.0, 5.0]
        (0,2,6),
        (1,4,40)
    ]
    for a, b, c in tests:
        roots = solve_quadratic(a, b, c)
        print(f"solve_quadratic({a}, {b}, {c}) -> {roots}")

if __name__ == "__main__":
    main()
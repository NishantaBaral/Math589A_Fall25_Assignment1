import math 
import cubic_solver as cubic_solver
import quartic_solver as quartic_solver 
import quadratic_solver as quadratic_solver
import two_ladders as two_ladders


def same_roots(actual, expected, tol=1e-12):
    A = sorted(map(complex, actual), key=lambda z: (z.real, z.imag))
    E = sorted(map(complex, expected), key=lambda z: (z.real, z.imag))
    assert len(A) == len(E)
    for a, e in zip(A, E):
        assert abs(a - e) <= tol


def test_quadratic_real_distinct():
    # x^2 - 4x + 3 = 0  ->  {1, 3}
    roots = quadratic_solver.solve_quadratic(1, -4, 3)
    same_roots(roots, [1, 3])


def test_cubic():
    # x^3 + 1 = 0 → {-1, 0.5±i*sqrt(3)/2}
    roots = cubic_solver.solve_cubic(1, 0, 0, 1)
    expected = [
        -1,
        0.5 + (math.sqrt(3)/2)*1j,
        0.5 - (math.sqrt(3)/2)*1j,
    ]
    same_roots(roots, expected)

def test_quartic():
    # (x-1)(x-2)(x-3)(x-4)
    roots = quartic_solver.solve_quartic(1, -10, 35, -50, 24)  # use your module's function name
    same_roots(roots, [1, 2, 3, 4])

def test_two_ladders():
    # L1=40, L2=30, h=10 → w ≈ 26.0328775442
    widths = two_ladders.two_ladders_width(40.0, 30.0, 10.0)
    assert any(math.isclose(w, 26.03287754423185, rel_tol=1e-9, abs_tol=1e-9) for w in widths)
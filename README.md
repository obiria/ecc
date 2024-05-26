Elliptic Curve Cryptography Library
Overview
This library provides basic operations and utility functions for working with elliptic curves and elliptic curve cryptography (ECC). It includes functionality for point addition, point doubling, and scalar multiplication on elliptic curves. The library is implemented in C and uses the GNU MP (GMP) library for handling large integers.

Features
Initialization and memory management for elliptic curves and points
Point addition, doubling, and subtraction on elliptic curves
Scalar multiplication using the double-and-add algorithm
Utility functions for testing commutativity, associativity, identity element, and point at infinity
Integration with the GMP library for arbitrary-precision arithmetic
Prerequisites
GMP library installed on your system. This library can be installed via package managers like Homebrew on macOS:

sh
Copy code
brew install gmp
Installation
Clone the repository:

sh
Copy code
git clone https://github.com/obiria/ecc/elliptic-curve-lib.git
cd elliptic-curve-lib
Compile the library:


sh
Copy code
g++ -I/opt/homebrew/Cellar/gmp/6.3.0/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib -lgmp -o main elliptic_curve.cpp
Usage
Initialization
To initialize a new elliptic curve and points on the curve:

c
Copy code
#include "elliptic_curve.h"

int main() {
    mpz_t a, b, p, n, x, y;
    mpz_init_set_str(a, "2", 10);
    mpz_init_set_str(b, "3", 10);
    mpz_init_set_str(p, "17", 10);
    mpz_init_set_str(n, "19", 10);
    mpz_init_set_str(x, "5", 10);
    mpz_init_set_str(y, "1", 10);

    Curve* curve = new_curve(a, b, p, n);
    Point* point = new_point(curve, x, y);

    // Use the point and curve...

    clear_point(point);
    clear_curve(curve);
    mpz_clears(a, b, p, n, x, y, NULL);

    return 0;
}
Point Operations
Perform point addition, doubling, and multiplication:

c
Copy code
Point* sum = point_add(point1, point2);
Point* doubled = point_double(point);
Point* product = point_multiplication(k, point);

// Free the memory when done
clear_point(sum);
clear_point(doubled);
clear_point(product);
Utility Functions
Test properties and operations on points:

c
Copy code
test_commutativity(point1, point2);
test_associativity(point1, point2, point3);
test_identity(point, point_at_infinity(curve));
test_point_operations(curve, base_point, k);
test_point_subtraction(point1, point2);
test_is_infinity(curve);
API Reference
Structures
Curve: Represents an elliptic curve with parameters a, b, p, and n.
Point: Represents a point on an elliptic curve with coordinates x, y, and a reference to the curve.
Functions
Curve* new_curve(const mpz_t a, const mpz_t b, const mpz_t p, const mpz_t n): Initialize a new elliptic curve.
void clear_curve(Curve* curve): Clear memory allocated for an elliptic curve.
Point* new_point(Curve* curve, const mpz_t x, const mpz_t y): Initialize a new point on the elliptic curve.
void clear_point(Point* point): Clear memory allocated for a point.
bool is_infinity(const Point* p): Check if a point is at infinity.
bool points_equal(const Point* p1, const Point* p2): Check if two points are equal.

#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>

// Curve represents an elliptic curve with parameters (a, b, p, n, G).
typedef struct {
    mpz_t a, b, p, n;  // Curve coefficients and parameters
} Curve;

// Point represents a point on an elliptic curve.
typedef struct {
    mpz_t x, y;  // Coordinates of the point
    Curve* curve; // The elliptic curve this point belongs to
} Point;

// Initialize a new point on the elliptic curve.
Point* new_point(Curve* curve, const mpz_t x, const mpz_t y);

// Clear memory for a point.
void clear_point(Point* point);

// Initialize a curve with given parameters.
Curve* new_curve(const mpz_t a, const mpz_t b, const mpz_t p, const mpz_t n);

// Clear memory for a curve.
void clear_curve(Curve* curve);

// Check if a point is at infinity.
bool is_infinity(const Point* p);

// Check if two points are equal.
bool points_equal(const Point* p1, const Point* p2);

// Point addition on the elliptic curve.
Point* point_add(const Point* p1, const Point* p2);

// Point doubling on the elliptic curve.
Point* point_double(const Point* p);

// Point multiplication on the elliptic curve using double-and-add algorithm.
Point* point_multiplication(const mpz_t k, const Point* p);

// Subtract one point from another on the elliptic curve.
Point* point_subtract(const Point* p1, const Point* p2);

// Function to create the point at infinity
Point* point_at_infinity(Curve* curve);

// Function to check if the addition is commutative: P + Q = Q + P
void test_commutativity(const Point* P, const Point* Q);


// Function to check if the addition is associative: (P + Q) + R = P + (Q + R)
void test_associativity(const Point* P, const Point* Q, const Point* R);

// Function to check if the identity element (infinity) works correctly: P + O = O + P = P
void test_identity(const Point* P, const Point* O) ;


// Test function to validate point addition, doubling, and multiplication
void test_point_operations(Curve* curve, const Point* G, const mpz_t k);
// Function to test point subtraction: P - Q = P + (-Q)
void test_point_subtraction(const Point* P, const Point* Q);


// Function to test if a point is at infinity
void test_is_infinity(Curve* curve);
#endif /* ELLIPTIC_CURVE_H */


// g++ -I/opt/homebrew/Cellar/gmp/6.3.0/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib -lgmp -o main elliptic_curve.cpp

// g++ -shared -o libelliptic_curve.so -fPIC elliptic_curve.cpp -I/opt/homebrew/Cellar/gmp/6.3.0/include -L/opt/homebrew/Cellar/gmp/6.3.0/lib  -lgmp
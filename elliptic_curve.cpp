#include "elliptic_curve.h"
#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

// // Curve represents an elliptic curve with parameters (a, b, p, n, G).
// typedef struct {
//     mpz_t a, b, p, n;  // Curve coefficients and parameters
// } Curve;

// // Point represents a point on an elliptic curve.
// typedef struct {
//     mpz_t x, y;  // Coordinates of the point
//     Curve* curve; // The elliptic curve this point belongs to
// } Point;

// Initialize a new point on the elliptic curve.
Point* new_point(Curve* curve, const mpz_t x, const mpz_t y) {
    Point* point = (Point*) malloc(sizeof(Point));
    mpz_init_set(point->x, x);
    mpz_init_set(point->y, y);
    point->curve = curve;
    return point;
}

// Clear memory for a point.
void clear_point(Point* point) {
    mpz_clear(point->x);
    mpz_clear(point->y);
    free(point);
}

// Initialize a curve with given parameters.
Curve* new_curve(const mpz_t a, const mpz_t b, const mpz_t p, const mpz_t n) {
    Curve* curve = (Curve*) malloc(sizeof(Curve));
    mpz_init_set(curve->a, a);
    mpz_init_set(curve->b, b);
    mpz_init_set(curve->p, p);
    mpz_init_set(curve->n, n);
    return curve;
}

// Clear memory for a curve.
void clear_curve(Curve* curve) {
    mpz_clear(curve->a);
    mpz_clear(curve->b);
    mpz_clear(curve->p);
    mpz_clear(curve->n);
    free(curve);
}

// Check if a point is at infinity.
bool is_infinity(const Point* p) {
    return mpz_cmp_ui(p->x, 0) == 0 && mpz_cmp_ui(p->y, 0) == 0;
}

// Check if two points are equal.
bool points_equal(const Point* p1, const Point* p2) {
    return mpz_cmp(p1->x, p2->x) == 0 && mpz_cmp(p1->y, p2->y) == 0;
}

// Point addition on the elliptic curve.
Point* point_add(const Point* p1, const Point* p2) {
    if (is_infinity(p1)) return new_point(p2->curve, p2->x, p2->y);
    if (is_infinity(p2)) return new_point(p1->curve, p1->x, p1->y);

    Curve* curve = p1->curve;
    mpz_t lambda, inv, x3, y3, temp;

    mpz_inits(lambda, inv, x3, y3, temp, NULL);

    if (points_equal(p1, p2)) {
        // Point doubling
        mpz_mul(lambda, p1->x, p1->x);              // lambda = x1^2
        mpz_mul_ui(lambda, lambda, 3);              // lambda = 3 * x1^2
        mpz_add(lambda, lambda, curve->a);          // lambda = 3 * x1^2 + a
        mpz_mul_ui(inv, p1->y, 2);                  // inv = 2 * y1
    } else {
        // Point addition
        mpz_sub(lambda, p2->y, p1->y);              // lambda = y2 - y1
        mpz_sub(inv, p2->x, p1->x);                 // inv = x2 - x1
    }

    mpz_invert(inv, inv, curve->p);                 // inv = inv(inv) mod p
    mpz_mul(lambda, lambda, inv);                   // lambda = lambda * inv
    mpz_mod(lambda, lambda, curve->p);              // lambda = lambda mod p

    mpz_mul(x3, lambda, lambda);                    // x3 = lambda^2
    mpz_sub(x3, x3, p1->x);                         // x3 = lambda^2 - x1
    mpz_sub(x3, x3, p2->x);                         // x3 = lambda^2 - x1 - x2
    mpz_mod(x3, x3, curve->p);                      // x3 = x3 mod p

    mpz_sub(y3, p1->x, x3);                         // y3 = x1 - x3
    mpz_mul(y3, y3, lambda);                        // y3 = lambda * (x1 - x3)
    mpz_sub(y3, y3, p1->y);                         // y3 = lambda * (x1 - x3) - y1
    mpz_mod(y3, y3, curve->p);                      // y3 = y3 mod p

    Point* result = new_point(curve, x3, y3);

    mpz_clears(lambda, inv, x3, y3, temp, NULL);

    return result;
}

// Point doubling on the elliptic curve.
Point* point_double(const Point* p) {
    if (is_infinity(p)) return new_point(p->curve, p->x, p->y);

    Curve* curve = p->curve;
    mpz_t lambda, inv, x3, y3;

    mpz_inits(lambda, inv, x3, y3, NULL);

    // lambda = (3 * x^2 + a) / (2 * y)
    mpz_mul(lambda, p->x, p->x);              // lambda = x^2
    mpz_mul_ui(lambda, lambda, 3);            // lambda = 3 * x^2
    mpz_add(lambda, lambda, curve->a);        // lambda = 3 * x^2 + a
    mpz_mul_ui(inv, p->y, 2);                 // inv = 2 * y
    mpz_invert(inv, inv, curve->p);           // inv = inv(2 * y) mod p
    mpz_mul(lambda, lambda, inv);             // lambda = (3 * x^2 + a) * inv
    mpz_mod(lambda, lambda, curve->p);        // lambda = lambda mod p

    // x3 = lambda^2 - 2 * x
    mpz_mul(x3, lambda, lambda);              // x3 = lambda^2
    mpz_submul_ui(x3, p->x, 2);               // x3 = lambda^2 - 2 * x
    mpz_mod(x3, x3, curve->p);                // x3 = x3 mod p

    // y3 = lambda * (x - x3) - y
    mpz_sub(y3, p->x, x3);                    // y3 = x - x3
    mpz_mul(y3, y3, lambda);                  // y3 = lambda * (x - x3)
    mpz_sub(y3, y3, p->y);                    // y3 = lambda * (x - x3) - y
    mpz_mod(y3, y3, curve->p);                // y3 = y3 mod p

    Point* result = new_point(curve, x3, y3);

    mpz_clears(lambda, inv, x3, y3, NULL);

    return result;
}

// Point multiplication on the elliptic curve using double-and-add algorithm.
Point* point_multiplication(const mpz_t k, const Point* p) {
    if (mpz_cmp_ui(k, 0) == 0 || is_infinity(p)) {
        return new_point(p->curve, p->x, p->y);
    }

    Point* result = new_point(p->curve, p->x, p->y);
    Point* addend = new_point(p->curve, p->x, p->y);

    mpz_t tmp_k;
    mpz_init_set(tmp_k, k);

    while (mpz_cmp_ui(tmp_k, 0) > 0) {
        if (mpz_odd_p(tmp_k)) {
            Point* temp = point_add(result, addend);
            clear_point(result);
            result = temp;
        }

        Point* temp = point_double(addend);
        clear_point(addend);
        addend = temp;

        mpz_fdiv_q_2exp(tmp_k, tmp_k, 1);  // tmp_k = tmp_k / 2
    }

    mpz_clear(tmp_k);
    clear_point(addend);

    return result;
}

// Subtract one point from another on the elliptic curve.
Point* point_subtract(const Point* p1, const Point* p2) {
    if (is_infinity(p2)) return new_point(p1->curve, p1->x, p1->y);

    mpz_t negY;
    mpz_init(negY);

    mpz_neg(negY, p2->y);
    mpz_mod(negY, negY, p1->curve->p);

    Point* negP2 = new_point(p2->curve, p2->x, negY);
    Point* result = point_add(p1, negP2);

    mpz_clear(negY);
    clear_point(negP2);

    return result;
}

// Function to create the point at infinity
Point* point_at_infinity(Curve* curve) {
    mpz_t zero;
    mpz_init_set_ui(zero, 0);
    Point* infinity = new_point(curve, zero, zero);
    mpz_clear(zero);
    return infinity;
}

// Function to check if the addition is commutative: P + Q = Q + P
void test_commutativity(const Point* P, const Point* Q) {
    Point* P_plus_Q = point_add(P, Q);
    Point* Q_plus_P = point_add(Q, P);
    assert(points_equal(P_plus_Q, Q_plus_P));
    clear_point(P_plus_Q);
    clear_point(Q_plus_P);
}

// Function to check if the addition is associative: (P + Q) + R = P + (Q + R)
void test_associativity(const Point* P, const Point* Q, const Point* R) {
    Point* P_plus_Q = point_add(P, Q);
    Point* P_plus_Q_plus_R = point_add(P_plus_Q, R);
    Point* Q_plus_R = point_add(Q, R);
    Point* P_plus_Q_plus_R_alt = point_add(P, Q_plus_R);
    assert(points_equal(P_plus_Q_plus_R, P_plus_Q_plus_R_alt));
    clear_point(P_plus_Q);
    clear_point(P_plus_Q_plus_R);
    clear_point(Q_plus_R);
    clear_point(P_plus_Q_plus_R_alt);
}

// Function to check if the identity element (infinity) works correctly: P + O = O + P = P
void test_identity(const Point* P, const Point* O) {
    Point* P_plus_O = point_add(P, O);
    Point* O_plus_P = point_add(O, P);
    assert(points_equal(P, P_plus_O));
    assert(points_equal(P, O_plus_P));
    clear_point(P_plus_O);
    clear_point(O_plus_P);
}


// Test function to validate point addition, doubling, and multiplication
void test_point_operations(Curve* curve, const Point* G, const mpz_t k) {
    // Point doubling: 2G = G + G
    Point* doubled_G = point_double(G);
    Point* G_plus_G = point_add(G, G);
    assert(points_equal(doubled_G, G_plus_G));
    clear_point(doubled_G);
    clear_point(G_plus_G);

    // Point multiplication: kG
    Point* kG = point_multiplication(k, G);
    mpz_t two;
    mpz_init_set_ui(two, 2);
    Point* twoG = point_multiplication(two, G);
    assert(points_equal(kG, twoG));
    clear_point(kG);
    clear_point(twoG);
    mpz_clear(two);
}

// Function to test point subtraction: P - Q = P + (-Q)
void test_point_subtraction(const Point* P, const Point* Q) {
    Point* P_minus_Q = point_subtract(P, Q);
    mpz_t negY;
    mpz_init(negY);
    mpz_neg(negY, Q->y);
    mpz_mod(negY, negY, P->curve->p);
    Point* negQ = new_point(Q->curve, Q->x, negY);
    Point* P_plus_negQ = point_add(P, negQ);
    assert(points_equal(P_minus_Q, P_plus_negQ));
    clear_point(P_minus_Q);
    clear_point(negQ);
    clear_point(P_plus_negQ);
    mpz_clear(negY);
}

// Function to test if a point is at infinity
void test_is_infinity(Curve* curve) {
    Point* infinity = point_at_infinity(curve);
    assert(is_infinity(infinity));

    mpz_t x, y;
    mpz_inits(x, y, NULL);
    mpz_set_str(x, "1", 10);
    mpz_set_str(y, "1", 10);

    Point* not_infinity = new_point(curve, x, y);
    assert(!is_infinity(not_infinity));

    clear_point(infinity);
    clear_point(not_infinity);
    mpz_clears(x, y, NULL);
}

// Example usage
int main() {
    mpz_t a, b, p, n;
    mpz_inits(a, b, p, n, NULL);

    mpz_set_str(a, "0", 10);
    mpz_set_str(b, "7", 10);
    mpz_set_str(p, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    mpz_set_str(n, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16);

    Curve* curve = new_curve(a, b, p, n);

    mpz_t x, y, k;
    mpz_inits(x, y, k, NULL);

    mpz_set_str(x, "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
    mpz_set_str(y, "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16);
    mpz_set_str(k, "2", 10);

    Point* G = new_point(curve, x, y);
    Point* result = point_multiplication(k, G);

    gmp_printf("Result: (%Zd, %Zd)\n", result->x, result->y);


    Point* O = point_at_infinity(curve);

    test_commutativity(G, G);
    test_associativity(G, G, G);
    test_identity(G, O);
    test_point_operations(curve, G, k);
    test_point_subtraction(G, G);
    test_is_infinity(curve);
    printf("All tests passed!\n");

    clear_point(G);
    clear_point(O);
    clear_curve(curve);

    mpz_clears(a, b, p, n, x, y, k, NULL);

    return 0;
}

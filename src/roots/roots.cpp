#include "roots.hpp"
#include <cmath>

const double TOL = 1e-6;
const int MAX_ITER = 1000000;


bool bisection(std::function<double(double)> f, double a, double b, double *root)

{
    if (f(a) * f(b) > 0) {
        return false;
    }

    double mid = a;

    for (int i = 0; i < MAX_ITER; i++)
    {
        mid = 0.5 * (a + b);

        if (std::abs(f(mid)) < TOL || (b - a) / 2.0 < TOL) {
            *root = mid;
            return true;
        }

        if (f(a) * f(mid) < 0) {
            b = mid;
        }
        else {
            a = mid;
        }
    }

    *root = mid;
    return true;
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root)
{
    if (f(a) * f(b) > 0) {
        return false;
    }

    double c = a;

    for (int i = 0; i < MAX_ITER; i++) {
        double fa = f(a);
        double fb = f(b);

        if (std::abs(fb - fa) < TOL) {
            return false;
        }

        c = (a * fb - b * fa) / (fb - fa); //Formula

        if (std::abs(f(c)) < TOL) {
            *root = c;
            return true;
        }

        if (fa * f(c) < 0) {
            b = c;
        }
        else {
            a = c;
        }
    }

    *root = c;
    return true;
}

bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g,
                    double a, double b, double c, double *root)
{

    double x = c;

    for (int i = 0; i < MAX_ITER; i++) {
        double fx = f(x);
        double gx = g(x);

        if (std::abs(gx) < TOL)
            return false;

        double x_new = x - fx / gx; //formula for next guess

        if (x_new < a || x_new > b) //check interval using new guess
            return false;

        if (std::abs(x_new - x) < TOL) { // lets the for loop stop if x_new within TOL
            *root = x_new;
            return true;
        }

        x = x_new;
    }

    *root = x;
    return true;
}

bool secant(std::function<double(double)> f, double a, double b, double c, double *root)
{

    double x0 = a;
    double x1 = b;
    double x2;

    for (int i = 0; i < MAX_ITER; i++)
    {
        double f0 = f(x0);
        double f1 = f(x1);

        if (std::abs(f1 - f0) < TOL) // Prevent division by zero
            return false;

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0); // formula for new estimate of x2

        if (std::abs(x2 - x1) < TOL)
        {
            *root = x2;
            return true;
        }

        x0 = x1;
        x1 = x2;
    }

    *root = x2;
    return true;
}
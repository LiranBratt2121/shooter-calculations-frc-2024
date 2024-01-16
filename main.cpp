#include <iostream>
#include <cmath>
#include <algorithm>

/*
    # This program calculates the desired velocity and angle for an object to land inside the target
    # during the FRC 2024 season!

    The program uses the Gauss-Jordan elimination algorithm to solve the system of equations.

    ## Equations:
    v0 * cos(θ) + v * cos(α) = 2Δx
    v0 * sin(θ) + v * sin(α) = 2Δy

    ## Constants:
    - Δx: Representing the distance difference from the shooter to the target:
    - Δy: Representing the height difference from the shooter to the target:
    - α:  Representing the angle in which the object should land at inside the target: 110deg
    - v0: Representing the initial velocity of the object:

    ## Variables:
    - v: Representing the desired velocity of the object:
    - θ: Representing the angle of the shooter:
*/

const int MAX_EQUATIONS = 2;
const double PI = 3.14159265358979323846;

const double TO_DEGREES = 180.0 / PI;
const int VELOCITY_INDEX = 0;
const int THETA_INDEX = 1;

void partial_pivot(double coefficients[MAX_EQUATIONS][MAX_EQUATIONS + 1], int n);
void back_substitute(double coefficients[MAX_EQUATIONS][MAX_EQUATIONS + 1], int n, double x[MAX_EQUATIONS]);

int main()
{
    const int n = MAX_EQUATIONS;

    // Constants and variables **TEMP**
    double dx = 3.5;
    double dy = 3;
    double v0 = 15;
    double alpha = 110;

    // Coefficient matrix initialization
    double coefficients[MAX_EQUATIONS][MAX_EQUATIONS + 1] = {
        {v0, cos(alpha) * TO_DEGREES, 2 * dx},
        {v0, sin(alpha) * TO_DEGREES, 2 * dy},
    };

    // Solution vector
    double x[MAX_EQUATIONS];

    // Partial pivot and back substitution
    partial_pivot(coefficients, n);
    back_substitute(coefficients, n, x);

    // Output
    double velocity = x[VELOCITY_INDEX];
    double theta = cos(x[THETA_INDEX]) * TO_DEGREES;

    std::cout << "Solution for the system:\n";
    std::cout << "v = " << velocity << std::endl;
    std::cout << "θ = " << theta << std::endl;

    return 0;
}

void partial_pivot(double coefficients[MAX_EQUATIONS][MAX_EQUATIONS + 1], int n)
{
    for (int i = 0; i < n; i++)
    {
        int pivot_row = i;

        for (int j = i + 1; j < n; j++)
        {
            if (std::abs(coefficients[j][i]) > std::abs(coefficients[pivot_row][i]))
            {
                pivot_row = j;
            }
        }

        if (pivot_row != i)
        {
            for (int j = i; j <= n; j++)
            {
                std::swap(coefficients[i][j], coefficients[pivot_row][j]);
            }
        }

        for (int j = i + 1; j < n; j++)
        {
            double factor = coefficients[j][i] / coefficients[i][i];

            for (int k = i; k <= n; k++)
            {
                coefficients[j][k] -= factor * coefficients[i][k];
            }
        }
    }
}

void back_substitute(double coefficients[MAX_EQUATIONS][MAX_EQUATIONS + 1], int n, double x[MAX_EQUATIONS])
{
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;

        for (int j = i + 1; j < n; j++)
        {
            sum += coefficients[i][j] * x[j];
        }

        x[i] = (coefficients[i][n] - sum) / coefficients[i][i];
    }
}

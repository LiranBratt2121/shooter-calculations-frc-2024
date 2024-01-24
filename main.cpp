#include <iostream>
#include <cmath>
#include <algorithm>
#include "time.h"
#include <unordered_map>

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


class SystemSolver
{
public:
    SystemSolver(double dx, double dy, double v0, double alpha);
    std::unordered_map<std::string, std::double_t> solveSystem();

private:
    const int MAX_EQUATIONS = 2;
    const double PI = 3.14159265358979323846;
    const double TO_DEGREES = 180.0 / PI;
    const int VELOCITY_INDEX = 0;
    const int THETA_INDEX = 1;

    double dx;
    double dy;
    double v0;
    double alpha;

    double coefficients[2][3];
    double x[2];

    void partial_pivot();
    void back_substitute();
};

SystemSolver::SystemSolver(double dx, double dy, double v0, double alpha)
    : dx(dx), dy(dy), v0(v0), alpha(alpha)
{
    coefficients[0][0] = v0;
    coefficients[0][1] = cos(alpha * TO_DEGREES);

    coefficients[1][0] = v0;
    coefficients[1][1] = sin(alpha * TO_DEGREES);
    coefficients[1][2] = 2 * dy;
}

void SystemSolver::partial_pivot()
{
    for (int i = 0; i < MAX_EQUATIONS; i++)
    {
        int pivot_row = i;

        for (int j = i + 1; j < MAX_EQUATIONS; j++)
        {
            if (std::abs(coefficients[j][i]) > std::abs(coefficients[pivot_row][i]))
            {
                pivot_row = j;
            }
        }

        if (pivot_row != i)
        {
            for (int j = i; j <= MAX_EQUATIONS; j++)
            {
                std::swap(coefficients[i][j], coefficients[pivot_row][j]);
            }
        }

        for (int j = i + 1; j < MAX_EQUATIONS; j++)
        {
            double factor = coefficients[j][i] / coefficients[i][i];

            for (int k = i; k <= MAX_EQUATIONS; k++)
            {
                coefficients[j][k] -= factor * coefficients[i][k];
            }
        }
    }
}

void SystemSolver::back_substitute()
{
    for (int i = MAX_EQUATIONS - 1; i >= 0; i--)
    {
        double sum = 0;

        for (int j = i + 1; j < MAX_EQUATIONS; j++)
        {
            sum += coefficients[i][j] * x[j];
        }

        x[i] = (coefficients[i][MAX_EQUATIONS] - sum) / coefficients[i][i];
    }
}

std::unordered_map<std::string, std::double_t> SystemSolver::solveSystem()
{
    partial_pivot();
    back_substitute();

    double velocity = x[VELOCITY_INDEX];
    double theta = cos(x[THETA_INDEX]) * TO_DEGREES;

    std::cout << "Solution for the system:\n";
    std::cout << "v = " << velocity << std::endl;
    std::cout << "θ = " << theta << std::endl;
    std::cout << "---------------------------------\n";

    std::unordered_map<std::string, std::double_t> results;

    results["v"] = velocity;
    results["theta"] = theta;

    return results;
}

int main()
{
    SystemSolver solver1(3.5, 3.0, 15.0, 110.0);
    auto result1 = solver1.solveSystem();

    std::cout << result1["v"] << std::endl;
    std::cout << result1["theta"] << std::endl;
    
    return 0;
}

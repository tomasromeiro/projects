/*******************************************************
 * Student_ID: 24022849
 * Finite-Difference Methods for Vanilla Option Pricing
 *
 * This program implements three separate finite-difference schemes:
 *   1. Explicit Scheme
 *   2. Fully Implicit Scheme
 *   3. Crank–Nicolson Scheme
 *
 * For the Crank–Nicolson scheme, the program computes two prices:
 *   a) Using a user-supplied theta value
 *   b) Using the special theta defined as: theta_special = 0.5 + (1/12)*((dS)^2/dt)
 *
 * Each method uses a switch statement to enforce the early exercise 
 * condition for American options.
 *
 * The asset-price grid is defined from 0 to Smax, with Smax set to 2 * max(S0, E).
 *
 * The program prints the option prices for each method.
 *  
 * Example input:
 *   S0 = 100, E = 100, T = 1, M = 500, N = 15000,
 *   Option type: C (Call) or P (Put)
 *   Exercise style: E (European) or A (American)
 *   sigma = 0.20, r = 0.05, q = 0
 * Example output (for an European Call Option with the parameters above):
 *   Explicit FD Price: 10.453360
 *   Implicit FD Price: 10.453226
 *   Crank-Nicolson FD Price (user theta = 0.500000): 10.453293
 *   Crank-Nicolson FD Price (special theta = 200.500000): 10.426359
 *
 * 
 * See comments on each exercise below:
 * 
 * 1. Dividends comments:
 * In this implementation, dividends are incorporated through the boundary conditions (and implicitly via the evolution of the option values). 
 * For example, for a call option, the upper boundary condition is set as:
 *   V(Smax, t) = Smax * exp(–q * (T – t)) – E * exp(–r * (T – t))
 * which discounts the asset price by exp(–q * (T – t)) to reflect the continuous dividend yield. 
 * For put options, the lower boundary condition is:
 *   V(0, t) = E * exp(–r * (T – t)).
 * Since our payoff function is applied uniformly across all time steps (i.e., there is one terminal payoff function), 
 * we conclude that the effect of dividends is captured solely via these boundary adjustments. 
 * This means that for calls, dividends reduce the effective asset value, thereby lowering call prices (and, conversely, increasing put prices).
 * 
 * 2. Special case discussion:
 * The special case adds a correction term to the standard Crank–Nicolson value of 0.5. 
 * The factor 1/12 comes from a truncation error analysis of the finite-difference discretization. 
 * This additional term adjusts the weighting between the implicit and explicit parts in order to reduce numerical oscillations and improve stability for certain grid choices.
 * In effect it nudges the scheme towards the implicit side. This can help reduce numerical oscillations and improve stability without fully sacrificing the second-order accuracy typically associated with Crank–Nicolson. 
 * The larger theta is, the closer the scheme is to being fully implicit, and the special theta attempts to balance stability and accuracy by moving slightly toward the implicit end of the spectrum.
 * 
 *******************************************************/

 #include <iostream>
 #include <vector>
 #include <cmath>
 #include <algorithm>
 #include <iomanip>
 using namespace std;
 
 /*
  * Linear interpolation to estimate the option price at S = S0.
  */
 double interpolateAtS0(const vector<double>& values, double S0, double dS) {
     double idx = S0 / dS;
     int i_low = static_cast<int>(floor(idx));
     if (fabs(idx - i_low) < 1e-12)
         return values[i_low];
     int i_high = i_low + 1;
     if (i_high >= (int)values.size())
         return values.back();
     double S_low = i_low * dS;
     double S_high = i_high * dS;
     return values[i_low] + (values[i_high] - values[i_low]) *
                            ((S0 - S_low) / (S_high - S_low));
 }
 
 /*
  * Terminal payoff for vanilla options.
  *   Call: max(S - E, 0)
  *   Put:  max(E - S, 0)
  */
 double terminalPayoff(double S, double E, char optType) {
     if (optType == 'C' || optType == 'c')
         return max(S - E, 0.0);
     else
         return max(E - S, 0.0);
 }
 
 /*
  * Intrinsic value (used for American options).
  */
 double intrinsicValue(double S, double E, char optType) {
     return terminalPayoff(S, E, optType);
 }
 
 /*
  * Set boundary conditions for vanilla options.
  * For a call:
  *   V(0,t) = 0,
  *   V(Smax,t) = Smax*exp(-q*(T-t)) - E*exp(-r*(T-t))
  * For a put:
  *   V(0,t) = E*exp(-r*(T-t)),
  *   V(Smax,t) = 0.
  */
 void setBoundaryConditions(vector<double>& vals, double Smax, double E,
                            double t, double T, char optType,
                            double r, double q) {
     if (optType == 'C' || optType == 'c') {
         vals[0] = 0.0;
         vals.back() = Smax * exp(-q * (T - t)) - E * exp(-r * (T - t));
     } else {
         vals[0] = E * exp(-r * (T - t));
         vals.back() = 0.0;
     }
 }
 
 /*
  * Explicit Finite-Difference Scheme
  */
 double explicitFD(double S0, double E, double r, double sigma, double T,
                   double Smax, int M, int N, char optType, char exerciseStyle,
                   double q) {
     double dS = Smax / M;
     double dt = T / N;
     vector<double> oldVals(M + 1), newVals(M + 1);
 
     // Terminal condition.
     for (int i = 0; i <= M; i++) {
         double S = i * dS;
         oldVals[i] = terminalPayoff(S, E, optType);
     }
 
     // Time marching backward.
     for (int n = 0; n < N; n++) {
         double tNext = (n + 1) * dt;
         setBoundaryConditions(newVals, Smax, E, tNext, T, optType, r, q);
         for (int i = 1; i < M; i++) {
             double A = 0.5 * dt * (sigma * sigma * i * i - (r - q) * i);
             double B = 1.0 - dt * (sigma * sigma * i * i + r);
             double C = 0.5 * dt * (sigma * sigma * i * i + (r - q) * i);
             newVals[i] = A * oldVals[i - 1] + B * oldVals[i] + C * oldVals[i + 1];
         }
         // Enforce early exercise if American.
         switch(exerciseStyle) {
             case 'A': case 'a':
                 for (int i = 1; i < M; i++) {
                     double S = i * dS;
                     newVals[i] = max(newVals[i], intrinsicValue(S, E, optType));
                 }
                 break;
             case 'E': case 'e':
             default:
                 break;
         }
         oldVals = newVals;
     }
     return interpolateAtS0(oldVals, S0, dS);
 }
 
 /*
  * Fully Implicit Finite-Difference Scheme
  */
 double implicitFD(double S0, double E, double r, double sigma, double T,
                   double Smax, int M, int N, char optType, char exerciseStyle,
                   double q) {
     double dS = Smax / M;
     double dt = T / N;
     vector<double> oldVals(M + 1), newVals(M + 1);
     vector<double> a(M + 1), b(M + 1), c(M + 1), d(M + 1);
 
     // Terminal condition.
     for (int i = 0; i <= M; i++) {
         double S = i * dS;
         oldVals[i] = terminalPayoff(S, E, optType);
     }
 
     // Time marching backward.
     for (int n = 0; n < N; n++) {
         double tNext = (n + 1) * dt;
         setBoundaryConditions(newVals, Smax, E, tNext, T, optType, r, q);
         for (int i = 1; i < M; i++) {
             a[i] = -0.5 * dt * (sigma * sigma * i * i - (r - q) * i);
             b[i] = 1.0 + dt * (sigma * sigma * i * i + r);
             c[i] = -0.5 * dt * (sigma * sigma * i * i + (r - q) * i);
             d[i] = oldVals[i];
         }
         d[1]   -= a[1] * newVals[0];
         d[M-1] -= c[M-1] * newVals[M];
         // Forward elimination.
         for (int i = 2; i < M; i++) {
             double m = a[i] / b[i - 1];
             b[i] -= m * c[i - 1];
             d[i] -= m * d[i - 1];
         }
         // Back substitution.
         newVals[M-1] = d[M-1] / b[M-1];
         for (int i = M - 2; i >= 1; i--) {
             newVals[i] = (d[i] - c[i] * newVals[i + 1]) / b[i];
         }
         // Enforce early exercise if American.
         switch(exerciseStyle) {
             case 'A': case 'a':
                 for (int i = 1; i < M; i++) {
                     double S = i * dS;
                     newVals[i] = max(newVals[i], intrinsicValue(S, E, optType));
                 }
                 break;
             case 'E': case 'e':
             default:
                 break;
         }
         oldVals = newVals;
     }
     return interpolateAtS0(oldVals, S0, dS);
 }
 
 /*
  * Crank–Nicolson Finite-Difference Scheme.
  *
  * This function uses the theta-method:
  *   (V_i^(n+1) - V_i^n)/dt = theta*L(V^(n+1)) + (1-theta)*L(V^n)
  *
  * The user-supplied theta is used to compute one CN price.
  * In addition, we compute a second CN price using the special theta:
  *
  *   theta_special = 0.5 + (1/12)*((dS)^2/dt)
  */
 double crankNicolsonFD(double S0, double E, double r, double sigma, double T,
                        double Smax, int M, int N, char optType, char exerciseStyle,
                        double q, double theta) {
     double dS = Smax / M;
     double dt = T / N;
     vector<double> oldVals(M + 1), newVals(M + 1);
     vector<double> a(M + 1), b(M + 1), c(M + 1), d(M + 1);
     
     // Terminal condition.
     for (int i = 0; i <= M; i++) {
         double S = i * dS;
         oldVals[i] = terminalPayoff(S, E, optType);
     }
     
     // Time marching backward.
     for (int n = 0; n < N; n++) {
         double tNext = (n + 1) * dt;
         setBoundaryConditions(newVals, Smax, E, tNext, T, optType, r, q);
         for (int i = 1; i < M; i++) {
             double A = 0.5 * dt * (sigma * sigma * i * i - (r - q) * i);
             double B = dt * (sigma * sigma * i * i + r);
             double C = 0.5 * dt * (sigma * sigma * i * i + (r - q) * i);
             
             a[i] = -theta * A;
             b[i] = 1 + theta * B;
             c[i] = -theta * C;
             d[i] = (1 - theta) * (A * oldVals[i - 1] - B * oldVals[i] + C * oldVals[i + 1])
                    + oldVals[i];
         }
         
         d[1]   -= a[1] * newVals[0];
         d[M-1] -= c[M-1] * newVals[M];
         
         // Solve the tridiagonal system: forward elimination.
         for (int i = 2; i < M; i++) {
             double m = a[i] / b[i - 1];
             b[i] -= m * c[i - 1];
             d[i] -= m * d[i - 1];
         }
         // Back substitution.
         newVals[M-1] = d[M-1] / b[M-1];
         for (int i = M - 2; i >= 1; i--) {
             newVals[i] = (d[i] - c[i] * newVals[i + 1]) / b[i];
         }
         
         // Enforce early exercise if American.
         switch(exerciseStyle) {
             case 'A': case 'a':
                 for (int i = 1; i < M; i++) {
                     double S = i * dS;
                     newVals[i] = max(newVals[i], intrinsicValue(S, E, optType));
                 }
                 break;
             case 'E': case 'e':
             default:
                 break;
         }
         
         oldVals = newVals;
     }
     return interpolateAtS0(oldVals, S0, dS);
 }
 
 /*
  * Main program.
  * It computes and displays option prices using:
  *   1. Explicit FD
  *   2. Fully Implicit FD
  *   3. Crank–Nicolson FD using a user-supplied theta
  *   4. Crank–Nicolson FD using the special theta
  */
 int main() {
     double S0, E, T;
     int M, N;
     char optType, exerciseStyle;
     double sigma, r, q;
     double thetaUser;  // user-supplied theta for CN
     
     // Input parameters.
     cout << "Enter current stock price (S0): ";
     cin >> S0;
     cout << "Enter strike price (E): ";
     cin >> E;
     cout << "Enter time to expiry (T in years): ";
     cin >> T;
     cout << "Enter number of asset price steps (M): ";
     cin >> M;
     cout << "Enter number of time steps (N): ";
     cin >> N;
     cout << "Enter option type (C for Call, P for Put): ";
     cin >> optType;
     cout << "Enter exercise style (E for European, A for American): ";
     cin >> exerciseStyle;
     cout << "Enter volatility (sigma, as a decimal): ";
     cin >> sigma;
     cout << "Enter risk-free rate (r, as a decimal): ";
     cin >> r;
     cout << "Enter dividend yield (q, as a decimal) [enter 0 if none]: ";
     cin >> q;
     
     // For the CN scheme, get a user-supplied theta.

     // See comments on the special theta in the comment section at the start of this script.

     cout << "Enter theta for CN scheme (typically 0.5 for standard CN): ";
     cin >> thetaUser;
     
     // Compute Smax automatically.
     double Smax = 2 * max(S0, E);
     
     // Compute option prices.
     double priceExplicit = explicitFD(S0, E, r, sigma, T, Smax, M, N,
                                       optType, exerciseStyle, q);
     double priceImplicit = implicitFD(S0, E, r, sigma, T, Smax, M, N,
                                       optType, exerciseStyle, q);
     double priceCNUser = crankNicolsonFD(S0, E, r, sigma, T, Smax, M, N,
                                          optType, exerciseStyle, q, thetaUser);
     
     // Compute special theta for CN.
     double dS = Smax / M;
     double dt = T / N;
     double thetaSpecial = 0.5 + (1.0 / 12.0) * ((dS * dS) / dt);
     double priceCNSpecial = crankNicolsonFD(S0, E, r, sigma, T, Smax, M, N,
                                             optType, exerciseStyle, q, thetaSpecial);
     
     cout << fixed << setprecision(6);
     cout << "\n------------------------------------------\n";
     cout << "Explicit FD Price: " << priceExplicit << "\n";
     cout << "Implicit FD Price: " << priceImplicit << "\n";
     cout << "Crank-Nicolson FD Price (user theta = " << thetaUser << "): " << priceCNUser << "\n";
     cout << "Crank-Nicolson FD Price (special theta = " << thetaSpecial << "): " << priceCNSpecial << "\n";
     cout << "------------------------------------------\n";
    
     return 0;
 }

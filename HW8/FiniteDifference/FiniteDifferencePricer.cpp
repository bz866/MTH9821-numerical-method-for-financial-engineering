#include "FiniteDifferencePricer.hpp"


// private compuational utilities

std::function<double(double, double, double, double, double)> FiniteDifferencePricer::compute_V_from_heat_equation_solution()
{
    return [&](double a, double x, double b, double tau, double u_x_tau)
    {
        return std::exp(-a * x - b * tau) * u_x_tau;
    };
}

std::tuple<double, std::size_t, double, double> 
FiniteDifferencePricer::compute_domain_params(const std::size_t M) const 
{
    double d_tau = tau_final_ / M;
    std::size_t N = std::floor( (x_r_ - x_l_) / std::sqrt(d_tau / alpha_temp_) );
    double d_x = (x_r_ - x_l_) / N;
    double alpha = d_tau / std::pow(d_x, 2);
    return std::make_tuple(d_tau, N, d_x,  alpha);
}

std::tuple<std::vector<double>, double> FiniteDifferencePricer::buildMesh(const std::size_t N, const double d_x)
{   
    // build x_mesh
    std::vector<double> x_mesh({x_l_});
    for (std::size_t i = 1; i < N; i++)
    { 
        x_mesh.emplace_back(x_mesh.back() + d_x);
    }

    // identify the interval containing x_compute = ln(S0/K), find i such that 
    // x_i <= x_compute < x_i+1
    double x_compute = std::log(S0_ / K_);
    auto interval_i_plus_1 = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), interval_i_plus_1) - 1;

    return std::make_tuple(x_mesh, interval_i);
}

std::vector<double> FiniteDifferencePricer::computeNextUMesh(
    const double tau, const double alpha, const std::vector<double> &x_mesh, const std::vector<double> &u_mesh)
{
    std::vector<double> next_u_mesh(x_mesh.size()+1);
    // left boundary
    next_u_mesh[0] = boundary_x_l_(tau);
    // middle values of u
    auto u_n_m_plus_1 = [&](const int i)->double
    {
        return alpha * u_mesh[i-1] + (1. - 2. * alpha) * u_mesh[i] + alpha * u_mesh[i+1];
    };
    for (std::size_t i = 1; i < x_mesh.size(); ++i)
    {
        next_u_mesh[i] = u_n_m_plus_1(i);
    }
    //  right boundary
    next_u_mesh[x_mesh.size()+1] = boundary_x_r_(tau);
    
    return next_u_mesh;
}

std::pair<std::vector<double>, std::vector<double>> FiniteDifferencePricer::buildUMeshFromFiniteDifference(
    const double alpha, const std::vector<double> &x_mesh, std::size_t M, double d_tau)
{
    // u_mesh initialization
    std::vector<double> u_mesh(x_mesh.size()+1);
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    std::vector<double> u_mesh_tau_plus_1(x_mesh.size()+1);
    // printVector(u_mesh);
    
    // moving forward from t_1 to t_{M-1}
    for (std::size_t i = 1; i < M; ++i)
    {
        double curr_tau = tau_final_ - d_tau * (M-i);
        u_mesh_tau_plus_1 = computeNextUMesh(curr_tau, alpha, x_mesh, u_mesh);
        std::swap(u_mesh, u_mesh_tau_plus_1);
        // printVector(u_mesh);
    }

    // moving forward to t_M
    u_mesh_tau_plus_1 = computeNextUMesh(tau_final_, alpha, x_mesh, u_mesh);
    // printVector(u_mesh_tau_plus_1);

    return std::make_pair(u_mesh, u_mesh_tau_plus_1);
}

// The function to decompose A and solve for the next u_mesh
// LU decomposition and solver for a tridiagonal matrix
std::vector<double> FiniteDifferencePricer::tridiagonalSolve(
    const std::vector<double>& l,
    const std::vector<double>& d,
    const std::vector<double>& u,
    const std::vector<double>& rhs)
{
    size_t n = d.size();
    std::vector<double> c(n), z(n), u_new(n);

    // Forward sweep (LU decomposition)
    c[0] = d[0];
    z[0] = rhs[0] / c[0];
    for (size_t i = 1; i < n; ++i)
    {
        c[i] = d[i] - l[i] * u[i-1] / c[i-1];
        z[i] = (rhs[i] - l[i] * z[i-1]) / c[i];
    }

    // Back substitution
    u_new[n-1] = z[n-1];
    for (int i = n-2; i >= 0; --i)
    {
        u_new[i] = z[i] - u[i] * u_new[i+1] / c[i];
    }

    return u_new;
}

std::vector<double> FiniteDifferencePricer::computeNextUMesh_BE(
    const double tau, const double alpha, const std::vector<double> &x_mesh, const std::vector<double> &u_mesh)
{
    std::vector<double> next_u_mesh(u_mesh.size());
    double lambda = alpha;

    std::vector<double> l(x_mesh.size(), -lambda);
    std::vector<double> d(x_mesh.size(), 1 + 2*lambda);
    std::vector<double> u(x_mesh.size(), -lambda);

    // Implement boundary conditions
    next_u_mesh[0] = boundary_x_l_(tau);
    next_u_mesh[x_mesh.size()] = boundary_x_r_(tau);

    return tridiagonalSolve(l, d, u, u_mesh);
}

std::pair<std::vector<double>, std::vector<double>> FiniteDifferencePricer::buildUMeshFromFiniteDifference_BE(
    const double alpha, const std::vector<double>& x_mesh, std::size_t M, double d_tau)
{
    // u_mesh initialization
    std::vector<double> u_mesh(x_mesh.size()+1);
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    std::vector<double> u_mesh_tau_plus_1(x_mesh.size()+1);
    printVector(u_mesh);

    // moving forward from t_1 to t_{M-1}
    for (std::size_t i = 1; i < M; ++i)
    {
        double curr_tau = tau_final_ - d_tau * (M-i);
        u_mesh_tau_plus_1 = computeNextUMesh_BE(curr_tau, alpha, x_mesh, u_mesh); // change to Backward Euler
        std::swap(u_mesh, u_mesh_tau_plus_1);
        printVector(u_mesh);
    }

    // moving forward to t_M
    u_mesh_tau_plus_1 = computeNextUMesh_BE(tau_final_, alpha, x_mesh, u_mesh); // change to Backward Euler
    printVector(u_mesh_tau_plus_1);

    return std::make_pair(u_mesh, u_mesh_tau_plus_1);
}


std::vector<double> FiniteDifferencePricer::approximateSecurityValue(
    const std::vector<double> &x_mesh, const std::vector<double>& u_mesh, const double d_x)
{
    std::vector<double> V_approximations;
    // identify the interval containing x_compute = ln(S0/K), find i such that 
    // x_i <= x_compute < x_i+1
    double x_compute = std::log(S0_ / K_);
    auto interval_i_plus_1 = std::upper_bound(x_mesh.cbegin(), x_mesh.cend(), x_compute);
    std::size_t interval_i = std::distance(x_mesh.cbegin(), interval_i_plus_1) - 1;

    // finite difference approximation
    double S_i = K_ * std::exp(x_mesh[interval_i]);
    double S_i_plus_1 = K_ * std::exp(x_mesh[interval_i + 1]);
    double V_i = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_i_plus_1 = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_approx = ((S_i_plus_1 - S0_) * V_i + (S0_ - S_i) * V_i_plus_1) / (S_i_plus_1 - S_i);
    V_approximations.emplace_back(V_approx);

    // linear interpolation approximation
    double u_approx = ((x_mesh[interval_i + 1] - x_compute) * u_mesh[interval_i] + (x_compute - x_mesh[interval_i]) * u_mesh[interval_i + 1]) / (d_x);
    double V_approx_2 = std::exp(-a_ * x_compute - b_ * tau_final_) * u_approx;
    V_approximations.emplace_back(V_approx_2);
     
    return V_approximations;
}

double FiniteDifferencePricer::computeRMSError(const std::vector<double>& x_mesh, const std::vector<double>& u_mesh) 
{    
    // Calculate vector of FD and BS option values
    auto V_mesh_approx_gen = [&](double x, double u)->double {
        return std::exp(-a_ * x - b_ * tau_final_) * u;
    };
    auto V_mesh_exact_gen = [&](double x)->double {
        EuropeanOption option(0., K_ * std::exp(x), K_, T_, sigma_, r_, q_);
        return option.Put();
    };
    
    std::vector<double> V_mesh_approx(x_mesh.size());
    std::vector<double> V_mesh_exact(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.cbegin(), V_mesh_approx.begin(), V_mesh_approx_gen);
    std::transform(x_mesh.cbegin(), x_mesh.cend(), V_mesh_exact.begin(), V_mesh_exact_gen);
    
    // Find RMS
    double error_sq = 0.;
    int error_count = 0;
    auto V_mesh_approx_it = V_mesh_approx.cbegin();
    auto V_mesh_exact_it = V_mesh_exact.cbegin();
    while (V_mesh_exact_it != V_mesh_exact.cend()) 
    {
        double V_BS = *V_mesh_exact_it;
        double V_FD = *V_mesh_approx_it;
        if (V_BS > 0.00001 * S0_) 
    {
            error_count++;
            error_sq += (V_BS - V_FD) * (V_BS - V_FD) / (V_BS * V_BS);
        }
        V_mesh_approx_it++;
        V_mesh_exact_it++;
    }
    double error_RMS = std::sqrt(error_sq / error_count);
    
    return error_RMS;
}

std::vector<double> FiniteDifferencePricer::computeGreeks(
    std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, 
    std::vector<double>& u_mesh_M_minus_1_th, double d_tau, double V_approx)
{    
    // DELTA & GAMMA
    // Get S and V of interest
    double S_i_minus_1 = K_ * std::exp(x_mesh[interval_i - 1]);
    double S_i = K_ * std::exp(x_mesh[interval_i]);
    double S_i_plus_1 = K_ * std::exp(x_mesh[interval_i + 1]);
    double S_i_plus_2 = K_ * std::exp(x_mesh[interval_i + 2]);
    
    double V_smaller = std::exp(-a_ * (x_mesh[interval_i - 1]) - b_ * tau_final_) * u_mesh[interval_i - 1];
    double V_small = std::exp(-a_ * (x_mesh[interval_i]) - b_ * tau_final_) * u_mesh[interval_i];
    double V_large = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * tau_final_) * u_mesh[interval_i + 1];
    double V_larger = std::exp(-a_ * (x_mesh[interval_i + 2]) - b_ * tau_final_) * u_mesh[interval_i + 2];
    
    // Find delta and gamma
    double delta = (V_large - V_small) / (S_i_plus_1 - S_i);
    double gamma = ((V_larger - V_large) / (S_i_plus_2 - S_i_plus_1) - (V_small - V_smaller) / (S_i - S_i_minus_1)) / (((S_i_plus_2 + S_i_plus_1) / 2.) - ((S_i + S_i_minus_1) / 2.));
    
    // THETA
    double dt = 2. * d_tau / (sigma_ * sigma_);
    
    // Get V at t = dt
    double V_small_prev = std::exp(-a_ * (x_mesh[interval_i]) - b_ * (tau_final_ - d_tau)) * u_mesh_M_minus_1_th[interval_i];
    double V_large_prev = std::exp(-a_ * (x_mesh[interval_i + 1]) - b_ * (tau_final_ - d_tau)) * u_mesh_M_minus_1_th[interval_i + 1];
    double V_approx_prev = ((S_i_plus_1 - S0_) * V_small_prev + (S0_ - S_i) * V_large_prev) / (S_i_plus_1 - S_i);
    double theta = (V_approx_prev - V_approx) / dt;
    
    return std::vector<double>({delta, gamma, theta});
}

std::vector<double> FiniteDifferencePricer::computeNextUMeshAmerican(
    const double tau, const double alpha, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh)
{
    std::vector<double> next_u_mesh(x_mesh.size());
    // left boundary
    next_u_mesh[0] = boundary_x_l_(tau);
    // middle values
    for (std::size_t i = 1; i < x_mesh.size() - 1; i++) {
        // Get the corresponding European option's value
        double euro_val = alpha * u_mesh[i - 1] + (1. - 2. * alpha) * u_mesh[i] + alpha * u_mesh[i + 1];
        // Find early exercise
        double early_ex_premium = 0.;        
        if (x_mesh[i] < 0.) {
            early_ex_premium = K_ * std::exp(a_ * x_mesh[i] + b_ * tau) * (1. - std::exp(x_mesh[i]));
        }
        // Compare and add to mesh
        next_u_mesh[i] = std::max(euro_val, early_ex_premium);
    }
    
    // Right boundary
    next_u_mesh[x_mesh.size()] = boundary_x_r_(tau);
    
    return next_u_mesh;
}

std::pair<std::vector<double>, std::vector<double>> FiniteDifferencePricer::buildUMeshFromFiniteDifferenceAmericanPut(
    const double alpha, const std::vector<double>& x_mesh, std::size_t M, double d_tau)
{
    // u_mesh initialization
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);
    std::vector<double> u_mesh_tau_plus_1(x_mesh.size());
    // printVector(u_mesh);
    
    // moving forward from t_1 to t_{M-1}
    for (std::size_t i = 1; i < M; ++i)
    {
        double curr_tau = d_tau * i;
        u_mesh_tau_plus_1 = computeNextUMeshAmerican(curr_tau, alpha, x_mesh, u_mesh);
        std::swap(u_mesh, u_mesh_tau_plus_1);
        // printVector(u_mesh);
    }

    // moving forward to t_M
    u_mesh_tau_plus_1 = computeNextUMeshAmerican(tau_final_, alpha, x_mesh, u_mesh);
    // printVector(u_mesh_tau_plus_1);

    return std::make_pair(u_mesh, u_mesh_tau_plus_1);
}

double FiniteDifferencePricer::varianceReductionAmericanPut(double V_approx, std::size_t M)
{
    // Get corresponding European option value for BS and FD
    FiniteDifferencePricer FDPricer(S0_, K_, T_, sigma_, r_, q_);
    FDPricer.setHyperparameters(0.45, 4);
    EuropeanOption BSPricer(0., S0_, K_, T_, sigma_, r_, q_);
    
    double FDEuro = FDPricer.priceEuropeanPut(M).front();
    double BSEuro = BSPricer.Put();
    
    // Adjust approximation by the pointwise difference
    return V_approx + BSEuro - FDEuro;
}


// boundary conditions setters
void FiniteDifferencePricer::setUpEuropeanCallBoundaryCondition(const std::size_t M) {}

void FiniteDifferencePricer::setUpEuropeanPutBoundaryCondition(const std::size_t M)
{
    // the boundary conditions for u(x, τ)
    boundary_tau_0_ = [=](double x)->double 
    {
        return K_ * std::exp(a_ * x) * std::max(1. - std::exp(x), 0.);
    };
    
    boundary_x_l_ = [=](double tau)->double 
    {
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
    };
    
    boundary_x_r_ = [](double tau)->double 
    {
        return 0.;
    };
}

void FiniteDifferencePricer::setUpAmericanCallBoundaryCondition(const std::size_t M) {}

void FiniteDifferencePricer::setUpAmericanPutBoundaryCondition(const std::size_t M) 
{
    // 2. Boundary conditions
    boundary_tau_0_ = [=](double x)->double 
    {
        return K_ * std::exp(a_ * x) * std::max(1. - std::exp(x), 0.);
    };
    
    boundary_x_l_ = [=](double tau)->double 
    {
        return K_ * std::exp(a_ * x_l_ + b_ * tau) * (std::exp(-2. * r_ * tau / (sigma_ * sigma_)) - std::exp(x_l_ - 2. * q_ * tau / (sigma_ * sigma_)));
    };
    
    boundary_x_r_ = [](double tau)->double 
    {
        return 0.;
    };
}

// public

void FiniteDifferencePricer::setHyperparameters(const double alpha_temp, const double M)
{
    alpha_temp_ = alpha_temp;
    M_ = M;
    flag_alpha_temp_defined_ = true;
    flag_M_defined_ = true;
}

void FiniteDifferencePricer::validateHyperparameters() const
{
    if (!flag_alpha_temp_defined_ || !flag_M_defined_)
    {
        throw("must define alpha_temp and M before using");
    }
}

FiniteDifferencePricer::FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q)
    : S0_(S0)
    , K_(K)
    , T_(T)
    , sigma_(sigma)
    , r_(r)
    , q_(q)
{
    double sigma_square = std::pow(sigma_, 2);
    // 1. computational domain
    tau_final_ = T_ * std::pow(sigma_, 2) / 2.;
    x_l_ = std::log(S0_ / K_) + (r - q - sigma_square / 2.) * T_ - 3 * sigma_ * std::sqrt(T_);
    x_r_ = std::log(S0_ / K_) + (r - q - sigma_square / 2.) * T_ + 3 * sigma_ * std::sqrt(T_);

    a_ = (r_ - q_) / sigma_square - .5;
    b_ = std::pow(a_ + 1., 2) + 2. * q_ / sigma_square;
}


// Pricers
std::vector<double> FiniteDifferencePricer::priceEuropeanPut(const std::size_t M)
{   
    validateHyperparameters();

    std::vector<double> res;

    double d_tau;
    std::size_t N;
    double d_x;
    double alpha;

    // 1. computational domain
    std::tie(d_tau, N, d_x, alpha) = compute_domain_params(M);
    // 2. boundary conditions
    setUpEuropeanPutBoundaryCondition(M);
    // 3. finite difference scheme
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = buildMesh(N, d_x);

    // 4. pointwise converage
    std::vector<double> u_mesh_M_minus_1_th, u_mesh;
    std::tie(u_mesh_M_minus_1_th, u_mesh) = buildUMeshFromFiniteDifference(alpha, x_mesh, M, d_tau);
    // printVector(u_mesh);
    std::vector<double> V_approximations = approximateSecurityValue(x_mesh, u_mesh, d_x);
    std::copy(V_approximations.cbegin(), V_approximations.cend(), std::back_inserter(res));

    // 5. RMS error
    double error_RMS = computeRMSError(x_mesh, u_mesh);
    res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = computeGreeks(interval_i, x_mesh, u_mesh, u_mesh_M_minus_1_th, d_tau, V_approximations[0]);
    std::copy(Greeks.cbegin(), Greeks.cend(), std::back_inserter(res));

    return res;   
}

std::vector<double> FiniteDifferencePricer::priceAmericanPut(const std::size_t M)
{
    validateHyperparameters();

    std::vector<double> res;

    double d_tau;
    std::size_t N;
    double d_x;
    double alpha;

    // 1. Computational domain
    std::tie(d_tau, N, d_x, alpha) = compute_domain_params(M);
    // 2. boundary conditions
    setUpAmericanPutBoundaryCondition(M);    
    // 3. Finite difference scheme
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = buildMesh(N, d_x);
    
    // 4. pointwise converage
    std::vector<double> u_mesh_M_minus_1_th, u_mesh;
    std::tie(u_mesh_M_minus_1_th, u_mesh) = buildUMeshFromFiniteDifferenceAmericanPut(alpha, x_mesh, M, d_tau);
    // printVector(u_mesh);
    std::vector<double> V_approximations = approximateSecurityValue(x_mesh, u_mesh, d_x);
    std::copy(V_approximations.cbegin(), V_approximations.cend(), std::back_inserter(res));
    
    // 6. Greeks
    std::vector<double> Greeks = computeGreeks(interval_i, x_mesh, u_mesh, u_mesh_M_minus_1_th, d_tau, V_approximations[0]);
    std::copy(Greeks.cbegin(), Greeks.cend(), std::back_inserter(res));
    
    // 7. Variance reduction
    res.push_back(varianceReductionAmericanPut(V_approximations[0], M));
    
    return res;
}

std::vector<double> FiniteDifferencePricer::compute_Sopt(
    const double alpha, const std::vector<double> &x_mesh, std::size_t M, double d_tau)
{    
    // u_mesh initialization
    std::vector<double> u_mesh(x_mesh.size());
    std::transform(x_mesh.cbegin(), x_mesh.cend(), u_mesh.begin(), boundary_tau_0_);

    std::vector<double> Sopt;
    for (std::size_t i = 1; i <= M; i++) 
    {
        // Record maximum exercise position
        std::size_t Nopt = 0;
        double tau = d_tau * i;
        
        std::vector<double> next_u_mesh(x_mesh.size());
        
        // Left boundary
        next_u_mesh[0] = boundary_x_l_(tau);
        
        // Middle values
        for (std::size_t i = 1; i < x_mesh.size() - 1; i++) {
            // Get the corresponding European option's value
            double euro_val = alpha * u_mesh[i - 1] + (1. - 2. * alpha) * u_mesh[i] + alpha * u_mesh[i + 1];
            
            // Find early exercise
            double early_ex_premium = 0.;
            
            if (x_mesh[i] < 0.) {
                early_ex_premium = K_ * std::exp(a_ * x_mesh[i] + b_ * tau) * (1. - std::exp(x_mesh[i]));
            }
            
            // Compare and add to mesh
            if (early_ex_premium >= euro_val) 
            {
                next_u_mesh[i] = early_ex_premium;
                if (early_ex_premium > 0) 
                {
                    Nopt = i;
                }
            } 
            else 
            {
                next_u_mesh[i] = euro_val;
            }
        }
        double S_small = K_ * std::exp(x_mesh[Nopt]);
        double S_large = K_ * std::exp(x_mesh[Nopt + 1]);
        Sopt.emplace_back((S_small + S_large) / 2.);
        
        // Right boundary
        next_u_mesh[x_mesh.size()] = boundary_x_r_(tau);
        
        u_mesh = std::move(next_u_mesh);
    }
    return Sopt;
}

std::pair<std::vector<double>, std::vector<double>> FiniteDifferencePricer::priceAmericanPutEarlyExercise(std::size_t M) 
{
    validateHyperparameters();

    std::vector<double> res;

    double d_tau;
    std::size_t N;
    double d_x;
    double alpha;
        
    // 1. Computational domain
    std::tie(d_tau, N, d_x, alpha) = compute_domain_params(M);
    // 2. boundary conditions
    setUpAmericanPutBoundaryCondition(M);    
    // 3. Finite difference scheme
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = buildMesh(N, d_x);
    // 4. pointwise converage
    // Advance M times while recording early exercise position
    std::vector<double> Sopt = compute_Sopt(alpha, x_mesh, M, d_tau);
    
    // Get array of t
    std::vector<double> t_mesh(16);
    std::iota(t_mesh.begin(), t_mesh.end(), 1.);
    auto convert_to_t = [&](double m)->double 
    {
        return T_ - (2. * m * d_tau) / std::pow(sigma_, 2);
    };
    std::transform(t_mesh.begin(), t_mesh.end(), t_mesh.begin(), convert_to_t);
    
    printVector(x_mesh);
    return std::make_pair(t_mesh, Sopt);
}


std::vector<double> FiniteDifferencePricer::priceEuropeanPut_BE(const std::size_t M)
{
    validateHyperparameters();

    std::vector<double> res;

    double d_tau;
    std::size_t N;
    double d_x;
    double alpha;

    // 1. computational domain
    std::tie(d_tau, N, d_x, alpha) = compute_domain_params(M);
    // 2. boundary conditions
    setUpEuropeanPutBoundaryCondition(M);
    // 3. finite difference scheme
    std::vector<double> x_mesh;
    std::size_t interval_i;
    std::tie(x_mesh, interval_i) = buildMesh(N, d_x);

    // 4. pointwise converge
    std::vector<double> u_mesh_M_minus_1_th, u_mesh;
    std::tie(u_mesh_M_minus_1_th, u_mesh) = buildUMeshFromFiniteDifference_BE(alpha, x_mesh, M, d_tau);
    std::vector<double> V_approximations = approximateSecurityValue(x_mesh, u_mesh, d_x);
    std::copy(V_approximations.cbegin(), V_approximations.cend(), std::back_inserter(res));

    // 5. RMS error
    double error_RMS = computeRMSError(x_mesh, u_mesh);
    res.push_back(error_RMS);
    
    // 6. Greeks
    std::vector<double> Greeks = computeGreeks(interval_i, x_mesh, u_mesh, u_mesh_M_minus_1_th, d_tau, V_approximations[0]);
    std::copy(Greeks.cbegin(), Greeks.cend(), std::back_inserter(res));

    return res;   
}

#include"mracs.h"

#define NUMTEST_THETA   5
#define NUMTEST_R       5
#define R0              10
#define R1              110

double two_point_correlation_Landy_Szalay(double r, double* s, double* s0);
double Proj_Value(double xx, double yy, double zz, double* s);
std::vector<double> return_and_print_3pcf(std::vector<double>& vec_r, std::vector<double>& vec_theta, double* s, double nphi, double np);

int main(){
    read_parameter();
    std::vector<Particle> p = read_in_DM_3vector(DataDirec);

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, 1);
    std::uniform_real_distribution<double> v(0,SimBoxL);

    std::vector<Particle> p0;
    for(size_t i = 0; i < p.size(); ++i) 
        p0.push_back({v(e), v(e), v(e), 1.});

    std::vector<double> vec_r, vec_theta, vec_s;
    std::vector<double> xi_r, xi_s, Q_norm;
    for(int i = 0; i < NUMTEST_R; ++i) 
        vec_r.push_back(R0 + (R1 - R0) * static_cast<double>(i)/NUMTEST_R);
    for(int i = 0; i < NUMTEST_THETA; ++i) 
        vec_theta.push_back(M_PI/2 * static_cast<double>(i+1)/NUMTEST_THETA);
    for(auto r : vec_r)
        for(auto theta : vec_theta){
            vec_s.push_back(2*r*sin(theta));
        }
    for(auto i : vec_r) std::cout << i << ", "; std::cout << std::endl; 
    for(auto i : vec_theta) std::cout << i/(M_PI/2) << ", "; std::cout << std::endl; 

    auto s  = sfc(p);
    auto s0 = sfc(p0);
    force_kernel_type(0);
    for(auto r : vec_r)
        xi_r.push_back(two_point_correlation_Landy_Szalay(r,s,s0));
    for(auto r : vec_s)
        xi_s.push_back(two_point_correlation_Landy_Szalay(r,s,s0));
    for(auto a : xi_r)
        for(auto b : xi_s)
            Q_norm.push_back(a*a + 2*a*b);

    auto d = new double[GridVol];
    #pragma omp parallel for
    for(size_t i = 0; i < GridVol; ++i) d[i] = (s[i]-s0[i]) / s0[i];

    // force_kernel_type(1);
    // auto w = wfc(Radius,0);
    // auto c = convol3d(s,w);
    // auto c0= convol3d(s0,w);

    // auto ds = new double[GridVol];
    // #pragma omp parallel for
    // for(size_t i = 0; i < GridVol; ++i) ds[i] = (c[i]-c0[i]) / c0[i];

    auto zeta_SS = return_and_print_3pcf(vec_r,vec_theta,d,10,10);

    std::cout << "Q: " << std::endl;
    for(int i = 0; i < NUMTEST_R; ++i){
        for(int j = 0; j < NUMTEST_THETA; ++j)
        {
            std::cout << zeta_SS[i* NUMTEST_THETA + j]/Q_norm[i* NUMTEST_THETA + j] << ", ";
        }
        std::cout << "\n---\n";
    }

}


double Proj_Value(double xx, double yy, double zz, double* s)
{
    double sum{0};
    int xxc = floor(xx), xxf = SampRate * (xx - xxc);      
    int yyc = floor(yy), yyf = SampRate * (yy - yyc);      
    int zzc = floor(zz), zzf = SampRate * (zz - zzc);
    
    for(int i = 0; i < phiSupport; ++i)
        for(int j = 0; j < phiSupport; ++j)
            for(int k = 0; k < phiSupport; ++k)
            {
                sum += s[((xxc-i) & (GridLen-1)) * GridLen * GridLen + ((yyc-j) & (GridLen-1)) * GridLen + ((zzc-k) & (GridLen-1))]
                        * phi[xxf + i*SampRate] * phi[yyf + j*SampRate] * phi[zzf + k*SampRate];
            }
    return sum;
}

double two_point_correlation_Landy_Szalay(double r, double* s, double* s0){
    auto w  = wfc(r,0);
    auto c  = convol3d(s,w);
    auto c0 = convol3d(s0, w);
    double xi = (inner_product(s,c,GridVol)-2*inner_product(s,c0,GridVol))/inner_product(s0,c0,GridVol)+1;
    delete[] w;
    delete[] c;
    delete[] c0;
    return xi;
}

std::vector<double> return_and_print_3pcf(std::vector<double>& vec_r, std::vector<double>& vec_theta, double* d, double nphi, double np)
{
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, 1);
    std::uniform_real_distribution<double> v(0,SimBoxL);
    std::vector<double> result;

    const int N_Phi_Twins = nphi;
    std::vector<double> dr_x(N_Phi_Twins);
    std::vector<double> dr_y(N_Phi_Twins);

    const int P_Sample = np; 
    const int64_t N_Sample = P_Sample * GridLen * GridLen;
    const double norm = N_Sample * N_Phi_Twins;

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
    for(int i = 0; i < vec_r.size(); ++i)
    {
        for(int j = 0; j < vec_theta.size(); ++j)
        {
            double pcf {0};
            double dr_xy = vec_r[i] * GridLen / SimBoxL * sin(vec_theta[j]);
            double dr_z  = vec_r[i] * GridLen / SimBoxL * cos(vec_theta[j]);

            for(int i = 0; i < N_Phi_Twins; ++i){
                dr_x[i] = (dr_xy * sin(i * M_PI / N_Phi_Twins));
                dr_y[i] = (dr_xy * cos(i * M_PI / N_Phi_Twins));}

            #pragma omp parallel for reduction (+:pcf)
            for(int i = 0; i < GridLen; ++i){
                for(int j = 0; j < GridLen; ++j){
                    for(int k = 0; k < P_Sample; ++k)
                    {
                        double x0 = i + u(e), y0 = j + u(e), z0 = v(e), a0, a1, a2, sum{0};
                        for(int ii = 0; ii < N_Phi_Twins; ++ii)
                        {
                            a1 = Proj_Value(x0 + dr_x[ii], y0 + dr_y[ii] ,z0, d);
                            a2 = Proj_Value(x0 - dr_x[ii], y0 - dr_y[ii] ,z0, d);
                            sum += a1 * a2;
                        }
                        a0  = Proj_Value(x0, y0, z0 - dr_z, d);
                        pcf += a0 * sum;
                    }
                }
            }
            result.push_back(pcf/norm);
            std::cout << pcf/norm << ", ";
        }
        std::cout << "\n---\n";
    }
    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference 3pcf       = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;

    return result;
}


/*
#include"mracs.h"

int main()
{
    read_parameter();
    auto p = read_in_DM_3vector(DataDirec);
    const int RS = 1;   // scale of random particles 
    const int RL = pow(p.size()*RS, 1/3.);
    const double norm = pow(RL, 3)/p.size();
    auto p0= generate_random_particle(RL,SimBoxL,0);

    const double Rmin{5}, Rmax{100}, NumR{10};
    const double Thetamin{0}, Thetamax{M_PI/2}, NumThta{10};
    std::vector<double> vec_radii, vec_theta, vec_pcf;


    auto s = sfc(p);
    auto s0= sfc(p0);
    force_kernel_type(1);
    auto w = wfc(Radius,0);
    auto c = convol3d(s,w);
    auto c0= convol3d(s0,w);
    delete[] w;

    
}
*/
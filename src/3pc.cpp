#include"mracs.h"

#define NUMTEST_THETA   1
#define NUMTEST_R       1
#define R0              1.
#define R1              51.


double Proj_Value(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport);

int main(){
    read_parameter();
    std::vector<Particle> p = read_in_TNG_3vector(DataDirec);

    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0, 1);
    std::uniform_real_distribution<double> v(0,GridLen);

    std::vector<Particle> p0;
    for(size_t i = 0; i < p.size(); ++i) p0.push_back({v(e), v(e), v(e), 1.});
    auto s  = sfc(p);
    auto s0 = sfc(p0);
    auto w  = wfc(Radius,0);
    auto c  = convol3d(s,w);
    auto c0 = convol3d(s0, w);

    std::vector<double> v_radii;
    std::vector<double> v_theta;
    for(int i = 0; i < NUMTEST_R; ++i) v_radii.push_back(R0 + (R1 - R0) * static_cast<double>(i)/NUMTEST_R);
    for(int i = 0; i < NUMTEST_THETA; ++i) v_theta.push_back(M_PI * static_cast<double>(i)/NUMTEST_THETA);
    for(auto i : v_radii) std::cout << i << ", "; std::cout << std::endl; 
    for(auto i : v_theta) std::cout << i/M_PI << ", "; std::cout << std::endl; 

    const int phiSupport = phi[phi.size() - 1] - phi[phi.size() - 2];
    std::vector<int> step(phiSupport);
    for(int i = 0; i < phiSupport; ++i) step[i] = i * SampRate;

    const int N_Phi_Twins = 2;
    std::vector<double> dr_x(N_Phi_Twins);
    std::vector<double> dr_y(N_Phi_Twins);

    const int P_Sample = 1; 
    const uint64_t N_Sample = P_Sample * GridLen * GridLen;
    const double ratio_v = 4/3.*M_PI*pow(Radius*GridLen/SimBoxL, 3)/GridNum;
    const double norm = pow(p.size()*ratio_v, 3)* N_Sample * N_Phi_Twins;
    

    std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
    for(int i = 0; i < v_radii.size(); ++i){
        for(int j = 0; j < v_theta.size(); ++j){
            double pcf {0};
            double dr_xy = v_radii[i] * GridLen / SimBoxL * sin(v_theta[j]);
            double dr_z  = v_radii[i] * GridLen / SimBoxL * cos(v_theta[j]);

            for(int i = 0; i < N_Phi_Twins; ++i){
                dr_x[i] = (dr_xy * sin(i * M_PI / N_Phi_Twins));
                dr_y[i] = (dr_xy * cos(i * M_PI / N_Phi_Twins));}

            #pragma omp parallel for reduction (+:pcf)
            for(int i = 0; i < GridLen; ++i){
                for(int j = 0; j < GridLen; ++j){
                    for(int k = 0; k < P_Sample; ++k){
                        double x0 = i + u(e), y0 = j + u(e), z0 = v(e), a0, a1, a2, a00, sum{0};
                        for(int ii = 0; ii < N_Phi_Twins; ++ii){
                            a1 = Proj_Value(x0 + dr_x[ii], y0 + dr_y[ii] ,z0, c, step, phiSupport);
                            a2 = Proj_Value(x0 - dr_x[ii], y0 - dr_y[ii] ,z0, c, step, phiSupport);
                            sum += a1 * a2;}
                        a0  = Proj_Value(x0, y0, z0 - dr_z, c , step, phiSupport);
                        a00 = Proj_Value(x0, y0, z0 - dr_z, c0, step, phiSupport);
                        pcf += (a0 - a00) * sum;}
                }
            }
            std::cout << pcf/norm + 2 << ", ";
        }
        std::cout << "\n---\n";
    }
    std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    std::cout << "Time difference 3pcf       = "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin1).count()
    << "[ms]" << std::endl;

}


double Proj_Value(double xx, double yy, double zz, double* s, std::vector<int> step, int phiSupport)
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
                        * phi[xxf + step[i]] * phi[yyf + step[j]] * phi[zzf + step[k]];
            }
    return sum;
}









/*
double two_point_correlation(double* s, fftw_complex* sc, int64_t n_particle, double r)
{
    if(!r) r = 0.001;
    auto w = wfc(r,0);
    auto c = convol_c2r(sc, w);
    double xi = inner_product(s,c,GridNum) * GridNum/pow(n_particle, 2) - 1;
    delete[] w;
    delete[] c;
    return xi;
}

    // force_kernel_type(0);
    // std::chrono::steady_clock::time_point begin0 = std::chrono::steady_clock::now();
    // for(auto i : v_radii){
    //     for(auto j : v_theta)
    //         Xi.push_back(two_point_correlation(s,sc,p.size(),2*i*sin(j)));
    //     Xi_r.push_back(two_point_correlation(s,sc,p.size(),i));
    // }
    // std::chrono::steady_clock::time_point end0 = std::chrono::steady_clock::now();


    auto p = read_in_TNG_3vector(DataDirec);
    for(int i = 0; i < 30; ++i) cout << p[i].x << " " << p[i].y << " " << p[i].z << " " << p[i].weight << endl;
    //for(auto i : p) cout << i.x << " " << i.y << " " << i.z << " " << i.weight << endl;
    vector<Particle>().swap(p);
int main(){
    read_parameter();
    auto g = read_in_Millennium_Run_galaxy_catalog(DataDirec);
    cout << g.size() << endl;
    std::cout << "before:" << std::endl;
    for(int i = 0; i < 10; ++i){
        std::cout << g[i].x << ", " << g[i].y << ", " << g[i].z << ", " << g[i].vx << ", " << g[i].vy << ", " << g[i].vz << ", " << g[i].BulgeMass <<"\n";
    }
    float temp;
    for(size_t i = 0; i < g.size(); ++i){
        temp = g[i].z + g[i].vz/100;
        g[i].z = temp - floor(temp/SimBoxL)*SimBoxL;
    }
    std::cout<< "after: " << std::endl;
    for(int i = 0; i < 10; ++i){
        std::cout << g[i].x << ", " << g[i].y << ", " << g[i].z << ", " << g[i].vx << ", " << g[i].vy << ", " << g[i].vz << ", " << g[i].BulgeMass <<"\n";
    }
    std::vector<float> d;
    //x >> y >> z >> vx >> vy >> vz >> Mag_u >> Mag_g >> Mag_r >> Mag_i >> Mag_z >>
            //BulgeMag_u >> BulgeMag_g >> BulgeMag_r >> BulgeMag_i >> BulgeMag_z >> StellarMass >>
            //BulgeMass >> ColdGas >> HotGas >> EjectedMass >> BlackHoleMass >> Sfr
    //for(auto a : g) {
    //    d.push_back(a.x);d.push_back(a.y);d.push_back(a.z);d.push_back(a.vx);d.push_back(a.vy);d.push_back(a.vz);d.push_back(a.Mag_u);d.push_back(a.Mag_g);
    //    d.push_back(a.Mag_r);d.push_back(a.Mag_i);d.push_back(a.Mag_z);d.push_back(a.BulgeMag_u);d.push_back(a.BulgeMag_g);d.push_back(a.BulgeMag_r);d.push_back(a.BulgeMag_i);d.push_back(a.BulgeMag_z);
    //    d.push_back(a.StellarMass);d.push_back(a.BulgeMass);d.push_back(a.ColdGas);d.push_back(a.HotGas);d.push_back(a.EjectedMass);d.push_back(a.BlackHoleMass);d.push_back(a.Sfr);
    //    }
    /*
    std::string ofname = "../simdata/croton_etal.ugriz.rsd.bin";
    std::ofstream ofs(ofname, std::ios_base::binary);
    if(!ofs){
        std::cout << "write open error" << std::endl;
    }
    unsigned int total = g.size();
    ofs.write((char*) &total +3, sizeof(char));
    ofs.write((char*) &total +2, sizeof(char));
    ofs.write((char*) &total +1, sizeof(char));
    ofs.write((char*) &total , sizeof(char));

    void* addr = &g[0];
    for(int i = 0; i < 23; ++i){
        for(int n = 0; n < g.size(); ++n)
            for(int j = 3; j >= 0; --j)
                ofs.write(((char*) addr) + (i + n*23)*4 +j , sizeof(char));
    }

}*/


/*
int main(){
    read_parameter();
    auto p = read_in_TNG_3vector(DataDirec);
    std::cout << p[0].x << ", " << p[0].y << ", " << p[0].z << std::endl;
    return 0;
}*/

/*
std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for(auto dr : v_radii){
        for(auto theta : v_theta)
        {
            double pcf {0};
            double dr_xy = dr * sin(theta);
            double dr_z  = dr * cos(theta);

            for(int i = 0; i < N_Phi_Twins; ++i){
                dr_x[i] = (dr_xy * sin(i * M_PI / N_Phi_Twins));
                dr_y[i] = (dr_xy * cos(i * M_PI / N_Phi_Twins));
            }

            #pragma omp parallel for reduction (+:pcf)
            for(int i = 0; i < GridLen; ++i){
                for(int j = 0; j < GridLen; ++j)
                {
                    double sum {0};
                    double x0,y0,z0;
                    double xx,yy,zz;
                    double a0,a1,a2;

                    x0 = i + u(e);
                    y0 = j + u(e);
                    z0 = v(e);

                    for(int Ti = 0; Ti < N_Phi_Twins; ++Ti){
                        xx = x0 + dr_x[Ti];
                        yy = y0 + dr_y[Ti];
                        zz = z0;
                        a1 = Proj_Value(xx,yy,zz,c,step,phiSupport);

                        xx = x0 - dr_x[Ti];
                        yy = y0 - dr_y[Ti];
                        zz = z0;
                        a2 = Proj_Value(xx,yy,zz,c,step,phiSupport);

                        sum += a1 * a2;
                    }
                    xx = x0;
                    yy = y0;
                    zz = z0 - dr_z;
                    a0 = Proj_Value(xx,yy,zz,c,step,phiSupport);

                    pcf += a0 * sum;
                }
            }

            pcf /= N_Sample * N_Phi_Twins * pow(p.size(),3) / pow(GridNum, 3) ;
            std::cout << pcf << ", ";
        }
        std::cout << "\n---\n";
    }
*/
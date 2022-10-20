#include"mracs.h"

using namespace std;


int main(){
    read_parameter();
    // std::vector<Galaxy> g = read_in_Millennium_Run_galaxy_catalog(DataDirec);
    // std::vector<Particle> p; for(Galaxy i : g) p.push_back({i.x, i.y, i.z, 1.});
    // std::vector<Galaxy>().swap(g);
    // auto p = read_in_TNG_3vector(DataDirec);
    //auto s = sfc(p);
    //
    //double* Pk_grid     = densityPowerGridFrame(s);
    //double* Pk_gridFree = densityPowerDWT(s);

    for(int i = 0; i < GridLen/2+1; ++i) std::cout << i/75.             << ", "; std::cout << '\n' << '\n';
    // for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_grid[i]            << ", "; std::cout << '\n' << '\n';
    // for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_gridFree[i]        << ", "; std::cout << '\n' << '\n';
    // for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_grid[i]*i*i*i      << ", "; std::cout << '\n' << '\n';
    // for(int i = 0; i < GridLen/2+1; ++i) std::cout << Pk_gridFree[i]*i*i*i  << ", "; std::cout << '\n';
}
/*
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
#include"mracs.h"
#include<list>

struct tuple
{
    int i;
    double x;
};
int main(){
    read_parameter();
    const int Mbin {4}; // default Mass split bin
    const int Nbin {4}; // else parameters split bin
    const int Npr {79}; // number of parameters
    std::list<int> exclid {0,1,2,3,4,5,6,7,8,14,31,32,35,36,57,63,69,70,71,77}; //exclusion index
    std::list<int> Forthin {0,1,2}; // will be substitude to Envir, LocDensity, Concentration
    
    std::cout << "reading halo catalog..\n";
    std::string ifn {"/data0/MDPL2/groups_130/cata_Mcut2e12.dat"};
    
    double a;
    std::vector<double> cata;
    std::ifstream ifs {ifn};
    while(ifs >> a) cata.push_back(a);

    //auto cata = read_in_float(ifn);
    
    if(cata.size()%Npr != 0){
        std::cout << "input error,abort\n";
        std::terminate();
    }
    auto dm = read_in_DM_3vector("/data0/MDPL2/dm_sub/dm_sub5e-4.bin");
    
    std::string ifnEnv {"output/envi_J10_GSR3_halo_Mcut2e12.txt"};
    int Npart = cata.size()/Npr;
    auto envi = envi_vector_readin(ifnEnv,Npart);

    std::vector<Particle> p(Npart);
    for(size_t i = 0; i < Npart; ++i){
        p[i] = {cata[i*Npr + 17], cata[i*Npr + 18], cata[i*Npr + 19], cata[i*Npr + 10]}; // x,y,z,Mvir
    }
    auto w = wfc(Radius,0);
    auto c = convol3d(sfc(p),w,true);
    delete[] w;

    auto locDens = project_value(c,p,true);

    
    std::vector<double> EnvirNode {0.5,1.5,2.5};
    std::vector<double> Mdata;
    for(auto x : p) Mdata.push_back(x.weight);
    auto Mnode = nodes_of_proto_sort(Mdata,Mbin);
    std::vector<tuple> ccc;
    
    for(int i = 0; i < Npr; ++i)
    {
        if(std::find(exclid.begin(),exclid.end(),i) == exclid.end() || std::find(Forthin.begin(),Forthin.end(),i) != Forthin.end()){
            std::cout << "param: ";
            std::vector<double> data(Npart);
            if(i == 0)      {std::cout << "Evi\n";   for(size_t j = 0; j < Npart; ++j) data[j] = envi[j];}
            else if(i == 1) {std::cout << "Loc\n";   for(size_t j = 0; j < Npart; ++j) data[j] = locDens[j];}
            else if(i == 2) {std::cout << "Con\n";   for(size_t j = 0; j < Npart; ++j) data[j] = cata[j*Npr + 11]/cata[j*Npr + 12];}
            else            {std::cout << i << "\n"; for(size_t j = 0; j < Npart; ++j) data[j] = cata[j*Npr + i];}
            
            std::vector<double> node;
            if(i == 0) node = EnvirNode;
            else       node = nodes_of_proto_sort(data,Nbin);

            std::vector<std::vector<Particle>*> vpts;

            if(i == 10) {
                for(int i = 0; i < Mbin; ++i) vpts.push_back(new std::vector<Particle>);
                for(size_t j = 0; j < Npart; ++j) 
                    vpts[classify_index(Mnode,Mdata[j])]->push_back({p[j]}); 
            }
            else {
                for(int i = 0; i < Mbin * Nbin; ++i) vpts.push_back(new std::vector<Particle>);
                for(size_t j = 0; j < Npart; ++j)
                    vpts[classify_index(Mnode,Mdata[j])*Nbin + classify_index(node,data[j])]->push_back({p[j]}); 
            }
            
            auto solve = optimal_solution(dm,vpts,Radius,true);

            ccc.push_back({i,solve[0]});

        }

    }
    //std::cout << "orignal:\n";
    //for(auto x : ccc) std::cout << x.i << ", "; std::cout << "\n";
    //for(auto x : ccc) std::cout << x.x << ", "; std::cout << "\n";

    // -----sorting and print-----
    std::vector<tuple> sorted;
    for(int i = 0; i < ccc.size(); ++i){
        double max = ccc[i].x;
        int index {i};
        for(int j = i; j < ccc.size(); ++j){
            if(ccc[j].x > max) {
                max = ccc[j].x;
                index = j;
            }
        }
        tuple tmp = {ccc[i].i,ccc[i].x};
        ccc[i] = {ccc[index].i,ccc[index].x};
        ccc[index] = tmp;
        
    }

    std::cout << "Sorted: (Nparam= " << ccc.size() << ")\n";
    std::cout << "[id]: r_m\n";
    for(int i = 0; auto x : ccc) {
        if(x.i == 0)       std::cout << "{Evi}";
        else if (x.i == 1) std::cout << "{Loc}";
        else if (x.i == 2) std::cout << "{Con}";
        else             std::cout << "[" << std::setw(2) << x.i << "]:";
        std::cout << std::setw(9) << x.x << ",    ";
        if(i%5 == 4) std::cout << std::endl;
        ++i;
    }std::cout << "\n";
   

}



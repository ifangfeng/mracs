#include"mracs.h"
#include<list>

struct tuple
{
    int i;
    double x;
};
int main(){
    read_parameter();
    const int Nbin {4}; // reconstruct bin
    const int Npr {79}; // number of parameters
    std::list<int> exid {0,2,3,4,5,6,7,8,14,31,32,35,36,57,63,69,70,71,77}; //exclusion index
    
    
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
    int Npart = cata.size()/Npr;
    std::vector<Particle> p(Npart);
    for(size_t i = 0; i < Npart; ++i){
        p[i] = {cata[i*Npr + 17], cata[i*Npr + 18], cata[i*Npr + 19], cata[i*Npr + 10]}; // x,y,z,Mvir
    }

    std::vector<tuple> ccc;
    
    for(int i = 0; i < Npr; ++i)
    {
        if(std::find(exid.begin(),exid.end(),i) == exid.end()){
            std::cout << "param: " << i << std::endl;
            std::vector<double> data(Npart);
            for(size_t j = 0; j < Npart; ++j){
                data[j] = cata[j*Npr + i];
            }
            auto node = nodes_of_proto_sort(data,Nbin);

            std::vector<std::vector<Particle>*> vpts;
            for(int i = 0; i < Nbin; ++i) vpts.push_back(new std::vector<Particle>);

            for(size_t j = 0; j < Npart; ++j){
                vpts[classify_index(node,data[j])]->push_back({p[j]}); 
            }
            auto solve = optimal_solution(dm,vpts,Radius,true);

            ccc.push_back({i,solve[0]});

        }

    }
    std::cout << "orignal:\n";
    for(auto x : ccc) std::cout << x.i << ", "; std::cout << "\n";
    for(auto x : ccc) std::cout << x.x << ", "; std::cout << "\n";

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
        std::cout << "[" << std::setw(2) << x.i << "]:" << std::setw(9) << x.x << ",    ";
        if(i%5 == 4) std::cout << std::endl;
        ++i;
    }
   

}



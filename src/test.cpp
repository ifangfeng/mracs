#include"mracs.h"
void trimming(std::vector<std::vector<double>*>& vpts){
    // -------size check before eigen solver-----
    bool EmptySize {false};
    const int ThdSize {1};

    for(auto x : vpts) 
        if (x->size() < ThdSize) 
            EmptySize = true;

    if(EmptySize){
        std::cout << "!EmptySize, some elements will be removed\n";
        std::cout << "---+bf\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
        for(int i = 0; i < vpts.size(); ++i) {
            if(vpts[i]->size() < ThdSize){
                delete vpts[i];
                vpts.erase(vpts.begin()+i);
                --i;
            }
        }
        std::cout << "---+af\n" << "size: ";
        std::cout << vpts.size() << "\n";for(auto x : vpts) std::cout << x->size() << ", "; std::cout << "\n";
    }
}
int main(){
    std::vector<std::vector<double>*> vpts;
    for(int i = 0; i < 4; ++i) vpts.push_back(new std::vector<double>);
    vpts[0]->push_back({1});
    vpts[2]->push_back({2});
    vpts[3]->push_back({3});

    trimming(vpts);
    std::cout <<"after trimming: "<< vpts.size() << ":\n";
    for(auto x : vpts) std::cout << x->size() << ", ";std::cout << std::endl;
}
#include<iostream>
#include<vector>

struct Bar {
    double vol, cur;
};

struct Node {
    double potential, current;
    int in_num = 0;
    int out_num = 0;
    std::vector<Bar> ins = {};
    std::vector<Bar> outs = {};
    
    void connect(Bar element, char side){
        if (side == 'o'){
            this->outs.push_back(element);
            ++out_num;
        } else {
            this->ins.push_back(element);
            ++in_num;
        }
    }

    void print(){
        std::cout << "This is a node with " << in_num << " ins and " << out_num << " outs.";
    }
};

int main() {
    Node a;
    Bar r1, r2, r3, r4;
    a.connect(r1, 'i');
    a.connect(r2, 'i');
    a.connect(r3, 'i');
    a.connect(r4, 'o');
    a.print();
    return 0;
}
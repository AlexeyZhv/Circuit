#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

struct Node;

struct Bar {
    private:
        Node* startNode;
        Node* endNode;
        std::string name;

    public:

        Bar(std::string name_) {
            startNode = nullptr;
            endNode = nullptr;
            name = name_;
        }

        Node* get_nodes() {
            return startNode, endNode;
        }

        Bar(Node* start, Node* end) {
            startNode = start;
            endNode = end;
        }

        void print_name() {
            std::cout << name << std::endl;
        }

        void connect_start(Node* node) {
            this->startNode = node;
        }

        void connect_end(Node* node) {
            this->endNode = node;
        }
};

struct Node {
    private:
        std::vector<Bar*> ins;
        std::vector<Bar*> outs;
        unsigned int in_num = 0;
        unsigned int out_num = 0;

    public:
        Node() {}

        void connect(Bar *&element, char side){
            if (side == 'o'){
                this->outs.push_back(element);
                element->connect_start(this);
                ++out_num;
            } else {
                this->ins.push_back(element);
                element->connect_end(this);
                ++in_num;
            }
        }

        void disconnect(Bar *&element){
            auto it_ins = std::find(ins.begin(), ins.end(), element);
            auto it_outs = std::find(outs.begin(), outs.end(), element);

            assert ((it_ins != ins.end() || it_outs != outs.end()) && "No such element connected to this node!");

            if (it_ins != ins.end()) {
                ins.erase(it_ins);
                element->connect_end(nullptr);
                --in_num;
            } else {
                outs.erase(it_outs);
                element->connect_start(nullptr);
                --out_num;
            }
        }

        void print(){
            std::cout << "This is a node with " << in_num << " ins and " << out_num << " outs.\n" << "ins:" << std::endl;

            for (size_t i = 0; i != ins.size(); ++i) {
                ins[i]->print_name();
            }
        }
};



int main() {
    Node a;
    Bar* r1 = new Bar("resistor");
    Bar* r2 = new Bar("capacitor");
    Bar* r3 = new Bar("battery");
    Bar* r4 = new Bar("diode");

    a.connect(r1, 'i');
    a.connect(r2, 'i');
    a.connect(r3, 'i');
    a.connect(r4, 'o');
    a.print();
    a.disconnect(r2);
    std::cout << std::endl;
    a.print();

    delete r1, r2, r3, r4;
    return 0;
}
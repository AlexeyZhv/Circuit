#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

struct Node;

struct Bar {
    private:
        double cur, vol;
        Node* startNode;
        Node* endNode;
        std::string name;

    public:

        Bar(std::string name_) {
            startNode = nullptr;
            endNode = nullptr;
            name = name_;
        }

        Node* get_startNode() {
            return startNode;
        }

        Node* get_endNode() {
            return endNode;
        }

        Bar(Node* start, Node* end) {
            startNode = start;
            endNode = end;
        }

        void set_vol(double vol_) {
            vol = vol_;

        }

        void set_cur(double cur_) {
            cur = cur_;
        }

        double get_vol() {
            return vol;
        }

        double get_cur() {
            return cur;
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

        virtual double vac(double voltage) { //returns current depending on voltage
            return 0;
        }
};

struct Node {
    private:
        std::vector<Bar*> ins;
        std::vector<Bar*> outs;
        unsigned int in_num = 0;
        unsigned int out_num = 0;
        double potential;

    public:
        void set_potential(double potential_){
            potential = potential_;
        }

        double get_potential() {
            return potential; 
        }

        std::vector<Bar*> get_ins() {
            return ins;
        }

        std::vector<Bar*> get_outs() {
            return outs;
        }

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

struct Resistor : Bar
{
    private:
        double res;

    public:
        Resistor(double res_, std::string name) : Bar(name) {
            res = res_;
        }

        void vac() {
            this->set_vol(this->get_startNode()->get_potential() - this->get_endNode()->get_potential());
            this->set_cur(this->get_vol() / this->res);
        }

        double get_res() {
            return res;
        }
};

struct Circuit {
    private:

        std::vector<Node*> nodes;
        std::vector<Bar*> bars;
    
    public:

        void print() {
            std::cout << "This is a circuit with following elements" << std::endl;
            for (size_t i = 0; i != bars.size(); ++i) {
                bars[i]->print_name();
            }
        }

        void add_node(Node*& node) {
            nodes.push_back(node);
        }

        void add_bar(Bar*& bar) {
            bars.push_back(bar);
        }

        void create () {
            
        }
};


int main() {
    Node* a = new Node;
    Bar* r1 = new Bar("resistor");
    Bar* r2 = new Bar("capacitor");
    Bar* r3 = new Bar("battery");
    Bar* r4 = new Bar("diode");

    a->connect(r1, 'i');
    a->connect(r2, 'i');
    a->connect(r3, 'i');
    a->connect(r4, 'o');
    a->print();
    a->disconnect(r2);
    std::cout << std::endl;
    a->print();

    Circuit c1;

    c1.add_node(a);
    c1.add_node(a);
    c1.add_node(a);
    c1.add_node(a);

    c1.print();

    delete r1, r2, r3, r4, a;
    return 0;
}
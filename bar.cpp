#include <vector>
#include <iostream>
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

    virtual void vac(double time, double step) { //returns current depending on voltage
    }

    virtual char get_type() {
        return 'N';
    }

    virtual double get_res() {
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

        void connect(Bar *element, char side){
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
                ins.shrink_to_fit();
            } else {
                outs.erase(it_outs);
                element->connect_start(nullptr);
                --out_num;
                outs.shrink_to_fit();
            }
        }

        void print(){
            std::cout << "This is a node with " << in_num << " ins and " << out_num << " outs.\n" << "ins:" << std::endl;
            for (size_t i = 0; i != ins.size(); ++i) {
                ins[i]->print_name();
            }
            std::cout << "outs:" << std::endl;
            for (size_t i = 0; i != outs.size(); ++i) {
                outs[i]->print_name();
            }
        }
};
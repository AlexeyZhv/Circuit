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

        virtual char get_type() {
            return 'N';
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

        char get_type() {
            return 'R';
        }
};

struct Voltage_source : Bar
{
    public:
        Voltage_source(double vol, std::string name) : Bar(name) {
            this->set_vol(vol);
        }

        char get_type() {
            return 'V';
        }

        void vac() {}
};

struct Wire : Bar
{
    public:
        Wire(double vol, std::string name) : Bar(name) {
            this->set_vol(0);
        }

        char get_type() {
            return 'V';
        }

        void vac() {}
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


// a function that solves a system of equations via Gauss algorythm

std::vector<double> gauss_solve(std::vector<std::vector<double>> matrix) {
    size_t checksum = 0;
    size_t matrix_height = matrix.size();
    size_t matrix_width = matrix_height + 1;
    std::vector<double> answer;

    //checking if matrix dimensions are correct

    for (size_t i = 0; i != matrix_height; ++i) {
        if (matrix[i].size() != matrix_width) {++ checksum;};
    }
    assert (checksum == 0 && "Wrong matrix dimensions!");

    //executing Gauss algorythm

    //forwards

    for (size_t j = 0; j != matrix_height - 1; ++j) {

        double coef = matrix[j][j];

        // swapping two rows if coef is 0
        std::vector<double> tmp;
        size_t counter = j;

        while (coef == 0) {
            ++counter;
            assert (counter != matrix_height && "System has 0 or infinite number of solutuons");
            coef = matrix[counter][j];
        }

        tmp = matrix[j];
        matrix[j] = matrix[counter];
        matrix[counter] = tmp;
        

        for (size_t i = j + 1; i != matrix_height; ++i) {
            double mul = matrix[i][j] / coef;

                for (size_t k = j; k != matrix_width; ++k) {
                    matrix[i][k] -= mul * matrix[j][k];
                }          
        }
    }

    //backwards

    for (size_t j = matrix_height - 1; j != 0; --j) {
        double coef = matrix[j][j];

        for (size_t i = j - 1; i != -1; --i) {
            double mul = matrix[i][j] / coef;

                for (size_t k = j; k != matrix_width; ++k) {
                    matrix[i][k] -= mul * matrix[j][k];
                }          
        }
    }

    for (size_t j = 0; j != matrix_height; ++j) {
        answer.push_back(matrix[j][matrix_width - 1] / matrix[j][j]);
    }  
    
    return answer;
}


// this is a test function for gauss_solve
void gauss_solve_test() {
    std::vector<std::vector<double>> matrix;
    std::vector<double> line1;
    std::vector<double> line2;
    std::vector<double> line3;

    std::vector<double> ans;

    line1.push_back(1); line1.push_back(2); line1.push_back(3); line1.push_back(4);
    line2.push_back(2); line2.push_back(4); line2.push_back(5); line2.push_back(10); 
    line3.push_back(3); line3.push_back(7); line3.push_back(4); line3.push_back(15);

    matrix.push_back(line1); 
    matrix.push_back(line2); 
    matrix.push_back(line3); 

    ans = gauss_solve(matrix);

    for (size_t i = 0; i != ans.size(); ++i) {
        std::cout << ans[i] << std::endl;
    }
}


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


    delete r1, r2, r3, r4, a;

    gauss_solve_test();
    return 0;
}
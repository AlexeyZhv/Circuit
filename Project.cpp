#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <unordered_map>

# define PI 3.1415926535897932384626433832795

// deletes all elements from a vector

template<typename T>
void delete_vector_elements(std::vector<T>& vec) {
    for (auto& elem : vec) {
        delete elem;
    }
    vec.clear();
}

void print_matrix(std::vector<std::vector<double>> matrix) {
    for (int i = 0; i != matrix.size(); ++i) {
        for (int j = 0; j != matrix[i].size(); ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

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

    //executing Gauss algorithm

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

        assert (coef != 0 && "System has 0 or infinite number of solutuons");

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


struct Node;

struct Bar {
private:
    double cur = 0;
    double vol = 0;
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

    virtual void print_name() {
        std::cout << name << std::endl;
    }

    void connect_start(Node* node) {
        this->startNode = node;
    }

    void connect_end(Node* node) {
        this->endNode = node;
    }

    virtual void vac(double time, double step) {}

    virtual char get_type() {
        return 'B';
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

struct Resistor : Bar
{
    private:
        double res;

    public:
        Resistor(double res_) : Bar("resistor") {
            res = res_;
        }

        Resistor(double res_, Node* start, Node* end) : Bar("resistor") {
            res = res_;
            start->connect(this, 'o'); 
            end->connect(this, 'i');
        }

        void vac(double time, double step) {
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

struct Current_source : Bar
{
    public:
        Current_source(double cur) : Bar("current source") {
            this->set_cur(cur);
        }

        Current_source(double cur, Node* start, Node* end) : Bar("current source") {
            this->set_cur(cur);
            start->connect(this, 'o');
            end->connect(this, 'i');
        }

        char get_type() {
            return 'I';
        }

        void vac(double time, double step) {
            this->set_vol(this->get_startNode()->get_potential() - this->get_endNode()->get_potential());
        }
};


struct Voltage_source : Bar
{
    public:
        Voltage_source(double vol) : Bar("voltage source") {
            this->set_vol(vol);
        }

        Voltage_source(double vol, Node* start, Node* end) : Bar("voltage source") {
            this->set_vol(vol);
            start->connect(this, 'o');
            end->connect(this, 'i');
        }


        char get_type() {
            return 'V';
        }

        void vac(double time, double step) {}
};


struct Sin_voltage_source : Voltage_source
{   
    private:
        double phase, phase_0, freq, ampl;
    public:
        Sin_voltage_source(double ampl_, double freq_, double phase_0_) : Voltage_source(ampl_ * sin(phase_0_)) {
            phase = phase_0_;
            phase_0 = phase_0_;
            ampl = ampl_;
            freq = freq_;
        }

        Sin_voltage_source(double ampl_, double freq_, double phase_0_, Node* start, Node* end) : Voltage_source(ampl_ * sin(phase_0_), start, end) {
            phase = phase_0_;
            phase_0 = phase_0_;
            ampl = ampl_;
            freq = freq_;
        }

        char get_type() {
            return 'V';
        }

        void vac(double time, double step) {
           this->set_vol(ampl * sin(phase));
           phase = phase_0 + 2 * PI * time * freq; 
        }
};


struct Wire : Voltage_source
{
    public:
        Wire() : Voltage_source(0) {
        }

        Wire(Node* start, Node* end) : Voltage_source(0, start, end) {
            this->set_vol(0);
        }


};


struct Capacitor : Voltage_source
{
    private:
        double cap, charge, charge_0;

    public:
        Capacitor(double cap_, double charge_0_) : Voltage_source(0) {
            cap = cap_;
            charge_0 = charge_0_;
            charge = charge_0;
        }

        Capacitor(double cap_, double charge_0_, Node* start, Node* end) : Voltage_source(0, start, end) {
            cap = cap_;
            charge_0 = charge_0_;
            charge = charge_0;
        }

        void vac(double time, double step) {
            charge -= this->get_cur() * step;
            this->set_vol(charge / cap);
        }

            void print_name() {
        std::cout << "capacitor" << std::endl;
        }   

};

struct Inductor : Current_source
{
    private:
        double ind, cur_0;
    public:
        Inductor(double ind_, double cur_0_) : Current_source(cur_0_) {
            ind = ind_;
            cur_0 = cur_0_;
            this->set_cur(cur_0);
        }

        Inductor(double ind_, double cur_0_, Node* start, Node* end) : Current_source(cur_0_, start, end) {
            ind = ind_;
            cur_0 = cur_0_;
            this->set_cur(cur_0);
        }

        char get_type() {
            return 'I';
        }

        void vac(double time, double step) {
            double cur_ = this->get_cur();
            this->set_vol(this->get_startNode()->get_potential() - this->get_endNode()->get_potential());
            this->set_cur(cur_ + this->get_vol() * step / ind);
        }

        void print_name() {
            std::cout << "inductor" << std::endl;
        }

};


// Oscilloscopes

struct Voltmeter : Current_source 
{
    public:
        Voltmeter() : Current_source(0) {}

        Voltmeter(Node* start, Node* end) : Current_source(0, start, end) {}

        void vac() {
            std::cout << this->get_vol() << std::endl; // FIXME Здесь нужен вывод в файл
        }
};


struct Ampermeter : Wire
{
    public:
        Ampermeter() : Wire() {}

        Ampermeter(Node* start, Node* end) : Wire(start, end) {}

        void vac() {
            std::cout << this->get_cur() << std::endl; // FIXME Здесь нужен вывод в файл
        }
};


struct Circuit {
    private:

        std::vector<Node*> nodes;
        std::vector<Bar*> resistors;
        std::vector<Bar*> cur_sources;
        std::vector<Bar*> wires;
    
    public:

        void print() {
            std::cout << "This is a circuit with following elements" << std::endl;
            for (size_t i = 0; i != resistors.size(); ++i) {
                resistors[i]->print_name();
            }
            for (size_t i = 0; i != wires.size(); ++i) {
                wires[i]->print_name();
            }
            for (size_t i = 0; i != cur_sources.size(); ++i) {
                cur_sources[i]->print_name();
            }
        }

        void add_node(Node*& node) {
            nodes.push_back(node);
        }

        void add_bar(Bar*& bar) {
            if (bar->get_type() == 'R') {
                resistors.push_back(bar);
            };

            if (bar->get_type() == 'V') {
                wires.push_back(bar);
            }

            if (bar->get_type() == 'I') {
                cur_sources.push_back(bar);
            }
        }

        void solve(double time, double step) {
            size_t height = nodes.size() + wires.size();
            size_t width = nodes.size() + wires.size() + 1;

            std::vector<std::vector<double>> matrix;
            std::vector<double> tmp;

            // setting first node's potential to zero to avoid uncertainty

            tmp.push_back(1);

            for (size_t i = 1; i != width; ++i) {
                tmp.push_back(0);
            }

            matrix.push_back(tmp);

            tmp.assign(width, 0);

            // adding equations for node currents to matrix

            for (size_t i = 1; i != nodes.size(); ++i) {
                std::vector<Bar*> ins_ = nodes[i]->get_ins();
                std::vector<Bar*> outs_ = nodes[i]->get_outs();

                // bars pointing inside the node

                for (size_t j = 0; j < ins_.size(); ++j) {
                    
                    // voltage sources and wires

                    if (ins_[j]->get_type() == 'V') {
                        size_t idx = nodes.size() + std::distance(wires.begin(), std::find(wires.begin(), wires.end(), ins_[j]));

                        tmp[idx] += 1;
                    } 

                    else{
                    // resistors

                    if (ins_[j]->get_type() == 'R') {

                        size_t idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), ins_[j]->get_startNode()));

                        assert (idx != nodes.size() && "no such node in this circuit!");

                        tmp[idx] += 1 / ins_[j]->get_res();
                        tmp[i] -= 1 / ins_[j]->get_res();
                    } 
                    else{
                    // current_sources

                    if (ins_[j]->get_type() == 'I') {
                        
                        // adding to total current of this node

                        tmp[width - 1] -= ins_[j]->get_cur(); // this element
                    } 
                    }
                    }
                }

                // bars pointing outside the node

                for (size_t j = 0; j < outs_.size(); ++j) {

                    // voltage sources and wires

                    if (outs_[j]->get_type() == 'V') {
                        size_t idx = nodes.size() + std::distance(wires.begin(), std::find(wires.begin(), wires.end(), outs_[j]));

                        tmp[idx] -= 1;
                    } 
                    else{
                    // resistors

                    if (outs_[j]->get_type() == 'R') {
                        size_t idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), outs_[j]->get_endNode()));

                        assert (idx != nodes.size() && "no such node in this circuit!");

                        tmp[idx] += 1 / outs_[j]->get_res();
                        tmp[i] -= 1 / outs_[j]->get_res();
                    }
                    else{
                    // current_sources

                    if (outs_[j]->get_type() == 'I') {

                        // adding to total current of this node

                        tmp[width - 1] += outs_[j]->get_cur(); // this element  
                    }
                    }  
                    }
                }
                matrix.push_back(tmp);
                tmp.assign(width, 0);
            }

            for (size_t i = 0; i != wires.size(); ++i) {

                size_t start_idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), wires[i]->get_startNode()));
                size_t end_idx = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), wires[i]->get_endNode()));

                tmp[start_idx] -= 1;
                tmp[end_idx] += 1;
                tmp[width - 1] = wires[i]->get_vol();

                matrix.push_back(tmp);
                tmp.assign(width, 0);
            }
        
            std::vector<double> solution = gauss_solve(matrix);

            for (size_t i = 0; i != nodes.size(); ++i){
                nodes[i]->set_potential(solution[i]);
            }

            for (size_t i = 0; i != wires.size(); ++i){
                wires[i]->set_cur(solution[i + nodes.size()]);
            }

            for (size_t i = 0; i != wires.size(); ++i) {
                wires[i]->vac(time, step);
            }

            for (size_t i = 0; i != resistors.size(); ++i) {
                resistors[i]->vac(time, step);
            }

            for (size_t i = 0; i != cur_sources.size(); ++i) {
                cur_sources[i]->vac(time, step);
            }

        matrix.clear();
        tmp.clear();
        }

        virtual ~Circuit() {
            delete_vector_elements<Node*>(nodes);
            delete_vector_elements<Bar*>(resistors);
            delete_vector_elements<Bar*>(cur_sources);
            delete_vector_elements<Bar*>(wires);
        }
};


struct Simulation
{
    private:
        double total_time, step;
        double time = 0;
        Circuit* circuit;
    public:
        Simulation(Circuit* circuit_, double total_time_, double step_) {
            total_time = total_time_;
            step = step_;
            circuit = circuit_;
        }

        void test_run(Bar* test_bar) {
            std::cout << "[";
            for (time = 0; time <= total_time; time += step) {
                circuit->solve(time, step);
                std::cout << test_bar->get_vol() << std::endl;
            }
            std::cout << "]";
            time = 0;
        }

        void run() {
            for (time = 0; time <= total_time; time += step) {
            circuit->solve(time, step);
            }
            time = 0;
        }

    virtual ~Simulation() {
        delete this->circuit;
    }
};

//void readDataAndCreateCircuit(const std::string& filename) {
//    std::ifstream inputFile(filename);
//
//    if (!inputFile.is_open()) {
//        std::cerr << "Unable to open file: " << filename << std::endl;
//        return;
//    }
//
//    std::vector<Node*> nodes;  // Вектор узлов
//
//    std::vector<CircuitElement> elements;
//    while (inputFile >> elementType) {
//        if (elementType == "0") {
//            break;
//        }
//    }
//}
//
//
//void addNodesToCircuit(Circuit& circuit, const std::vector<Node*>& nodes) {
//    for (const auto& node : nodes) {
//        circuit.add_node(node);
//    }
//}
//
//void addElementsToCircuit(Circuit& circuit, const std::vector<Node*>& nodes, const std::vector<CircuitElement>& elements) {
//    for (const auto& element : elements) {
//        if (element.type == "Resistor") {
//            Bar* resistor = new Resistor(element.value, nodes[element.node1], nodes[element.node2]);
//            circuit.add_bar(resistor);
//        }
//    }
//}


// FIXME Элементы цепи могут принимать разное значение аргументов, поэтому надо сделать так чтобы их можно было задавать разное количество.
void readData(std::ifstream& file, std::vector<std::string>& types, std::vector<double>& values, std::vector<std::string>& nodes1, std::vector<std::string>& nodes2) {
    std::string line;

    while (getline(file, line)) {
        std::istringstream ss(line);
        std::string type;
        double value;
        std::string node0;
        std::string node1;
        std::string node2;

        ss >> type;
        if (type == "0") {
            break;
        }
        ss.ignore();

        ss >> value >> node0 >> node1 >> node2;

        types.push_back(type);
        values.push_back(value);
        nodes1.push_back(node1);
        nodes2.push_back(node2);
    }
}

// Функция для создания узлов и элементов цепи на основе данных из файлов nodes1, nodes2, types и values
void createCircuitFromData(std::vector<std::string>& nodes1, std::vector<std::string>& nodes2,
    std::vector<std::string>& types, std::vector<double>& values, Circuit* circuit) {
    std::vector<Node*> nodes;
    std::map<std::string, Node*> node_map;

    for (size_t i = 0; i < nodes1.size(); ++i) {
        if (node_map.find(nodes1[i]) == node_map.end()) {
            Node* new_node = new Node();
            nodes.push_back(new_node);
            node_map[nodes1[i]] = new_node;
            circuit->add_node(new_node);
        }
        if (node_map.find(nodes2[i]) == node_map.end()) {
            Node* new_node = new Node();
            nodes.push_back(new_node);
            node_map[nodes2[i]] = new_node;
            circuit->add_node(new_node);
        }
    }

    for (size_t i = 0; i < types.size(); ++i) {
        if (types[i] == "Resistor") {
            Bar* new_resistor = new Resistor(values[i], node_map[nodes1[i]], node_map[nodes2[i]]);
            circuit->add_bar(new_resistor);
        }
        else if (types[i] == "Voltage_source") {
            Bar* new_voltage_source = new Voltage_source(values[i], node_map[nodes1[i]], node_map[nodes2[i]]);
            circuit->add_bar(new_voltage_source);
        }
        else if (types[i] == "Inductor") {
            Bar* new_inductor = new Inductor(values[i], 0, node_map[nodes1[i]], node_map[nodes2[i]]);
            circuit->add_bar(new_inductor);
        }
        else if (types[i] == "Current_source") {
            Bar* new_voltage_source = new Current_source(values[i], node_map[nodes1[i]], node_map[nodes2[i]]);
            circuit->add_bar(new_voltage_source);
        }
        else if (types[i] == "Capacitor") {
            Bar* new_inductor = new Capacitor(values[i], 0, node_map[nodes1[i]], node_map[nodes2[i]]);
            circuit->add_bar(new_inductor);
        }
        // Добавляем другие типы элементов цепи...
    }
}

int main() {
    
    std::ifstream file("input.txt");
    std::vector<std::string> types;
    std::vector<double> values;
    std::vector<std::string> nodes1;
    std::vector<std::string> nodes2;

    if (file.is_open()) {
        readData(file, types, values, nodes1, nodes2);
        file.close();

        //тест на вывод
        for (size_t i = 0; i < types.size(); ++i) {
            std::cout << "Type: " << types[i] << ", Value: " << values[i] << ", Node1: " << nodes1[i] << ", Node2: " << nodes2[i] << std::endl;
        }
    }

    std::cout << '\n';

    Circuit* c1 = new Circuit;

    createCircuitFromData(nodes1, nodes2, types, values, c1);

    c1->print();

    // a->print();
    // b->print();

    Simulation* sim1 = new Simulation(c1, 0.05, 0.0001);

    //c1->add_bar(v1);
    //c1->add_bar(l1);

    //c1->add_bar(r1);
    //c1->add_bar(r2);
    //c1->add_bar(r3);
    //c1->add_bar(r4);
    //c1->add_bar(r5);

    //c1->add_node(a);
    //c1->add_node(b);
    //c1->add_node(c);
    //c1->add_node(d);
    //c1->add_node(e);

    //sim1->test_run(l1);

    // c1->solve(0, 0);

    // std::cout << r5->get_cur() << std::endl;
    // std::cout << r2->get_cur() << std::endl;
    // std::cout << v1->get_cur() << std::endl;

    // d->print();

    sim1->run();

    delete sim1;

    return 0;
}
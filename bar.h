#ifndef BAR_H
#define BAR_H

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

#endif
#include <unordered_map>

// Stores data structure for different types of particle-based methods (atoms, coarse-grained molecules, etc)

struct Atom {
    int id, type;
    double x, y, z, ke, pe, sxx, syy, szz, sxy, sxz, syz;
    int bulk_coordination, coordination;
    std::string neighbor_list;
    double q, dqx, dqy, dqz;
    double Fxx, Fxy, Fxz, Fyx, Fyy, Fyz, Fzx, Fzy, Fzz;

    Atom() : id(0), type(1), x(0.0), y(0.0), z(0.0), 
             ke(0.0), pe(0.0), 
             sxx(0.0), syy(0.0), szz(0.0), sxy(0.0), sxz(0.0), syz(0.0),
             bulk_coordination(0), coordination(0),
             neighbor_list(""),
             q(0.0), dqx(0.0), dqy(0.0), dqz(0.0),
             Fxx(0.0), Fxy(0.0), Fxz(0.0), Fyx(0.0), Fyy(0.0), Fyz(0.0), Fzx(0.0), Fzy(0.0), Fzz(0.0) {}
};
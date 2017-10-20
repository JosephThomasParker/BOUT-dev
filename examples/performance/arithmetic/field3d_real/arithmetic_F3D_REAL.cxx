/*
 * Timing of arithmetic operations
 *
 */

#include <bout/physicsmodel.hxx>

#include <bout/expr.hxx>

#include <chrono>

typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
typedef std::chrono::duration<double> Duration;
using namespace std::chrono;

class Arithmetic : public PhysicsModel {
protected: 
  int init(bool restarting) {
    
    BoutReal a = 1.0;
    BoutReal b = 2.0;
    Field3D c = 3.0;

    Field3D result1, result2, result3, result4, result5;

    // Using Field methods (classic operator overloading)
    
    result1 = a + b * c;

    SteadyClock start1 = steady_clock::now();
    result1 = a + b * c;
    Duration elapsed1 = steady_clock::now() - start1;
    
    // Using C loops
    result2.allocate();
    BoutReal *rd = &result2(0,0,0);
    BoutReal *cd = &c(0,0,0);
    SteadyClock start2 = steady_clock::now();
    for(int i=0, iend=(mesh->LocalNx*mesh->LocalNy*mesh->LocalNz)-1; i != iend; i++) {
      *rd = a + b*(*cd);
      rd++;
      cd++;
    }
    Duration elapsed2 = steady_clock::now() - start2;
    
    // Template expressions
    SteadyClock start3 = steady_clock::now();
    result3 = eval3D(add(a, mul(b,c)));
    Duration elapsed3 = steady_clock::now() - start3;
    
    // Range iterator
    result4.allocate();
    SteadyClock start4 = steady_clock::now();
    for(auto i : result4)
      result4[i] = a + b * c[i];
    Duration elapsed4 = steady_clock::now() - start4;

    // SingleDataIterator
    result5.allocate();
    SteadyClock start5 = steady_clock::now();
    for(SingleDataIterator i = result5.sdi_region(RGN_ALL); !i.done(); ++i)
      result5(i) = a + b * c(i);
    Duration elapsed5 = steady_clock::now() - start5;
    
    output << "TIMING\n======\n";
    output << "Fields: " << elapsed1.count() << endl;
    output << "C loop: " << elapsed2.count() << endl;
    output << "Templates: " << elapsed3.count() << endl;
    output << "Range For: " << elapsed4.count() << endl;
    output << "SingleDataIterator For: " << elapsed5.count() << endl;
    
    return 1;
  }
};

BOUTMAIN(Arithmetic);

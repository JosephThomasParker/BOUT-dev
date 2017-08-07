/*
 * Testing performance of an iterator over the mesh
 *
 */

#include <bout.hxx>

#include <chrono>
#include <iostream>
#include <iterator>
#include <time.h>
#include <omp.h>

// A simple iterator over a 3D set of indices
class MeshIterator
  : public std::iterator< std::forward_iterator_tag, Indices > {
public:
  /// Constructor. This would set ranges. Could depend on thread number
  MeshIterator() : x(0), y(0), z(0), xstart(0), ystart(0), zstart(0) {
    xend = mesh->LocalNx-1;
    yend = mesh->LocalNy-1;
    zend = mesh->LocalNz;
  }

  MeshIterator(int x, int y, int z) : x(x), y(y), z(z), xstart(0), ystart(0), zstart(0) {
    xend = mesh->LocalNx-1;
    yend = mesh->LocalNy-1;
    zend = mesh->LocalNz;
  }

  /// The index variables, updated during loop
  int x, y, z;

  /// Increment operators
  MeshIterator& operator++() { next(); return *this; }
  MeshIterator& operator++(int) { next(); return *this; }

  // Comparison operator
  bool operator!=(const MeshIterator& rhs) const {
    return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
  }

  // Dereference operator
  Indices operator*() {
    return {x, y, z};
  }

  /// Checks if finished looping. Is this more efficient than
  /// using the more idiomatic it != MeshIterator::end() ?
  bool isDone() const {
    return x > xend;
  }

private:
  int xstart, xend;
  int ystart, yend;
  int zstart, zend;

  /// Advance to the next index
  void next() {
    z++;
    if(z > zend) {
      z = zstart;
      y++;
      if(y > yend) {
        y = ystart;
        x++;
      }
    }
  }
};

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  Field3D a = 1.0;
  Field3D b = 2.0;

  Field3D result;
  result.allocate();

  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;
  
  // A single loop over block data
  
  BoutReal *ad = &a(0,0,0);
  BoutReal *bd = &b(0,0,0);
  BoutReal *rd = &result(0,0,0);
  
///aaa///\  // Loop over data so first test doesn't have a disadvantage from caching
///aaa///\  for(int i=0;i<10;++i) {
///aaa///\#pragma omp parallel for
///aaa///\    for(int j=0;j<mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;++j) {
///aaa///\      rd[j] = ad[j] + bd[j];
///aaa///\    }
///aaa///\  }
///aaa///\  
///aaa///\  SteadyClock start1 = steady_clock::now();
///aaa///\  int len = mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;
///aaa///\  for(int i=0;i<10;++i) {
///aaa///\#pragma omp parallel for
///aaa///\    for(int j=0;j<len;++j) {
///aaa///\      rd[j] = ad[j] + bd[j];
///aaa///\    }
///aaa///\  }
///aaa///\  Duration elapsed1 = steady_clock::now() - start1;
///aaa///\
///aaa///\  // Nested loops over block data
///aaa///\  SteadyClock start2 = steady_clock::now();
///aaa///\///  for(int x=0;x<10;x++) {
///aaa///\///    for(int i=0;i<mesh->LocalNx;++i) {
///aaa///\///      for(int j=0;j<mesh->LocalNy;++j) {
///aaa///\///#pragma ivdep
///aaa///\///#pragma omp parallel for
///aaa///\///        for(int k=0;k<mesh->LocalNz;++k) {
///aaa///\///          result(i,j,k) = a(i,j,k) + b(i,j,k);
///aaa///\///        }
///aaa///\///      }
///aaa///\///    }
///aaa///\///  }
///aaa///\  Duration elapsed2 = steady_clock::now() - start2;
///aaa///\
///aaa///\  // MeshIterator over block data
///aaa///\  SteadyClock start3 = steady_clock::now();
///aaa///\  for(int x=0;x<10;x++) {
///aaa///\    for(MeshIterator i; !i.isDone(); ++i){
///aaa///\      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
///aaa///\    }
///aaa///\  }
///aaa///\  Duration elapsed3 = steady_clock::now() - start3;
///aaa///\
///aaa///\  // DataIterator using begin(), end()
///aaa///\  SteadyClock start4 = steady_clock::now();
///aaa///\  for(int x=0;x<10;x++) {
///aaa///\#pragma omp parallel
///aaa///\    {
///aaa///\    for(DataIterator i = result.beginDI(), rend=result.endDI(); i != rend; ++i){
///aaa///\    //for(DataIterator i = result.iterator(); !i.done(); ++i){
///aaa///\      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
///aaa///\    }
///aaa///\    }
///aaa///\}
///aaa///\  Duration elapsed4 = steady_clock::now() - start4;
///aaa///\
///aaa///\  // DataIterator with done()
///aaa///\  SteadyClock start5 = steady_clock::now();
///aaa///\  for(int x=0;x<10;x++) {
///aaa///\    //for(DataIterator i = begin(result); !i.done() ; ++i){
///aaa///\#pragma omp parallel
///aaa///\    {
///aaa///\    for(DataIterator i = result.iterator(); !i.done(); ++i){
///aaa///\      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
///aaa///\    }
///aaa///\    }
///aaa///\  }
///aaa///\  Duration elapsed5 = steady_clock::now() - start5;
///aaa///\
///aaa///\  // Range based for DataIterator with indices
///aaa///\///  SteadyClock start6 = steady_clock::now();
///aaa///\///  for(int x=0;x<10;x++) {
///aaa///\///    for(auto i : result){
///aaa///\///      result(i.x,i.y,i.z) = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
///aaa///\///    }
///aaa///\///  }
///aaa///\///  Duration elapsed6 = steady_clock::now() - start6;
///aaa///\
///aaa///\  // Range based for with single index
///aaa///\  SteadyClock start7 = steady_clock::now();
///aaa///\///  for(int x=0;x<10;x++) {
///aaa///\///#pragma ivdep
///aaa///\///#pragma omp parallel
///aaa///\///{
///aaa///\///#pragma omp single
///aaa///\///  {
///aaa///\///    for(auto i : result){
///aaa///\///    //for(auto &i : result){
///aaa///\///#pragma omp task
///aaa///\///    {
///aaa///\///      result[i] = a[i] + b[i];
///aaa///\///    }
///aaa///\///    }
///aaa///\///  }
///aaa///\///}
///aaa///\///}
///aaa///\  Duration elapsed7 = steady_clock::now() - start7;
///aaa///\
///aaa///\  // Range based for with single index
///aaa///\  SteadyClock start8 = steady_clock::now();
///aaa///\  for(int x=0;x<10;x++) {
///aaa///\#pragma ivdep
///aaa///\#pragma omp parallel
///aaa///\    {
///aaa///\    for(auto &i : result.region(RGN_ALL)){
///aaa///\      //output << i.x << " " << i.y << " " << i.z <<  "\n";
///aaa///\      result[i] = a[i] + b[i];
///aaa///\    }
///aaa///\    }
///aaa///\  }
///aaa///\  Duration elapsed8 = steady_clock::now() - start8;
  
  // Range based DataIterator 
  SteadyClock start9 = steady_clock::now();
  for (int x=0;x<10;++x) {
#pragma ivdep
#pragma omp parallel
    {
    //for (const auto &i : result) {
    //for(SingleDataIterator i = result.Siterator(); !i.done(); ++i){
    //for(SingleDataIterator i = result.sdi_region(RGN_ALL); !i.done(); ++i){
    for(SingleDataIterator i = result.sdi_region_all(); !i.done(); ++i){
    //for(SingleDataIterator i = result.sdi_region(RGN_NOY); !i.done(); ++i){
    //for(SingleDataIterator i = result.sdi_region(RGN_NOX); !i.done(); ++i){
      //output << "inside loop " << i.x ;
///      if( omp_get_thread_num() == 1 ){
///	//output << i.x << " " << i.y << " " << i.z << " " << omp_get_thread_num() <<  "\n";
///        output << i.x << " " << omp_get_thread_num() <<  "\n";
///      }
      //output << i.nx << " " << i.ny << " " << i.nz <<  "\n";
      //if( omp_get_thread_num() == 1 ){
      //	output << i.i << " " << i.x << " " << i.y << " " << i.z <<  "\n";
      //}
      output << "Performing iteration: " << i.i << ", with count: " << i.icount << "\n";
      result(i) = a(i) + b(i);
      //output << i.icount << "\n"; 
      //result[i.x] = a[i.x] + b[i.x];
    }
    }
  }
  Duration elapsed9 = steady_clock::now() - start9;
///aaa///\
///aaa///\  // DataIterator over fields
///aaa///\///  SteadyClock start10 = steady_clock::now();
///aaa///\///  for(int x=0;x<10;x++)
///aaa///\///    //for(DataIterator d = DataIterator(0,mesh->LocalNx,0,mesh->LocalNy,0,mesh->LocalNz); !d.done(); ++d)
///aaa///\///    for(DataIterator d = std::begin(result) ; !d.done(); ++d)
///aaa///\///    //for(DataIterator d = DataIterator(0,mesh->LocalNx,0,mesh->LocalNy,0,mesh->LocalNz); d <= d.end(); ++d)
///aaa///\///      result[d] = a[d] + b[d];
///aaa///\///  Duration elapsed10 = steady_clock::now() - start10;
///aaa///\  // DataIterator over fields
///aaa///\  SteadyClock start10 = steady_clock::now();
///aaa///\  for(int x=0;x<10;x++)
///aaa///\#pragma omp parallel
///aaa///\{
///aaa///\    for(DataIterator d = result.iterator(); !d.done(); d++)
///aaa///\      result[d] = a[d] + b[d];
///aaa///\}
///aaa///\  Duration elapsed10 = steady_clock::now() - start10;
///aaa///\
///aaa///\  
///aaa///\  output << "TIMING\n======\n";
///aaa///\  output << "C loop                     : " << elapsed1.count() << std::endl;
///aaa///\  output << "----- (x,y,z) indexing ----" << std::endl;
///aaa///\  output << "Nested loops               : " << elapsed2.count() << std::endl;
///aaa///\  output << "MeshIterator               : " << elapsed3.count() << std::endl;
///aaa///\  output << "DataIterator (begin/end)   : " << elapsed4.count() << std::endl;
///aaa///\  output << "DataIterator (begin/done)  : " << elapsed5.count() << std::endl;
///aaa///\///  output << "C++11 range-based for      : " << elapsed6.count() << std::endl;
///aaa///\  output << "------ [i] indexing -------" << std::endl;
///aaa///\  output << "Single index               : " << elapsed7.count() << std::endl;
///aaa///\  output << "Three indices              : " << elapsed8.count() << std::endl;
///aaa///\  output << "C++11 Range-based for      : " << elapsed9.count() << std::endl;
///aaa///\  output << "DataIterator (begin/done)  : " << elapsed10.count() << std::endl;
  BoutFinalise();
  return 0;
}

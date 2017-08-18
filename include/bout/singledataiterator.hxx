/*


 */

#ifndef __SINGLEDATAITERATOR_H__
#define __SINGLEDATAITERATOR_H__

#include <iterator>
#include <iostream>
#include "unused.hxx"
#include <output.hxx>

#ifdef _OPENMP
#include <omp.h>
int SDI_spread_work(int num_work, int thread, int max_thread);
#endif

/*!
 * Set of indices - DataIterator is dereferenced into these
 */
struct SIndices {
  int i; // index of array
  int j; // count index
  int nx;
  int ny;
  int nz;
};

#define DI_GET_END ((void *) NULL)

/*!
 * Provides range-based iteration over indices. 
 * If OpenMP is enabled, then this divides work between threads.
 * 
 * This is used mainly to loop over the indices of fields,
 * and provides convenient ways to index 
 * 
 * Example
 * -------
 * 
 * Start,end values for (x,y,z) can be specified directly:
 * 
 *     for(d = DataIterator(xs, xe, ys, ye, zs, ze); !d.done(); ++d) {
 *       // print index
 *       output.write("%d,%d,%d\n", d.x, d.y, d.z);
 *       // Index into a Field3D variable 'f'
 *       output.write("Value = %e\n", f[d]);
 *     }
 * 
 * Usually DataIterator is used to loop over fields. Field3D::begin()
 * and Field3D::end() return DataIterator objects:
 * 
 *     Field3D f(0.0); // Initialise field
 *     for(auto i : f) { // Loop over all indices, including guard cells
 *       f[i] = i.x; // Indexing using DataIterator
 *     }
 * 
 */
class SingleDataIterator
  : public std::iterator<std::bidirectional_iterator_tag, SIndices> {
private:
  /*!
   * This initialises OpenMP threads if enabled, and
   * divides iteration index ranges between threads
   */
#ifdef _OPENMP
  void omp_init(bool end);
#endif
  void make_region(int xs,int xe,
                         int ys,int ye,
			 int zs,int ze,
			 int nx, int ny, int nz);
  void idx_to_xyz(int i);

public:
  std::vector<int> make_region_all(int xs,int xe,
                         int ys,int ye,
			 int zs,int ze,
			 int nx, int ny, int nz);
  /*!
   * Constructor. This sets index ranges.
   * If OpenMP is enabled, the index range is divided
   * between threads using the omp_init method.
   */ 
  SingleDataIterator(int xs, int xe,
	       int ys, int ye,
	       int zs, int ze,
	       int nx, int ny, int nz, std::vector<int>& rgn) : 
#ifndef _OPENMP
    x(xs), y(ys), z(zs),
    xstart(xs),   ystart(ys),   zstart(zs),
    xmin(xstart), ymin(ystart), zmin(zstart),
    xend(xe),     yend(ye),     zend(ze),
    xmax(xend),   ymax(yend),   zmax(zend),
    i((x*ny+y)*nz+z), istart((xs*ny+ys)*nz+zs), imin(istart), iend((xe*ny+ye)*nz+ze), imax(iend),
    icount(0),
    icountstart(0),
    icountend((xe-xs+1)*(ye-ys+1)*(ze-zs+1)),
    nx(nx), ny(ny), nz(nz),
#else
    xmin(xs),     ymin(ys),     zmin(zs),
    xmax(xe),     ymax(ye),     zmax(ze),
    imin((xs*ny+ys)*nz+zs), imax((xe*ny+ye)*nz+ze), 
    icount(0),
    icountstart(0),
    icountend((xe-xs+1)*(ye-ys+1)*(ze-zs+1)),
    nx(nx), ny(ny), nz(nz),
#endif
    rgn(rgn),
    isEnd(false)
  {
#ifdef _OPENMP
    omp_init(false);
#endif
  }

  /*!
   * set end();
   * use as DataIterator(int,int,int,int,int,int,DI_GET_END);
   */
  SingleDataIterator(int xs, int xe,
	       int ys, int ye,
	       int zs, int ze,
	       int nx, int ny, int nz, std::vector<int>& rgn, void* UNUSED(dummy)) : 
#ifndef _OPENMP
    x(xs), y(ys), z(zs),
    xstart(xs),   ystart(ys),   zstart(zs),
    xmin(xstart), ymin(ystart), zmin(zstart),
    xend(xe),     yend(ye),     zend(ze),
    xmax(xend),   ymax(yend),   zmax(zend),
    i((x*ny+y)*nz+z), istart((xs*ny+ys)*nz+zs), imin(istart), iend((xe*ny+ye)*nz+ze), imax(iend),
    icount(0),
    icountstart(0),
    icountend((xe-xs+1)*(ye-ys+1)*(ze-zs+1)),
#else
    xmin(xs),     ymin(ys),     zmin(zs),
    xmax(xe),     ymax(ye),     zmax(ze),
    imin((xs*ny+ys)*nz+zs), imax((xe*ny+ye)*nz+ze), 
    icount(0),
    icountstart(0),
    icountend((xe-xs+1)*(ye-ys+1)*(ze-zs+1)),
#endif
    rgn(rgn),
    isEnd(true)
  {
#ifdef _OPENMP
    omp_init(true);
#endif
    //rgn = new int [icountend];
    //make_region(xstart,xend,ystart,yend,zstart,zend,nx,ny,nz);
    
    next();
  }
  
  /*!
   * The index variables, updated during loop
   * Should make these private and provide getters?
   */
  int i, icount;
  int icountstart, icountend;
  int x, y, z;
  int nx, ny, nz;
  const std::vector<int>& rgn;

  /// Pre-increment operator. Use this rather than post-increment when possible
  SingleDataIterator& operator++() { next(); return *this; }
  
  /// Post-increment operator
  SingleDataIterator operator++(int) { SingleDataIterator tmp(*this); next(); return tmp; }
  
  /// Pre-decrement operator
  SingleDataIterator& operator--() { prev(); return *this; }
  
  /// Post-decrement operator
  SingleDataIterator operator--(int) { SingleDataIterator tmp(*this); prev(); return tmp; }

  /// Comparison operator. Most common use is in for loops
  inline bool operator!=(const SingleDataIterator& rhs) const {
    //return  !(x == rhs.x && y == rhs.y && z == rhs.z);
    if (rhs.isEnd){
      return !this->done();
    } else {
      return  !(i == rhs.i);
    }
  }
  
  /*!
   * Dereference operators
   * These are needed because the C++11 for loop
   * dereferences the iterator
   */
  SingleDataIterator& operator*() {
    return *this;
  }
  
  /*!
   * Const dereference operator. 
   * Needed because C++11 for loop dereferences the iterator
   */
  const SingleDataIterator& operator*() const {
    return *this;
  }

  /*!
   * Add an offset to the index for general stencils
   */
///  const SIndices offset(int dx, int dy, int dz) const {
///    if (dz>0){
///      int zp=rgn[icount]%nz;
///      for (int j=0;j<dz;++j)
///        zp=(zp == nz-1 ? 0 : zp+1);
///      return { rgn[icount] + ny*nz*dx + nz*dy + zp , nx, ny, nz };
///    } else {
///      int zm=rgn[icount]%nz;
///      for (;dz!= 0;++dz)
///        zm = (zm == 0 ? nz-1 : zm-1);
///      return { rgn[icount] + ny*nz*dx + nz*dy + zm , nx, ny, nz};
///    }
///  }
  
  /*
   * Shortcuts for common offsets, one cell
   * in each direction.
   */
  
///  /// The index one point +1 in x
///  //const SIndices xp() const { return { i + ny*nz , nx, ny, nz}; }
///  const SIndices xp() const { return { rgn[icount] + ny*nz , nx, ny, nz}; }
///  /// The index one point -1 in x
///  //const SIndices xm() const { return { i - ny*nz , nx, ny, nz}; }
///  const SIndices xm() const { return { rgn[icount] - ny*nz , nx, ny, nz}; }
///  /// The index one point +1 in y
///  //const SIndices yp() const { return { i + nz , nx, ny, nz}; }
///  const SIndices yp() const { return { rgn[icount] + nz , nx, ny, nz}; }
///  /// The index one point -1 in y
///  //const SIndices ym() const { return { i - nz , nx, ny, nz}; }
///  const SIndices ym() const { return { rgn[icount] - nz , nx, ny, nz}; }
///  /// The index one point +1 in z. Wraps around zend to zstart
///  //const SIndices zp() const { return { (i+1)%nz == 0 ? i-nz+1 : i+1 , nx, ny, nz}; }
///  const SIndices zp() const { return { (rgn[icount]+1)%nz == 0 ? rgn[icount]-nz+1 : rgn[icount]+1 , nx, ny, nz}; }
///  /// The index one point -1 in z. Wraps around zstart to zend
///  //const SIndices zm() const { return { i%nz == 0 ? i+nz-1 : i-1 , nx, ny, nz }; }
///  const SIndices zm() const { return { rgn[icount]%nz == 0 ? rgn[icount]+nz-1 : rgn[icount]-1 , nx, ny, nz}; }

  /*!
   * Resets DataIterator to the start of the range
   */
  void start() {
    i = istart;
  }

  /*!
   * Resets DataIterator to the start of the range
   */
  void begin() {
    i = istart;
  }

  /*!
   * Sets DataIterator to one index past the end of the range
   */ 
  void end() {
    i = iend;
    next();
  }

  /*!
   * Checks if finished looping. Is this more efficient than
   * using the more idiomatic it != DataIterator::end() ?
   */
  bool done() const {
    return icount >= icountend ;
///#ifndef _OPENMP
///    //return (i > iend) || (i < istart);
///#else //_OPENMP
///    return (i > iend) || (i < istart);
///    //return (x == xend) || x > xend || (x <= xstart)  ;
///#endif //_OPENMP
  }
  
private:
  SingleDataIterator(); // Disable null constructor

  //const int nx, ny, nz;
#ifndef _OPENMP
  const int istart;
  const int xstart, ystart, zstart;
#else
  int istart;
  int xstart, ystart, zstart;
#endif

  int imin;
  int xmin, ymin, zmin;

#ifndef _OPENMP
  const int iend;
  const int xend, yend, zend;
#else
  int iend;
  int xend, yend, zend;
#endif

  int imax;
  int xmax, ymax, zmax;

  const bool isEnd;
  /// Advance to the next index
  void next() {
    icount++;
    //i = rgn[icount];
    //idx_to_xyz(i);
  }

  /// Rewind to the previous index
  void prev() {
    icount--;
    //i = rgn[icount];
  }
};

/*!
 * Specifies a range of indices which can be iterated over
 * and begin() and end() methods for range-based for loops
 * 
 * Example
 * -------
 *
 * Index ranges can be defined manually:
 *
 *     IndexRange r(0, 10, 0, 20, 0, 30);
 *     
 * then iterated over using begin() and end()
 *
 *     for( DataIterator i = r.begin(); i != r.end(); i++ ) {
 *       output.write("%d,%d,%d\n", i.x, i.y, i.z);
 *     }
 *
 * or the more convenient range for loop:
 *
 *     for( auto i : r ) {
 *       output.write("%d,%d,%d\n", i.x, i.y, i.z);
 *     }
 *
 * A common use for this class is to loop over
 * regions of a field:
 *
 *     Field3D f(0.0);
 *     for( auto i : f.region(RGN_NOBNDRY) ) {
 *       f[i] = 1.0;
 *     }
 * 
 * where RGN_NOBNDRY specifies a region not including
 * boundary/guard cells. 
 */
struct SIndexRange {
  int istart, nx, ny, nz;
  int iend = nx*ny*nz-1;
  int xstart, xend;
  int ystart, yend;
  int zstart, zend;
  std::vector<int>& rgn;
  
  const SingleDataIterator begin() const {
    return SingleDataIterator(xstart, xend, 
                              ystart, yend,
                              zstart, zend,
                              nx, ny, nz, rgn);
  }
  const SingleDataIterator end() const {
    return SingleDataIterator(xstart, xend, 
			      ystart, yend,
			      zstart, zend,
			      nx, ny, nz, rgn, DI_GET_END);
  }
};

inline void SingleDataIterator::idx_to_xyz(int i){
  // function for debugging
  // print x,y,z for a given icount

  // i = (x*ny+y)*nz+z
  output << "i = " << i << ", x = " << ((i/nz)/ny) << ", y = " << (i/nz)%ny << ", z = " << (i%nz) << "\n";
};

inline std::vector<int> SingleDataIterator::make_region_all(int xstart,int xend,
                                            int ystart,int yend,
		                            int zstart,int zend,
		                            int nx,int ny,int nz){
  // Make an array of indices corresponding to a region.

  int j=0;
  int x = xstart;
  int y = ystart;
  int z = zstart;
  std::vector<int> rgn;
  
///  if( omp_get_thread_num() == 1 ){
///    output<<"starting xyz:"<<x<<" "<<y<<" "<<z<<" "<<xend<<" "<<yend<<" "<<zend<<"\n";
///  }

  bool done = false;
  j=-1;
  while( !done ){
      j++;
      rgn[j] = (x*ny+y)*nz+z;
///      if( omp_get_thread_num() == 1 ){
///	output << rgn[j] << " " << j << ", xy index: " << rgn[j]/nz << ", x index: "<< ((rgn[j]/nz)/ny) << ", y index = " << (rgn[j]/nz)%ny <<  ", z index: " << rgn[j]%nz  << ", routine's x,y,z " << x << " " << y << " " << z << ", xend, yend, zend:"<<xend<<" "<<yend<<" "<<zend<< "\n" <<std::flush;
///      }
      if(x == xend && y == yend && z == zend){
	icountend = j+1;
	done = true;
      }
      ++z;
      if(z > zmax) {
	z = zmin;
	++y;
	if(y > ymax) {
	  y = ymin;
	  ++x;
	}
      }
    }
  return rgn;
};

///inline void SingleDataIterator::make_region(int xstart,int xend,
///                                            int ystart,int yend,
///		                            int zstart,int zend,
///		                            int nx,int ny,int nz){
///  // Make an array of indices corresponding to a region.
///
///  int j=0;
///  int x = xstart;
///  int y = ystart;
///  int z = zstart;
///  
//////  if( omp_get_thread_num() == 1 ){
//////    output<<"starting xyz:"<<x<<" "<<y<<" "<<z<<" "<<xend<<" "<<yend<<" "<<zend<<"\n";
//////  }
///
///  bool done = false;
///  j=-1;
///  while( !done ){
///      j++;
///      rgn[j] = (x*ny+y)*nz+z;
//////      if( omp_get_thread_num() == 1 ){
//////	output << rgn[j] << " " << j << ", xy index: " << rgn[j]/nz << ", x index: "<< ((rgn[j]/nz)/ny) << ", y index = " << (rgn[j]/nz)%ny <<  ", z index: " << rgn[j]%nz  << ", routine's x,y,z " << x << " " << y << " " << z << ", xend, yend, zend:"<<xend<<" "<<yend<<" "<<zend<< "\n" <<std::flush;
//////      }
///      if(x == xend && y == yend && z == zend){
///	icountend = j+1;
///	done = true;
///      }
///      ++z;
///      if(z > zmax) {
///	z = zmin;
///	++y;
///	if(y > ymax) {
///	  y = ymin;
///	  ++x;
///	}
///      }
///    }
///};

#ifdef _OPENMP
inline int SDI_spread_work(int work,int cp,int np){
  // Spread work between threads. If number of points do not
  // spread evenly between threads, put the remaining "rest"
  // points on the threads 0 ... rest-1.
  int pp=work/np;
  int rest=work%np;
  int result=pp*cp;
  if (rest > cp){
    result +=cp;
  } else {
    result +=rest;
  }
  return result;
};

inline void SingleDataIterator::omp_init(bool end){
  // In the case of OPENMP we need to calculate the range
  int threads=omp_get_num_threads();
  if (threads > 1){
    int ny=ymax-ymin+1;
    int nz=zmax-zmin+1;
    int work  = (xmax-xmin+1)*ny*nz;
    int current_thread = omp_get_thread_num();
    int begin = SDI_spread_work(work,current_thread,threads);
    int end   = SDI_spread_work(work,current_thread+1,threads);
    //icountend   = end - begin;
    icountend   = end;
    icount      = begin;
///#pragma omp critical
///    {
///    output << "thread: "<<omp_get_thread_num()<<", begin: "<<begin<<", end: "<<end<<", icountend: "<<icountend<<"\n";
///    }
    --end;
    zend   = (end   % nz) + zmin;
    zstart = (begin % nz) + zmin;
    end   /= nz;
    begin /= nz;
    yend   = (end   % ny) + ymin;
    ystart = (begin % ny) + ymin;
    end   /= ny;
    begin /= ny;
    xend   = end;
    xstart = begin;
  } else {
    zstart = zmin;
    zend   = zmax;
    ystart = ymin;
    yend   = ymax;
    xstart = xmin;
    xend   = xmax;
  }
  if (!end){
    x=xstart;
    y=ystart;
    z=zstart;
  } else {
    x=xend;
    y=yend;
    z=zend;
  }
};
#endif

#endif // __SINGLEDATAITERATOR_H__

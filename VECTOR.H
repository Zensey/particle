//
// vector.h
// 3D vector & plane class
//

#ifndef __vector_h__
#define __vector_h__

// Includes

#include <math.h>

// Utility




typedef long double REAL;

static REAL maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


inline REAL RandFlt(REAL min = 0.0f, REAL max = 1.0f)
{
  return min + (max - min) * rand() / (REAL)RAND_MAX;
}

// Vector

struct VECTOR {

  REAL x, y, z;

  VECTOR() {}

  VECTOR(REAL _x, REAL _y, REAL _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }

  
  friend inline VECTOR operator +(const VECTOR& a, const VECTOR & b) 
  {
    return VECTOR( a.x + b.x, a.y + b.y, a.z + b.z );
  }
  
  friend inline VECTOR & operator +=(VECTOR & a, const VECTOR & b)
  {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
  }
  
  friend inline VECTOR operator -(const VECTOR& a, const VECTOR & b) 
  {
    return VECTOR( a.x - b.x, a.y - b.y, a.z - b.z );
  }
  
  friend inline VECTOR & operator -=(VECTOR & a, const VECTOR & b)
  {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
  }
  
  friend inline VECTOR operator *(const VECTOR& a, REAL s)
  {
    return VECTOR( a.x * s, a.y * s, a.z * s );
  }
  
  friend inline REAL DotProd(const VECTOR & a, const VECTOR & b)
  {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  
  friend inline REAL Length2(const VECTOR & a)
  {
    return DotProd( a, a );
  }
  
  friend inline REAL Length(const VECTOR & a)
  {
    return (REAL)sqrtl( Length2( a ) );
  }
  
  friend inline REAL Distance(const VECTOR & a, const VECTOR & b)
  {
    return Length( a - b );
  }
  
  friend inline VECTOR Normalize(const VECTOR & a )
  {
    return a * (1.0f / Length( a ));
  }

  friend void Null(VECTOR & a )
  {
    a.x = 0;
    a.y = 0;
    a.z = 0;
  }
  
};

// Plane

struct PLANE {
  
  VECTOR n;
  REAL k;
  
  PLANE() {}
  
  PLANE(const VECTOR & _n, REAL _k)
  {
    n = _n;
    k = _k;
  }
  
  friend inline REAL DotProd(const VECTOR & v, const PLANE & p)
  {
    return DotProd( v, p.n ) - p.k;
  }
  
};

#endif



struct ATOMS
{

  float charge;
  VECTOR pos[175+175];
};

struct PARTICLES
{
  REAL mass, charge;
  VECTOR *vel;
  VECTOR *pos;
  bool *finished;

  PARTICLES (int n)
  {
    vel = new VECTOR[n];
    pos = new VECTOR[n];
    finished = new bool[n];
	for (int i = 0; i < n; i++) finished[i] = FALSE;
  }
  int DESTROY()
  {
    delete vel;
    delete pos;

	return 0;

  }
};

class CPARTICLES
{
 public:
  int N;
  float mass, charge;
  VECTOR *vel;
  VECTOR *pos;
  bool *finished;

  CPARTICLES (int n)
  {
    N = n;
    vel = new VECTOR[n];
    pos = new VECTOR[n];
    finished = new bool[n];
    int i;
    for (i = 0; i < n; i++) finished[i] = false;

  }
  int DESTROY()
  {
    delete vel;
    delete pos;

	return 0;
  }
};


















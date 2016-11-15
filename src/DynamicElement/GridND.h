
#ifndef GRIDND_H
#define GRIDND_H

#include "Grid.h"

#ifndef TOROIDAL_BOUNDARY
#define TOROIDAL_BOUNDARY
#endif

#define CLAMP(A, AMIN, AMAX) A = max((AMIN), min((AMAX), (A)))

template<uint T_dim, class T_type>
class CGridND
{
public:
	CGridND()
	{
		m_size = 0;
	}

	CGridND(Vec<T_dim, uint> size)
	{
		m_size = size;
		m_data.Allocate(productOfAll(m_size));
	}

	~CGridND() {}

	void Allocate(Vec<T_dim, uint> size)
	{
		m_size = size;
		m_data.Allocate(productOfAll(m_size));
	}

	void Allocate(uint i, uint j)
	{
		assert(T_dim == 2);
		m_size[0] = i;
		m_size[1] = j;
		m_data.Allocate(productOfAll(m_size));
	}

	void Allocate(uint i, uint j, uint k)
	{
		assert(T_dim == 3);
		m_size[0] = i;
		m_size[1] = j;
		m_size[2] = k;
		m_data.Allocate(productOfAll(m_size));
	}

	/// Copy (constructor)
	CGridND(const CGridND& a)
	{
		m_size = a.m_size;
		m_data = a.m_data;
	}

	/// Copy (affectation)
	const CGridND& operator = (const CGridND& a)
	{
		m_size = a.m_size;
		m_data = a.m_data;
		return (*this);
	}

	/// Fill array with a given value
	void Fill(T_type value_to_fill_with)
	{
		m_data.Fill(value_to_fill_with);
	}

	// Fill array with a given array
	void Fill(const T_type* data_to_fill_with)
	{
		m_data.Fill(data_to_fill_with);
	}

	/// Read only access
	const T_type& Get(const Vec<T_dim,uint>& idx) const
	{
		return (At(idx));
	}

	/// Read only access
	const T_type& At(const Vec<T_dim,uint>& idx) const
	{
		return (m_data[addressOf(idx,m_size)]);
	}

	/// Read/write access
	T_type& Set(const Vec<T_dim,uint>& idx)
	{
		return (At(idx));
	}

	/// Read/write access
	T_type& At(const Vec<T_dim,uint>& idx)
	{
		return (m_data[addressOf(idx,m_size)]);
	}

	/// [] operator with Tuple uint
	const T_type& operator[](const Vec<T_dim,uint>& idx) const
	{
		return At(idx);
	}

	/// [] operator with Tuple int
	const T_type& operator[](const Vec<T_dim,int>& idx) const
	{
		return At(Vec<T_dim,uint>(idx));
	}

	T_type& operator[](const Vec<T_dim,int>& idx)
	{
		return At(Vec<T_dim,uint>(idx));
	}

	T_type& operator[](const Vec<T_dim,uint>& idx)
	{
		return At(idx);
	}

	const T_type& operator()(uint i, uint j) const
	{
		assert(T_dim == 2);
		Vec<T_dim,uint> idx;
		idx[0] = i;
		idx[1] = j;
		return At(idx);
	}

	T_type& operator()(uint i, uint j)
	{
		assert(T_dim == 2);
		Vec<T_dim,uint> idx;
		idx[0] = i;
		idx[1] = j;
		return At(idx);
	}

	// Linear interpolation, range of x and y: [0, 1)
	T_type Sample(float x, float y) const
	{
		assert(T_dim == 2);
		Vec<T_dim,uint> idx;
		float nx = x * float(m_size[0]);
		float ny = y * float(m_size[1]);
		int ni1 = floor(nx);
		int nj1 = floor(ny);
		float wx2 = nx - ni1;
		float wy2 = ny - nj1;
		float wx1 = 1.0 - wx2;
		float wy1 = 1.0 - wy2;
		ni1 = (ni1 + m_size[0]) % m_size[0];
		nj1 = (nj1 + m_size[1]) % m_size[1];
		int ni2 = (ni1 + 1 + m_size[0]) % m_size[0];
		int nj2 = (nj1 + 1 + m_size[1]) % m_size[1];
		idx[0] = ni1;
		idx[1] = nj1;
		T_type v11 = At(idx) * wx1 * wy1;
		idx[0] = ni1;
		idx[1] = nj2;
		T_type v12 = At(idx) * wx1 * wy2;
		idx[0] = ni2;
		idx[1] = nj1;
		T_type v21 = At(idx) * wx2 * wy1;
		idx[0] = ni2;
		idx[1] = nj2;
		T_type v22 = At(idx) * wx2 * wy2;
		return (v11 + v12 + v21 + v22);
	}

	const T_type& operator()(uint i, uint j, uint k) const
	{
		assert(T_dim == 3);
		Vec<T_dim,uint> idx;
		idx[0] = i;
		idx[1] = j;
		idx[2] = k;
		return At(idx);
	}

	T_type& operator()(uint i, uint j, uint k)
	{
		assert(T_dim == 3);
		Vec<T_dim,uint> idx;
		idx[0] = i;
		idx[1] = j;
		idx[2] = k;
		return At(idx);
	}

	/// Array size as a tuple
	Vec<T_dim,uint> GetSize() const
	{
		return (m_size);
	}

	/// empty?
	bool Empty() const
	{
		return (m_data.Empty());
	}

	/// Erase
	void Erase()
	{
		m_data.Erase();
		m_size = 0;
	}

	/// Raw pointer
	const T_type *Raw() const {return (m_data.Raw());}
	T_type       *Raw()       {return (m_data.Raw());}

	/// Array size of data
	uint sizeOfData() const
	{
		return (m_data.sizeOfData());
	}

protected:
	Vec<T_dim, uint> m_size; // resolution
	CGrid<T_type> m_data;
};

typedef CGridND<2,float> CGrid2Df;
typedef CGridND<3,float> CGrid3Df;
typedef CGridND<2,Vec3f> CGrid2D3f;
typedef CGridND<3,Vec3f> CGrid3D3f;

#endif // GRIDND_H

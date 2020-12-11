#ifndef GRID_H
#define GRID_H

#include "main.h"
#include "vec.h"

template <class T_type>
class CGrid
{
public:
    CGrid()
    {
        m_size = 0;
        m_allocSize = 0;
        m_data = NULL;
    }

    CGrid(uint size)
    {
        m_size = 0;
        m_allocSize = 0;
        m_data = NULL;
        Allocate(size);
    }

    CGrid(const std::vector<T_type>& v)
    {
        m_size = 0;
        m_allocSize = 0;
        m_data = NULL;
        Allocate(v.size());
        for (uint n = 0; n < v.size(); n++)
        {
            m_data[n] = v[n];
        }
    }

    ~CGrid()
    {
        DELETE_ARRAY(m_data);
    }

    /// Erase the array
    void Erase()
    {
        if (m_data != NULL)
        {
            delete[](m_data);
            m_data = NULL;
            m_allocSize = 0;
            m_size = 0;
        }
    }

    /// Allocate the array
    virtual void Allocate(uint size_to_allocate)
    {
        Erase();
        m_allocSize = size_to_allocate;
        m_size = size_to_allocate;
        m_data = new T_type[m_allocSize];
        // init array
        for (uint n = 0; n < int(m_allocSize); n++)
        {
            m_data[n] = 0;
        }
    }

    /// Copy (constructor)
    CGrid(const CGrid& a)
    {
        m_size = 0;
        m_allocSize = 0;
        m_data = NULL;
        if (a.size() > 0)
        {
            allocate(a.size());
            for (uint n = 0; n < a.size(); n++)
            {
                m_data[n] = a[n];
            }
        }
    }

    /// Copy (affectation)
    const CGrid& operator=(const CGrid& a)
    {
        if (a.Size() > m_allocSize) // do not erase if size is sufficient
        {
            Erase();
        }
        if (!a.Empty()) // skip if 'a' is empty
        {
            if (Empty()) // only allocate if necessary
            {
                Allocate(a.Size());
            }
            else
            {
                m_size = a.Size(); // make size equal, keep allocated memory
            }
            for (uint n = 0; n < m_size; n++)
            {
                m_data[n] = a[n];
            }
        }
        return (*this);
    }

    /// Fill array with a given value
    void Fill(const T_type& value_to_fill_with)
    {
        for (uint n = 0; n < m_size; n++)
        {
            m_data[n] = value_to_fill_with;
        }
    }

    // Fill array with a given array
    void Fill(const T_type* data_to_fill_with)
    {
        for (uint n = 0; n < m_size; n++)
        {
            m_data[n] = data_to_fill_with[n];
        }
    }

    /// Resize array
    void Truncate(uint new_size)
    {
        m_size = new_size;
    }

    /// Read only access
    const T_type& operator[](uint n) const
    {
        return (m_data[n]);
    }

    /// Read/write access
    T_type& operator[](uint n)
    {
        return (m_data[n]);
    }

    /// Grid size
    uint Size() const
    {
        return (m_size);
    }

    /// Grid allocated size
    uint AllocatedSize() const
    {
        return (m_allocSize);
    }

    /// Empty?
    bool Empty() const
    {
        return (m_size == 0);
    }

    /// Raw pointer
    const T_type* Raw() const { return (m_data); }
    T_type* Raw() { return (m_data); }

private:
    uint m_size;
    uint m_allocSize;
    T_type* m_data;
};

#endif // GRID_H

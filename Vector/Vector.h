#pragma once
 
#include <assert.h>
#include <vector>
#include <memory>
#include <algorithm>


namespace MaxLib {
namespace Vector {

// Creates a new vector based on items within an input vector
//  Usage:
//       std::vector<typeB> Bs = VectorCopy<typeA, typeB>(As, [](typeA& from) {
//           return typeB(from.a);
//       });
template<typename T1, typename T2>
std::vector<T2> VectorCopy(std::vector<T1>& from, std::function<T2(T1&)> cb = [](T1& item) { return item; })
{
    std::vector<T2> returnItems;
    for(T1& item : from) {
        returnItems.emplace_back(cb(item));
    }
    return returnItems;
}

template<typename T1, typename T2>
std::vector<T2> VectorCopy(const std::vector<T1>& from, std::function<T2(const T1&)> cb = [](const T1& item) { return item; })
{
    std::vector<T2> returnItems;
    for(const T1& item : from) {
        returnItems.emplace_back(cb(item));
    }
    return returnItems;
}


// Don't use this directly, use Vector_Ptrs, Vector_SharedPtrs etc. instead
// Smart_Ptr is either unique_ptr or shared_ptr
template <typename T, template <typename> class Smart_Ptr>
class Vector_SmartPtrs {
public:

    // Constructor with no parameters
    Vector_SmartPtrs() = default;
    
    // Constructor with initialiser list
    Vector_SmartPtrs(std::initializer_list<T> items)
    {
        m_Items.reserve(items.size());
        // fill `m_Items` from the initializer list by creating a smart ptr from each element
        std::transform(items.begin(), items.end(), std::back_inserter(m_Items), [](const T& item) {
            if constexpr (std::is_same_v<Smart_Ptr<T>, std::unique_ptr<T>>) {
                return make_unique<T>(item);
            } else if constexpr (std::is_same_v<Smart_Ptr<T>, std::shared_ptr<T>>) {
                return make_shared<T>(item);
            } else {
                assert(0 && "Unknown Type");
            }
        });
    };

    // Marked virtual so that any derived classes' destructor also gets called 
    virtual ~Vector_SmartPtrs() = default;
    
	// Adds a Polymorphic Item (and returns a pointer to it)
	// usage: T* ptr = v.Addp<Child>()
    template <typename U, typename... Args>
    T* Addp(Args&&... args) {
        // Forward args to make_unique
        if constexpr (std::is_same_v<Smart_Ptr<T>, std::unique_ptr<T>>) {
            m_Items.emplace_back(make_unique<U>(std::forward<Args>(args)...));
        } else if constexpr (std::is_same_v<Smart_Ptr<T>, std::shared_ptr<T>>) {
            m_Items.emplace_back(make_shared<U>(std::forward<Args>(args)...));
        } else {
            assert(0 && "Unknown Type");
        }
        // Return a reference
        return static_cast<U*>(m_Items.back().get());
    }
    
	// Adds an Item (and returns a pointer to it)
	// usage: T* ptr = v.Addp()     (equivelent to v.Addp<Parent>())
    template <typename... Args>
    T* Addp(Args&&... args) {
        // Forward to Add<U>
        return Addp<T>(std::forward<Args>(args)...);
    }
    
	// Adds a preconstructed Item (and returns a pointer to it)
    T* Addp(Smart_Ptr<T> item) {
        // pass item to container
        m_Items.emplace_back(std::move(item));
        // Return a reference
        return static_cast<T*>(m_Items.back().get());
    }
    
	// Adds a Polymorphic Item (and returns a reference to it)
	// usage: T& ref = v.Add<Child>()
    template <typename U, typename... Args>
    T& Add(Args&&... args) {
        // Forward args to make_unique & return a reference
        return static_cast<U&>(*Addp<U>(std::forward<Args>(args)...));
    }
    
	// Adds an Item (and returns a reference to it)
	// usage: T& ref = v.Add()  (equivelent to v.Add<Parent>())
    template <typename... Args>
    T& Add(Args&&... args) {
        // Forward to Add<U>
        return Add<T>(std::forward<Args>(args)...);
    }

	// Adds a preconstructed Item (and returns a reference to it)
    T& Add(Smart_Ptr<T> item) {
        // pass item to make_unique
        return static_cast<T&>(*Addp(std::move(item)));
    }
    
    // Remove item from vector
    void Remove(size_t index) { 
        Assert_IsValid(index);
        m_Items.erase(std::next(m_Items.begin(), index));
    }

    
    // Access last item in vector
    T& Back()                                           { return (*this)[Size() - 1]; } // non-const
    const T& Back() const                               { return (*this)[Size() - 1]; } // const
    
    // Access item in vector
    T& operator[](size_t index)                         { Assert_IsValid(index); return *m_Items[index]; } // non-const
    const T& operator[](size_t index) const             { Assert_IsValid(index); return *m_Items[index]; } // const
    
    template<typename U>
    U* CastItem(size_t index)                           { Assert_IsValid(index); return dynamic_cast<U*>(m_Items[index].get()); } // non-const
    template<typename U>
    const U* CastItem(size_t index) const               { Assert_IsValid(index); return dynamic_cast<U*>(m_Items[index].get()); } // const
    
    // Copy smart pointer
    Smart_Ptr<T> CopyItem(size_t index)                 { Assert_IsValid(index); 
        
        if constexpr (std::is_same_v<Smart_Ptr<T>, std::shared_ptr<T>>) {
            return m_Items[index]; 
        } else {
            assert(0 && "Type must be shared_ptr");
        }
    } 

    // Swaps items[n] and items[n_Next]
    void SwapItems(size_t n, size_t n_Next) {
        Assert_IsValid(n); 
        Assert_IsValid(n_Next);
        std::swap(m_Items[n], m_Items[n_Next]);
    }
    
    // Gets the number of elements in the vector
    size_t Size() const                         { return m_Items.size(); }
    // Gets the number of uses of shared pointer
    size_t UseCount(size_t index) const         { 
        Assert_IsValid(index); 
        if constexpr (std::is_same_v<Smart_Ptr<T>, std::unique_ptr<T>>) {
            return 1;
        } else if constexpr (std::is_same_v<Smart_Ptr<T>, std::shared_ptr<T>>) {
            return m_Items[index].use_count();
        } else {
            assert(0 && "Unknown Type");
        }
    }
    
    // Enable range base loop
    auto begin()         { return m_Items.begin(); }
    auto begin() const   { return m_Items.begin(); }
    auto end()           { return m_Items.end(); }
    auto end() const     { return m_Items.end(); }
        
protected:
    // The container
    std::vector<Smart_Ptr<T>> m_Items;
    
    // Validity function ( 0 < x < size)
    void Assert_IsValid(int index)const        { assert(index > -1 && index < (int)Size()); }
};


// Don't use this directly, use Vector_Ptrs, Vector_SharedPtrs etc. instead
// Smart_Ptr is either unique_ptr or shared_ptr
template <typename T, template <typename> class Smart_Ptr>
class Vector_SelectableSmartPtrs : public Vector_SmartPtrs<T, Smart_Ptr>
{
public:

    // Constructor with no parameters
    Vector_SelectableSmartPtrs() = default;
    
    // Constructor with initialiser list
    Vector_SelectableSmartPtrs(std::initializer_list<T> items) 
        : m_CurrentIndex(items.size() - 1), Vector_SmartPtrs<T, Smart_Ptr>(items) {}
    
	// Adds a Polymorphic Item (and returns a reference to it)
	// usage: v.Add<sub_class>()
    template <typename U, typename... Args>
    T& Add(Args&&... args) {
        // Set current item to last in list
        m_CurrentIndex = Vector_SmartPtrs<T, Smart_Ptr>::m_Items.size();
        // Forward args to base class Add<U> & return a reference
        return Vector_SmartPtrs<T, Smart_Ptr>::template Add<U>(std::forward<Args>(args)...);
    }
    
	// Adds an Item (and returns a reference to it)
	// usage: v.Add() (equivelent to v.Add<base_class>())
    template <typename... Args>
    T& Add(Args&&... args) {
        // Forward to Add<U>
        return Add<T>(std::forward<Args>(args)...);
    }

    // Remove item from vector
    void Remove(size_t index) {
        Vector_SmartPtrs<T, Smart_Ptr>::Remove(index);
        if (m_CurrentIndex != -1 && (int)index <= m_CurrentIndex)
            m_CurrentIndex--;
    }
    
    // Access current item in vector
    T& CurrentItem()                            { return *CastCurrentItem<T>(); } // non-const
    const T& CurrentItem() const                { return *CastCurrentItem<T>(); } // const
    
    // Retrieves item in vector and converts it to derived class type: CurrentItem<DerivedClass>()
    template<typename U>
    U* CastCurrentItem()                        { Vector_SmartPtrs<T, Smart_Ptr>::Assert_IsValid(m_CurrentIndex); return dynamic_cast<U*>(Vector_SmartPtrs<T, Smart_Ptr>::m_Items[m_CurrentIndex].get()); } // non-const
    template<typename U>
    const U* CastCurrentItem() const            { Vector_SmartPtrs<T, Smart_Ptr>::Assert_IsValid(m_CurrentIndex); return dynamic_cast<U*>(Vector_SmartPtrs<T, Smart_Ptr>::m_Items[m_CurrentIndex].get()); } // const
    
    // Set the current item by index
    void SetCurrentIndex(int index)             { Assert_IsValidOrUnset(index); m_CurrentIndex = index; }
    
    // Remove the current item
    void RemoveCurrent()                        { Vector_SmartPtrs<T, Smart_Ptr>::Assert_IsValid(m_CurrentIndex); Remove(m_CurrentIndex); }
    
    // Checks if there is a current item
    bool HasItemSelected()                      { return m_CurrentIndex != -1; }
    
    // Get the current item index
    int CurrentIndex() const                    { return m_CurrentIndex; }

private:
    // Current item index
    int m_CurrentIndex = -1;
    
    // Validity functions 
    void Assert_IsValidOrUnset(int index)       { assert(index >= -1 && (index < (int)Vector_SmartPtrs<T, Smart_Ptr>::Size())); }  
};


// A vector wrapper for unique ptrs, simplifies polymorphism
// Note: The destructor of the Base Class must be marked virtual, otherwise the derived classes' constructor won't get called too
/*  Usage:
        Vector_Ptrs<Shape> v;
        
        v.Add("Shape");
        v.Add<Point>("Point", 10.0f, 20.0f );
        
        Shape& s1 = v[0]; 
                
        if(Point* p1 = v.CastItem<Point>(1)) {
            // Success
        }
*/
template <typename T>
using Vector_Ptrs = Vector_SmartPtrs<T, std::unique_ptr>;
template <typename T>
using Vector_SharedPtrs = Vector_SmartPtrs<T, std::shared_ptr>;

// A vector wrapper for unique ptrs, simplifies polymorphism and enables selecting a current item
// Note: The destructor of the Base Class must be marked virtual, otherwise the derived classes' constructor won't get called too
/*  Usage:
        Vector_SelectablePtrs<Shape> v;
        
        v.Add("Shape");
        v.Add<Point>("Point", 10.0f, 20.0f );
        
        Shape& s1 = v[0]; 
                
        if(Point* p1 = v.CastItem<Point>(1)) {
            // Success
        }
         
        Shape& s2 = v.CurrentItem();
        
        if(Point* p2 = v.CastCurrentItem<Point>()) {
            // Success
        }
*/
template <typename T>
using Vector_SelectablePtrs = Vector_SelectableSmartPtrs<T, std::unique_ptr>;
template <typename T>
using Vector_SelectableSharedPtrs = Vector_SelectableSmartPtrs<T, std::shared_ptr>;



} // end namespace Vector
} // end namespace MaxLib

#pragma once
 
#include <assert.h>
#include <functional>
#include <vector>
#include <memory>
#include <algorithm>
#include "Vector_Impl.h"


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
class Vector_Ptrs : public Vector_SmartPtrs<T, std::unique_ptr>
{
    // Constructor from child
    public: using Vector_SmartPtrs<T, std::unique_ptr>::Vector_SmartPtrs;
}; 
 
template <typename T>
class Vector_SharedPtrs : public Vector_SmartPtrs<T, std::shared_ptr>
{
    // Constructor from child
    public: using Vector_SmartPtrs<T, std::shared_ptr>::Vector_SmartPtrs;
}; 
 
 
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
class Vector_SelectablePtrs : public Vector_SelectableSmartPtrs<T, std::unique_ptr>
{
    // Constructor from child
    public: using Vector_SelectableSmartPtrs<T, std::unique_ptr>::Vector_SelectableSmartPtrs;
}; 

template <typename T>
class Vector_SelectableSharedPtrs : public Vector_SelectableSmartPtrs<T, std::shared_ptr>
{
    // Constructor from child
    public: using Vector_SelectableSmartPtrs<T, std::shared_ptr>::Vector_SelectableSmartPtrs;
}; 
 

} // end namespace Vector
} // end namespace MaxLib

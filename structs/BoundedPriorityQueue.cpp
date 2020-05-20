//
// Created by pavel on 23.09.2019.
//
#include "BoundedPriorityQueue.h"

template<class T, class Compare>
size_t knn_queue::BoundedPriorityQueue<T, Compare>::size() const { return elements_.size(); }

template<class T, class Compare>
const T &knn_queue::BoundedPriorityQueue<T, Compare>::operator[](size_t index) const { return elements_[index]; }

template<class T, class Compare>
const T &knn_queue::BoundedPriorityQueue<T, Compare>::back() const { return elements_.back(); }

template<class T, class Compare>
void knn_queue::BoundedPriorityQueue<T, Compare>::push(const T &val) {
    auto it = std::find_if(std::begin(elements_), std::end(elements_),
                           [&](const T &element) { return Compare()(val, element); });
    elements_.insert(it, val);

    if (elements_.size() > bound_)
        elements_.resize(bound_);
}

template<class T, class Compare>
knn_queue::BoundedPriorityQueue<T, Compare>::BoundedPriorityQueue(size_t bound)  : bound_(bound) {
    elements_.reserve(bound + 1);
};

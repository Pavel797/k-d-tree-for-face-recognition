//
// Created by pavel on 23.09.2019.
//

#ifndef FACESDATABASE_BOUNDEDPRIORITYQUEUE_H
#define FACESDATABASE_BOUNDEDPRIORITYQUEUE_H

#include <functional>
namespace knn_queue {
    template<class T, class Compare = std::less<T>>
    class BoundedPriorityQueue {
    public:

        BoundedPriorityQueue() = delete;

        explicit BoundedPriorityQueue(size_t bound);

        void push(const T &val);

        const T &back() const;

        const T &operator[](size_t index) const;

        size_t size() const;

    private:
        size_t bound_;
        std::vector<T> elements_;
    };

}

#endif //FACESDATABASE_BOUNDEDPRIORITYQUEUE_H

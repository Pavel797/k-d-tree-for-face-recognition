//
// Created by pavel on 22.09.2019.
//

#ifndef FACESDATABASE_KDTREE_H
#define FACESDATABASE_KDTREE_H

#include <vector>
#include <numeric>
#include <algorithm>
#include <exception>
#include <functional>
#include <cmath>
#include <stack>

#include "KDNode.h"
#include "FaceEmbedding.h"
#include "BoundedPriorityQueue.h"

namespace kdt {
    class KDTree {
    public:
        explicit KDTree(size_t dim) : root_(nullptr), dim(dim) {};

        explicit KDTree(const std::vector<embedding::FaceEmbedding> &points) : root_(nullptr) {
            assert(!points.empty()); // mast be is not empty

            dim = points[0].length();
            for (auto &i: points) {
                if (i.length() != dim) {
                    throw Exception();
                }
            }

            build(points);
        }

        ~KDTree() { clear(); }

        void build(const std::vector<embedding::FaceEmbedding> &points) {
            clear();

            points_ = points;

            std::vector<int> indices(points.size());
            std::iota(std::begin(indices), std::end(indices), 0);

            root_ = buildRecursive(indices.data(), (int) points.size(), 0);
        }

        void clear() {
            clearRecursive(root_);
            root_ = nullptr;
            points_.clear();
        }

        [[nodiscard]] bool validate() const {
            try {
                validateRecursive(root_, 0);
            }
            catch (const Exception &) {
                return false;
            }

            return true;
        }

        int nnSearch(const embedding::FaceEmbedding &query, double *minDist = nullptr) const {
            int guess;
            double _minDist = std::numeric_limits<double>::max();

            nnSearchRecursive(query, root_, &guess, &_minDist);

            if (minDist)
                *minDist = _minDist;

            return guess;
        }

        [[nodiscard]] std::vector<int> knnSearch(const embedding::FaceEmbedding &query, int k) const {
            KnnQueue queue(k);
            knnSearchRecursive(query, root_, queue, k);

            std::vector<int> indices(queue.size());
            for (size_t i = 0; i < queue.size(); i++)
                indices[i] = queue[i].second;

            return indices;
        }

        [[nodiscard]] std::vector<int> radiusSearch(const embedding::FaceEmbedding &query, double radius) const {
            std::vector<int> indices;
            radiusSearchRecursive(query, root_, indices, radius);
            return indices;
        }

    private:

        class Exception : public std::exception {
            using std::exception::exception;
        };

        using KnnQueue = knn_queue::BoundedPriorityQueue<std::pair<double, int>>;

        struct StDump {
            int *index_arr;
            int npoints;
            size_t depth;

            KDNode *node;
            size_t id_ax;

            StDump(int *index_arr, size_t npoints, size_t depth, KDNode *node, size_t id_ax) :
                    index_arr(index_arr), npoints(npoints), depth(depth), node(node), id_ax(id_ax) {}
        };

        std::pair<int, KDNode *> buildNode(int *index_arr, int depth, int npoints) {
            const size_t axis = depth % dim;
            const int mid = (npoints - 1) / 2;

            std::nth_element(index_arr, index_arr + mid, index_arr + npoints, [&](int lhs, int rhs) {
                return points_[lhs][axis] < points_[rhs][axis];
            });

            auto *node = new KDNode();
            node->idx = index_arr[mid];
            node->axis = axis;

            return std::make_pair(mid, node);
        }

        KDNode *_build(int *index_arr, int nps, size_t dp) {
            std::stack<StDump> st;

            auto *rnode = new KDNode();
            rnode->idx = index_arr[(nps - 1) / 2];
            rnode->axis = dp % dim;

            std::cout << nps << "  " << std::endl;

            st.push(StDump(index_arr, nps, dp, rnode, -1));

            while (true) {
                int st_size = st.size();

                StDump dump = st.top();
                st.pop();

                auto resBuild = buildNode(index_arr, dump.depth, dump.npoints);

                if (dump.depth == st.top().depth) {
                    StDump dump2 = st.top();
                    auto resBuild2 = buildNode(index_arr, dump.depth, dump.npoints);

                    if (resBuild2.first > 0) {
                        st.push(StDump(index_arr, resBuild2.first, dump.depth + 1, resBuild2.second, 0));
                    }
                    if (int(dump.npoints - resBuild2.first - 1) > 0) {
                        st.push(StDump(index_arr + resBuild2.first + 1, dump.npoints - resBuild2.first - 1,
                                       dump.depth + 1, resBuild2.second, 1));
                    }
                }


                st.push(dump);
                if (resBuild.first > 0) {
                    st.push(StDump(index_arr, resBuild.first, dump.depth + 1, resBuild.second, 0));
                }
                if (int(dump.npoints - resBuild.first - 1) > 0) {
                    st.push(StDump(index_arr + resBuild.first + 1, dump.npoints - resBuild.first - 1, dump.depth + 1,
                                   resBuild.second, 1));
                }

                if (st_size == st.size()) {
                    break;
                }
            }

            bool is_t2_not_null = false;
            StDump t1 = st.top(), t2 = st.top();
            st.pop();

            std::cout << "0 " << t1.depth << std::endl;

            if (!st.empty() && t1.depth == st.top().depth) {
                t2 = st.top();
                st.pop();
                std::cout << "1 " << t2.depth << std::endl;
                is_t2_not_null = true;
            }

            while (!st.empty()) {
                st.top().node->next[t1.id_ax] = t1.node;
                if (is_t2_not_null) {
                    st.top().node->next[t2.id_ax] = t2.node;
                }

                t1 = st.top();
                std::cout << "0 " << t1.depth << std::endl;
                st.pop();

                if (!st.empty() && t1.depth == st.top().depth) {
                    t2 = st.top();
                    st.pop();
                    std::cout << "1 " << t2.depth << std::endl;
                    is_t2_not_null = true;
                } else {
                    is_t2_not_null = false;
                }
            }

            return rnode;
        }

        KDNode *buildRecursive(int *indices, size_t npoints, size_t depth) {
            if (npoints <= 0)
                return nullptr;

            const size_t axis = depth % dim;
            const size_t mid = (npoints - 1) / 2;

            std::nth_element(indices, indices + mid, indices + npoints, [&](int lhs, int rhs) {
                return points_[lhs][axis] < points_[rhs][axis];
            });

            auto *node = new KDNode();
            node->idx = indices[mid];
            node->axis = axis;

            node->next[0] = buildRecursive(indices, mid, depth + 1);
            node->next[1] = buildRecursive(indices + mid + 1, npoints - mid - 1, depth + 1);

            return node;
        }

        /** @brief Clears k-d tree recursively.
        */
        static void clearRecursive(KDNode *node) {
            if (node == nullptr)
                return;

            if (node->next[0])
                clearRecursive(node->next[0]);

            if (node->next[1])
                clearRecursive(node->next[1]);

            delete node;
        }

        /** @brief Validates k-d tree recursively.
        */
        void validateRecursive(const KDNode *node, int depth) const {
            if (node == nullptr)
                return;

            const int axis = node->axis;
            const KDNode *node0 = node->next[0];
            const KDNode *node1 = node->next[1];

            if (node0 && node1) {
                if (points_[node->idx][axis] < points_[node0->idx][axis])
                    throw Exception();

                if (points_[node->idx][axis] > points_[node1->idx][axis])
                    throw Exception();
            }

            if (node0)
                validateRecursive(node0, depth + 1);

            if (node1)
                validateRecursive(node1, depth + 1);
        }

        static double distance(const embedding::FaceEmbedding &p, const embedding::FaceEmbedding &q) {
            double dist = 0;
            for (size_t i = 0; i < p.length(); i++)
                dist += (p[i] - q[i]) * (p[i] - q[i]);
            return sqrt(dist);
        }

        struct StackSearchDump {
            const embedding::FaceEmbedding &query;
            const KDNode *node;
            int *guess;
            double *minDist;

            StackSearchDump(const embedding::FaceEmbedding &query, const KDNode *node, int *guess,
                            double *minDist) : query(query), node(node), guess(guess), minDist(minDist) {}
        };

        void nnSearchRecursive(const embedding::FaceEmbedding &query, const KDNode *node, int *guess,
                               double *minDist) const {
            if (node == nullptr)
                return;

            const embedding::FaceEmbedding &train = points_[node->idx];

            const double dist = distance(query, train);
            if (dist < *minDist) {
                *minDist = dist;
                *guess = node->idx;
            }

            const int axis = node->axis;
            const int dir = query[axis] < train[axis] ? 0 : 1;
            nnSearchRecursive(query, node->next[dir], guess, minDist);

            const double diff = fabs(query[axis] - train[axis]);
            if (diff < *minDist)
                nnSearchRecursive(query, node->next[!dir], guess, minDist);
        }

        void nnSearchLiner(const embedding::FaceEmbedding &_query, const KDNode *_node, int *_guess,
                           double *_minDist) const {
            std::stack<StackSearchDump> st;
            st.push(StackSearchDump(_query, _node, _guess, _minDist));

            while (!st.empty()) {
                StackSearchDump dump = st.top();
                st.pop();

                const embedding::FaceEmbedding &train = points_[dump.node->idx];

                const double dist = distance(dump.query, train);
                if (dist < *dump.minDist) {
                    *dump.minDist = dist;
                    *dump.guess = dump.node->idx;
                }

                const int axis = dump.node->axis;
                const int dir = dump.query[axis] < train[axis] ? 0 : 1;
                if (dump.node->next[dir] != nullptr) {
                    st.push(StackSearchDump(dump.query, dump.node->next[dir], dump.guess, dump.minDist));
                }

                const double diff = fabs(dump.query[axis] - train[axis]);
                if (diff < *dump.minDist && dump.node->next[!dir] != nullptr)
                    st.push(StackSearchDump(dump.query, dump.node->next[!dir], dump.guess, dump.minDist));
            }
        }

        void
        knnSearchRecursive(const embedding::FaceEmbedding &query, const KDNode *node, KnnQueue &queue, int k) const {
            if (node == nullptr)
                return;

            const embedding::FaceEmbedding &train = points_[node->idx];

            const double dist = distance(query, train);
            queue.push(std::make_pair(dist, node->idx));

            const int axis = node->axis;
            const int dir = query[axis] < train[axis] ? 0 : 1;
            knnSearchRecursive(query, node->next[dir], queue, k);

            const double diff = fabs(query[axis] - train[axis]);
            if ((int) queue.size() < k || diff < queue.back().first)
                knnSearchRecursive(query, node->next[!dir], queue, k);
        }

        void radiusSearchRecursive(const embedding::FaceEmbedding &query, const KDNode *node, std::vector<int> &indices,
                                   double radius) const {
            if (node == nullptr)
                return;

            const embedding::FaceEmbedding &train = points_[node->idx];

            const double dist = distance(query, train);
            if (dist < radius)
                indices.push_back(node->idx);

            const int axis = node->axis;
            const int dir = query[axis] < train[axis] ? 0 : 1;
            radiusSearchRecursive(query, node->next[dir], indices, radius);

            const double diff = fabs(query[axis] - train[axis]);
            if (diff < radius)
                radiusSearchRecursive(query, node->next[!dir], indices, radius);
        }

        KDNode *root_;
        std::vector<embedding::FaceEmbedding> points_;
        size_t dim;
    };
} // kdt

#endif //FACESDATABASE_KDTREE_H

//
// Created by pavel on 23.09.2019.
//

#ifndef FACESDATABASE_KDNODE_H
#define FACESDATABASE_KDNODE_H

#include <cstddef>

struct KDNode
{
    int idx;
    KDNode* next[2];
    size_t axis;

    KDNode();
};

#endif //FACESDATABASE_KDNODE_H

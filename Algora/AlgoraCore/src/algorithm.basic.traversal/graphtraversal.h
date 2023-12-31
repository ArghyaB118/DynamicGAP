/**
 * Copyright (C) 2013 - 2019 : Kathrin Hanauer
 *
 * This file is part of Algora.
 *
 * Algora is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Algora is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Algora.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact information:
 *   http://algora.xaikal.org
 */

#ifndef GRAPHTRAVERSAL_H
#define GRAPHTRAVERSAL_H

#include "algorithm/propertycomputingalgorithm.h"
#include "graph/digraph.h"
#include "graph/vertex.h"
#include "graph/graph_functional.h"

namespace Algora {

template<typename PropertyType, bool reverseArcDirection = false, bool ignoreArcDirection = false>
class GraphTraversal : public PropertyComputingAlgorithm<DiGraph::size_type, PropertyType>
{
public:
    explicit GraphTraversal(bool computeValues = true)
        : PropertyComputingAlgorithm<DiGraph::size_type, PropertyType>(computeValues),
          startVertex(nullptr),
          onVertexDiscovered(vertexTrue), onArcDiscovered(arcTrue),
          vertexStopCondition(vertexFalse), arcStopCondition(arcFalse)
    { }

    virtual ~GraphTraversal() { }

    void setStartVertex(const Vertex *v) {
        startVertex = v;
    }

    void onVertexDiscover(const VertexPredicate &vFun) {
        onVertexDiscovered = vFun;
    }

    void onArcDiscover(const ArcPredicate &aFun) {
        onArcDiscovered = aFun;
    }

    void setVertexStopCondition(const VertexPredicate &vStop) {
        vertexStopCondition = vStop;
    }

    void setArcStopCondition(const ArcPredicate &aStop) {
        arcStopCondition = aStop;
    }

    virtual DiGraph::size_type numVerticesReached() const = 0;

    // DiGraphAlgorithm interface
public:
    virtual bool prepare() override
    {
        return PropertyComputingAlgorithm<DiGraph::size_type, PropertyType>::prepare()
                && ( startVertex == nullptr
                     || (this->diGraph->containsVertex(startVertex) && startVertex->isValid()));
    }


protected:
    const Vertex *startVertex;

    VertexPredicate onVertexDiscovered;
    ArcPredicate onArcDiscovered;
    VertexPredicate vertexStopCondition;
    ArcPredicate arcStopCondition;
};

}

#endif // GRAPHTRAVERSAL_H

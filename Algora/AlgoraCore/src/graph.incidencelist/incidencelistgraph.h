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

#ifndef INCIDENCELISTGRAPH_H
#define INCIDENCELISTGRAPH_H

#include "graph/digraph.h"

namespace Algora {

class IncidenceListVertex;
class IncidenceListGraphImplementation;
template<typename T>
class ModifiableProperty;

class IncidenceListGraph : public DiGraph
{
public:
    explicit IncidenceListGraph(GraphArtifact *parent = nullptr);
    virtual ~IncidenceListGraph() override;

    // copying
    IncidenceListGraph(const IncidenceListGraph &other,
                                     ModifiableProperty<GraphArtifact*> *otherToThisVertices = nullptr,
                                     ModifiableProperty<GraphArtifact*> *otherToThisArcs = nullptr,
                                     ModifiableProperty<GraphArtifact*> *thisToOtherVertices = nullptr,
                                     ModifiableProperty<GraphArtifact*> *thisToOtherArcs = nullptr
                       );
    IncidenceListGraph &operator=(const IncidenceListGraph &other);
    IncidenceListGraph &assign(const IncidenceListGraph &other,
                                     ModifiableProperty<GraphArtifact*> *otherToThisVertices = nullptr,
                                     ModifiableProperty<GraphArtifact*> *otherToThisArcs = nullptr,
                                     ModifiableProperty<GraphArtifact*> *thisToOtherVertices = nullptr,
                                     ModifiableProperty<GraphArtifact*> *thisToOtherArcs = nullptr
                               );

    // moving
    IncidenceListGraph(IncidenceListGraph &&other);
    IncidenceListGraph &operator=(IncidenceListGraph &&other);

    // Graph interface
public:
    virtual Vertex *addVertex() override;
    virtual void removeVertex(Vertex *v) override;
    virtual bool containsVertex(const Vertex *v) const override;
    virtual Vertex *getAnyVertex() const override;
    virtual IncidenceListVertex *vertexAt(size_type i) const;

    virtual void mapVerticesUntil(const VertexMapping &vvFun, const VertexPredicate &breakCondition) override;

    virtual bool isEmpty() const override;
    virtual size_type getSize() const override;

    virtual void clear() override;
    virtual void clearAndRelease();
    virtual void clearOrderedly();

    // DiGraph interface
public:
    DiGraph *createReversedGraph(PropertyMap<GraphArtifact *> &map) const override;
    virtual Arc *addArc(Vertex *tail, Vertex *head) override;
    virtual MultiArc *addMultiArc(Vertex *tail, Vertex *head, size_type size) override;
    virtual void removeArc(Arc *a) override;
    virtual bool containsArc(const Arc *a) const override;
    virtual Arc *findArc(const Vertex *from, const Vertex *to) const override;
    virtual size_type getNumArcs(bool multiArcsAsSimple) const override;

    virtual size_type getOutDegree(const Vertex *v, bool multiArcsAsSimple) const override;
    virtual size_type getInDegree(const Vertex *v, bool multiArcsAsSimple) const override;
    virtual bool isSource(const Vertex *v) const override;
    virtual bool isSink(const Vertex *v) const override;

    virtual void mapArcsUntil(const ArcMapping &avFun, const ArcPredicate &breakCondition) override;
    virtual void mapOutgoingArcsUntil(const Vertex *v, const ArcMapping &avFun, const ArcPredicate &breakCondition) override;
    virtual void mapIncomingArcsUntil(const Vertex *v, const ArcMapping &avFun, const ArcPredicate &breakCondition) override;

public:
    void bundleParallelArcs();
    void unbundleParallelArcs();

    void reserveVertexCapacity(size_type n);
    void reserveArcCapacity(size_type n);

    void activateVertex(Vertex *v, bool activateIncidentArcs);
    void deactivateVertex(Vertex *v);
    void activateArc(Arc *a);
    void deactivateArc(Arc *a);

protected:
    IncidenceListVertex *recycleOrCreateIncidenceListVertex();
    IncidenceListVertex *createIncidenceListVertex();
    Arc *recycleOrCreateArc(IncidenceListVertex *tail, IncidenceListVertex *head);

    // either un-hide function of base class or leave it hidden and disable clang warning
    // using DiGraph::createArc;
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#endif
    Arc *createArc(IncidenceListVertex *tail, IncidenceListVertex *head);
#ifdef __clang__
#pragma clang diagnostic pop
#endif

private:
    IncidenceListGraphImplementation *impl;
};

}

#endif // INCIDENCELISTGRAPH_H

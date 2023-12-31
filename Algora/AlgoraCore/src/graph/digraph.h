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

#ifndef DIGRAPH_H
#define DIGRAPH_H

#include "graph.h"
#include "arc.h"
#include "weightedarc.h"
#include "parallelarcsbundle.h"

#include "graph.visitor/arcvisitor.h"

namespace Algora {

class Vertex;
template <typename T>
class PropertyMap;

class DiGraph : public Graph
{
public:
    explicit DiGraph(GraphArtifact *parent = nullptr);
    virtual ~DiGraph() override { }

    // copying
    DiGraph(const DiGraph &other);
    DiGraph &operator=(const DiGraph &other);

    // moving
    DiGraph(DiGraph &&other) = default;
    DiGraph &operator=(DiGraph &&other) = default;

    virtual DiGraph *createReversedGraph(PropertyMap<GraphArtifact*> &) const { return nullptr; }

    // Arcs
    virtual Arc *addArc(Vertex *tail, Vertex *head) = 0;
    virtual MultiArc *addMultiArc(Vertex *tail, Vertex *head, size_type size) = 0;
    virtual void removeArc(Arc *a) = 0;
    virtual bool containsArc(const Arc *a) const = 0;
    virtual Arc *findArc(const Vertex *from, const Vertex *to) const = 0;

    virtual size_type getOutDegree(const Vertex *v, bool multiArcsAsSimple) const = 0;
    virtual size_type getInDegree(const Vertex *v, bool multiArcsAsSimple) const = 0;
    virtual size_type getDegree(const Vertex *v, bool multiArcsAsSimple) const {
        return getOutDegree(v, multiArcsAsSimple) + getInDegree(v, multiArcsAsSimple);
    }
    virtual bool isSource(const Vertex *v) const {
        return getInDegree(v, true) == 0;
    }
    virtual bool isSink(const Vertex *v) const {
        return getOutDegree(v, true) == 0;
    }
    virtual bool isIsolated(const Vertex *v) const {
        return isSource(v) && isSink(v);
    }
    virtual size_type getNumArcs(bool multiArcsAsSimple) const;

    virtual void onArcAdd(void *id, const ArcMapping &avFun);
    virtual void onArcRemove(void *id, const ArcMapping &avFun);
    virtual void removeOnArcAdd(void *id);
    virtual void removeOnArcRemove(void *id);

    // Accomodate visitors
    virtual void acceptArcVisitor(ArcVisitor *aVisitor) {
        mapArcs(aVisitor->getVisitorFunction());
    }
    virtual void acceptOutgoingArcVisitor(const Vertex *v, ArcVisitor *aVisitor) {
        mapOutgoingArcs(v, aVisitor->getVisitorFunction());
    }
    virtual void acceptIncomingArcVisitor(const Vertex *v, ArcVisitor *aVisitor) {
        mapIncomingArcs(v, aVisitor->getVisitorFunction());
    }
    virtual void mapArcs(const ArcMapping &avFun) {
        mapArcsUntil(avFun, arcFalse);
    }
    virtual void mapOutgoingArcs(const Vertex *v, const ArcMapping &avFun) { mapOutgoingArcsUntil(v, avFun, arcFalse); }
    virtual void mapIncomingArcs(const Vertex *v, const ArcMapping &avFun) { mapIncomingArcsUntil(v, avFun, arcFalse); }
    virtual void mapIncidentArcs(const Vertex *v, const ArcMapping &avFun) { mapOutgoingArcs(v, avFun); mapIncomingArcs(v, avFun); }

    virtual void mapArcsUntil(const ArcMapping &avFun, const ArcPredicate &breakCondition) = 0;
    virtual void mapOutgoingArcsUntil(const Vertex *v, const ArcMapping &avFun, const ArcPredicate &breakCondition) = 0;
    virtual void mapIncomingArcsUntil(const Vertex *v, const ArcMapping &avFun, const ArcPredicate &breakCondition) = 0;

    virtual void clear() override;

    // GraphArtifact interface
public:
    virtual std::string toString() const override;

protected:
   Observable<Arc*> observableArcGreetings;
   Observable<Arc*> observableArcFarewells;

   void greetArc(Arc *a) { observableArcGreetings.notifyObservers(a); }
   void dismissArc(Arc *a) { observableArcFarewells.notifyObservers(a); }

    virtual Arc *createArc(Vertex *tail, Vertex *head) {
        return new Arc(tail, head, this);
    }

    virtual Arc *createArc(Vertex *tail, Vertex *head, id_type id) {
        return new Arc(tail, head, id, this);
    }

    virtual MultiArc *createMultiArc(Vertex *tail, Vertex *head, size_type size) {
        return new WeightedArc(tail, head, size, this);
    }

    virtual MultiArc *createMultiArc(Vertex *tail, Vertex *head, size_type size, id_type id) {
        return new WeightedArc(tail, head, size, id, this);
    }

    void invalidateArc(Arc *a) {
        a->invalidate();
    }
};

}

#endif // DIGRAPH_H

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from toverboom.intersections import *
import numpy as np
import functools
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import collections
import random
import itertools

def interpolateBezier( points, steps=10, t=None):
    points = tuple(points)
    if len(points)==3:
        mapper = lambda t,p: (1-t)**2 * p[0] + 2*(1-t)*t*p[1] + t**2*p[2]
    elif len(points)==4:
        mapper = lambda t,p: (np.power( (1-t),3)*p[0] +\
         3* np.power((1-t),2) *t *p[1] +\
         3*(1-t)*np.power(t,2)*p[2] +\
         np.power(t,3)*p[3])
    else:
        raise Exception('Can only interpolate cubic and quadratic splines (3 or 4 parameters, got: %s'  % str(points))

    if t is not None:
        return   mapper(t, [q[0] for q in points]), mapper(t, [q[1] for q in points])
    xGen = ( mapper(t, [q[0] for q in points]) for t in np.linspace(0, 1, steps) )
    yGen = ( mapper(t, [q[1] for q in points]) for t in np.linspace(0, 1, steps) )

    return zip(xGen, yGen)

def interpolateBezierAngle(points, t, ds=0.001):
    x0, y0 = interpolateBezier(points, t=t-ds)
    x1, y1 = interpolateBezier(points, t=t+ds)

    angs = np.arctan2( y1-y0, x0-x1)
    if type(angs)==np.float64:
        if t==0 or t==1:
            angs = np.pi
    else:
        angs[ (t==0) | (t==1) ] = np.pi
    return angs




def format_x_axis_labels(ax, divider=10):
    def d(x, pos):
        return str(int(x/divider))
    formatter = FuncFormatter(d)
    ax.xaxis.set_major_formatter(formatter)


""" Remove spines from axis """
def despine(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False)

class LineageGraph():

    def __init__(self, graph=None):
        self.graph = graph
        self.radiusAttribute = None
        self.defaultRadius=None
        self.xDistance=1
        self.radiusMultiplier=1
        self.lineageSizes=None
        self.lineageExistences=None
        self.trellisOrder=None
        self.verticalSpacing=None
        self.minRadius=0.1
        self.initRadius=0.2


    def getInitRadius(self, override=None):
        if override is not None:
            return override
        if self.minRadius is not None:
            return self.initRadius
        else:
            raise ValueError('No initRadius set')

    def getMinRadius(self, override=None):
        if override is not None:
            return override
        if self.minRadius is not None:
            return self.minRadius
        else:
            raise ValueError('No minRadius set')

    def getTrellisOrder(self, override=None):
        if override is not None:
            return override
        if self.trellisOrder is not None:
            return self.trellisOrder
        else:
            raise ValueError('No trellisOrder set')

    def getLineageSizes(self, override=None):
        if override is not None:
            return override
        if self.lineageSizes is not None:
            return self.lineageSizes
        else:
            raise ValueError('No lineageSizes set')

    def getVerticalSpacing(self, override=None):
        if override is not None:
            return override
        if self.verticalSpacing is not None:
            return self.verticalSpacing
        else:
            raise ValueError('No vertical spacing set')

    def getRadiusAttribute(self, override=None):
        if override is not None:
            return override
        if self.radiusAttribute is not None:
            return self.radiusAttribute
        else:
            raise ValueError('No radiusAttribute set')

    def getXDistance(self, override=None):
        if override is not None:
            return override
        if self.xDistance is not None:
            return self.xDistance
        else:
            raise ValueError('No xDistance set')

    def getDefaultRadius(self, override=None):
        if override is not None:
            return override
        if self.defaultRadius is not None:
            return self.defaultRadius
        else:
            raise ValueError('No defaultRadius set')

    def getGraph(self, override=None):
        if override is not None:
            return override
        if self.graph is not None:
            return self.graph
        else:
            raise ValueError('No graph set')

    def getRadiusMultiplier(self, override=None):
        if override is not None:
            return override
        if self.radiusMultiplier is not None:
            return self.radiusMultiplier
        else:
            raise ValueError('No radiusMultiplier set')


    def getNodeRadius(self, node,**kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        minRadius = self.getMinRadius(kwargs.get('minRadius'))
        radiusAttribute = self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))

        return max(minRadius, graph.nodes[node].get(radiusAttribute, defaultRadius)*radiusMultiplier )

    """Obtain the two radii for an edge from source to sink"""
    def getEdgeRadii(self, source, sink, **kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        minRadius = self.getMinRadius(kwargs.get('minRadius'))
        radiusAttribute = self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        initRadius = self.getInitRadius(kwargs.get('initRadius'))

        r1 =  self.getNodeRadius(sink, **kwargs)

        if source[0]==sink[0]:
            r0 =  self.getNodeRadius(source, **kwargs)
        else:
            r0 = initRadius

        return r0,r1

    def getSizePerLineage(self, **kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        radiusAttribute = self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        minRadius = self.getMinRadius(kwargs.get('minRadius'))

        lineageSizes = {}
        for (lineage, timePoint) in graph:
            current = self.getNodeRadius(
                                         (lineage, timePoint),
                                        **kwargs)
            lineageSizes[lineage] = max( lineageSizes.get(lineage,0), current)
        return lineageSizes

    def getLineageExistences(self,**kwargs):
        graph = self.getGraph(kwargs.get('graph'))

        lineageExistences = {}
        for (lineage, timePoint) in graph:
            if not lineage in lineageExistences:
                lineageExistences[lineage] = {'start':timePoint, 'end':timePoint }
            else:
                if lineageExistences[lineage]['start']>timePoint:
                    lineageExistences[lineage]['start']=timePoint
                if lineageExistences[lineage]['end']<timePoint:
                    lineageExistences[lineage]['end']=timePoint
        return lineageExistences

    def getTrellisCoordinates(self, **kwargs):
        trellisOrder = self.getTrellisOrder(kwargs.get('trellisOrder'))
        lineageSizes=self.getLineageSizes(kwargs.get('lineageSizes'))
        verticalSpacing = self.getVerticalSpacing(kwargs.get('verticalSpacing'))
        currentY = 0
        trellisCoordinates = {}
        for lineage in trellisOrder:
            currentY+= lineageSizes[lineage]
            trellisCoordinates[lineage] = currentY
            currentY+= lineageSizes[lineage]+ verticalSpacing
        return trellisCoordinates

    def plotTrellis(self, ax, **kwargs):
        trellisOrder = self.getTrellisOrder(kwargs.get('trellisOrder'))
        lineageSizes=self.getLineageSizes(kwargs.get('lineageSizes'))
        lineageExistences = self.getLineageExistences(lineageExistences)
        for lineage in trellisOrder:
            y = trellisCoordinates[lineage]
            xStart = lineageExistences[lineage]['start']*xDistance
            xEnd = lineageExistences[lineage]['end']*xDistance
            radius = lineageSizes[lineage]
            ax.plot([xStart,xEnd],[y, y], 'k' ,lw=radius/2)

    def getNodeCoordinates(self, **kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        trellisCoordinates = kwargs.get('trellisCoordinates')
        if trellisCoordinates is None:
            trellisCoordinates = self.getTrellisCoordinates(**kwargs)
        xDistance = self.getXDistance(kwargs.get('xDistance'))

        coordinates = {}
        for node in graph:
            (lineage, timePoint) = node
            y = trellisCoordinates[lineage]
            x = timePoint*xDistance
            coordinates[node] = (x,y)
        return coordinates




    def getEdgeCoordinates(self, **kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        radiusAttribute= self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        nodeCoordinates= self.getNodeCoordinates(**kwargs)
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        edgeCoordinates= {}
        for source,sink,d in graph.edges(data=True):
            sourceLineage, sourceTimePoint = source
            sinkLineage, sinkTimePoint = sink
            x0,y0 = nodeCoordinates[source]
            x1,y1 =  nodeCoordinates[sink]
            edgeCoordinates[(source,sink)] = (x0,y0, x1, y1)
        return edgeCoordinates

    def plotEdges(self, ax, offsetCentroid=False, **kwargs):
        plotArgs = kwargs.get('plotArgs',{})
        for (source,sink),(x0,y0, x1, y1) in self.getEdgeCoordinates(**kwargs).items():
            if offsetCentroid:
                (x0,y0),(x1,y1) = self.getNodeConnectionCoordinate(source,sink,**kwargs)
            if kwargs.get('bezier'):
                xData = []
                yData = []
                for x,y in self.interpolateEdge(x0,y0, x1, y1, **kwargs):
                    xData.append(x)
                    yData.append(y)
                ax.plot(xData,yData,**plotArgs)
            else:
                ax.plot([x0,x1],[y0,y1],**plotArgs)

    def getSegmentWidthControlPoints(self, r0,r1, source, sink, **kwargs):
        #@todo: this ignores the distance between the points!
        radiusAttribute= self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        initRadius=self.getInitRadius(kwargs.get('initRadius'))
        r0,r1 = self.getEdgeRadii(source, sink, **kwargs)


        return [
                (r0,0),
                (r0+(r1-r0)*(0.1),0),
                (r0+(r1-r0)*(0.9),0)  ,
                (r1,0) ]

    def getNodeConnectionCoordinate(self, source, sink,**kwargs):

        graph = self.getGraph(kwargs.get('graph'))
        trellisOrder = self.getTrellisOrder(kwargs.get('trellisOrder'))
        lineageSizes=self.getLineageSizes(kwargs.get('lineageSizes'))
        verticalSpacing=self.getVerticalSpacing(kwargs.get('verticalSpacing'))
        trellisCoordinates = self.getTrellisCoordinates( **kwargs)
        radiusAttribute= self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        nodeCoordinates= self.getNodeCoordinates(**kwargs)
        initRadius=self.getInitRadius(kwargs.get('initRadius'))

        x0,y0 = nodeCoordinates[source]
        x1,y1 = nodeCoordinates[sink]

        r0 = self.getNodeRadius(source, **kwargs)
        r1 = self.getNodeRadius(sink, **kwargs)

        if source[0]!=sink[0]: #Decide if the clone comes out from the center or the side:
            if initRadius is None:
                raise ValueError('Initradius is required')
            r1 = initRadius
            if y0>y1:
                y0-= r0 - r1
            else:
                y0+= r0 - r1

        return((x0,y0),(x1,y1))

    def getSegmentOutline(self, source, sink, wavyness,**kwargs):
                  #Interpolate the bezier curve:

        graph = self.getGraph(kwargs.get('graph'))
        trellisOrder = self.getTrellisOrder(kwargs.get('trellisOrder'))
        lineageSizes=self.getLineageSizes(kwargs.get('lineageSizes'))
        verticalSpacing=self.getVerticalSpacing(kwargs.get('verticalSpacing'))
        trellisCoordinates = self.getTrellisCoordinates( **kwargs)
        radiusAttribute= self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        nodeCoordinates= self.getNodeCoordinates(**kwargs)
        initRadius=self.getInitRadius(kwargs.get('initRadius'))
        widthControlPoints = kwargs.get('widthControlPoints')
        stepCount = kwargs.get('stepCount')

        # Obtains start and end radius:
        r0,r1 = self.getEdgeRadii(source, sink, **kwargs)

        # The width of the curve is described by a beziercurve too (single dimension cubic)
        if widthControlPoints is None:
            widthControlPoints = self.getSegmentWidthControlPoints(r0,r1,source,sink,**kwargs)

        pathForward  =[]
        pathReverse = []

        #x0,y0 = nodeCoordinates[source]
        #x1,y1 = nodeCoordinates[sink]
        (x0,y0),(x1,y1) = self.getNodeConnectionCoordinate(source,sink,**kwargs)
        controlPoints = self.getEdgeInternalControlPoints(x0,y0, x1, y1, wavyness=wavyness)


        for bi,(bx,by) in enumerate(self.interpolateEdge(x0,y0, x1, y1,wavyness=wavyness,stepCount=stepCount)):
            t = (bi/(stepCount-1))
            angle = interpolateBezierAngle(controlPoints,t )
            currentWidth =  interpolateBezier(widthControlPoints,t=t)[0]  # r0 + (r1-r0)*t

            pathForward.append((bx+np.sin(angle)*currentWidth, by+np.cos(angle)*currentWidth))
            pathReverse.append((bx-np.sin(angle)*currentWidth,by-np.cos(angle)*currentWidth))
            #pathForward.append( f"{'M' if bi==0 else 'L'}{bx+np.sin(angle)*currentWidth} {by+np.cos(angle)*currentWidth}" )
            #pathReverse.append( f"L{bx-np.sin(angle)*currentWidth} {by-np.cos(angle)*currentWidth}")

        path= pathForward + list(reversed(pathReverse))
        return path

    def interpolateSegmentOuterEdge(self, source, sink, wavyness, timeAxis, **kwargs):

        graph = self.getGraph(kwargs.get('graph'))
        trellisOrder = self.getTrellisOrder(kwargs.get('trellisOrder'))
        lineageSizes=self.getLineageSizes(kwargs.get('lineageSizes'))
        verticalSpacing=self.getVerticalSpacing(kwargs.get('verticalSpacing'))
        trellisCoordinates = self.getTrellisCoordinates( **kwargs)
        radiusAttribute= self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        xDistance=self.getXDistance(kwargs.get('xDistance'))
        nodeCoordinates= self.getNodeCoordinates(**kwargs)
        initRadius=self.getInitRadius(kwargs.get('initRadius'))
        widthControlPoints=kwargs.get('widthControlPoints')
        kwargs['wavyness'] = wavyness
        # Obtains start and end radius:
        r0,r1 = self.getEdgeRadii(source, sink, **kwargs)

        # The width of the curve is described by a beziercurve too (single dimension cubic)
        if widthControlPoints is None:
            widthControlPoints = self.getSegmentWidthControlPoints(r0,r1,source,sink,**kwargs)

        pathForward  =[]
        pathReverse = []

        #x0,y0 = nodeCoordinates[source]
        #x1,y1 = nodeCoordinates[sink]
        (x0,y0),(x1,y1) = self.getNodeConnectionCoordinate(source,sink,**kwargs)
        controlPoints = self.getEdgeInternalControlPoints(x0,y0, x1, y1,**kwargs)

        timeStart ,timeEnd= source[1],sink[1]

        # Interpolate the center bezier curce
        xyCenter = np.array(
            [interpolateBezier(t=(t-timeStart)/(timeEnd-timeStart), points=controlPoints)
             for t in timeAxis])

        for timePoint in timeAxis:
            timeRatio = (timePoint-timeStart)/(timeEnd-timeStart) # value between 0 and 1
            bx,by = interpolateBezier(t=timeRatio, points=controlPoints)
            angle = interpolateBezierAngle(controlPoints,timeRatio )
            currentWidth =  interpolateBezier(widthControlPoints,t=timeRatio)[0]
            pathForward.append((bx+np.sin(angle)*currentWidth, by+np.cos(angle)*currentWidth))
            pathReverse.append((bx-np.sin(angle)*currentWidth,by-np.cos(angle)*currentWidth))


        return pathForward, pathReverse

    """
    This function maps x[y..y] values  (timeAxis, yStart,yEnd) onto a sub segment between source and sink
    yStart and yEnd should be iterables of floats between 0 and 1
    """
    def interpolateSubSegment(self, source, sink,  timeAxis, yStarts, yEnds, **kwargs):

        graph = self.getGraph(kwargs.get('graph'))
        wavyness = kwargs.get('wavyness')
        trellisOrder = self.getTrellisOrder(kwargs.get('trellisOrder'))
        lineageSizes=self.getLineageSizes(kwargs.get('lineageSizes'))
        verticalSpacing=self.getVerticalSpacing(kwargs.get('verticalSpacing'))
        trellisCoordinates = self.getTrellisCoordinates( **kwargs)
        radiusAttribute= self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        xDistance=self.getXDistance(kwargs.get('xDistance'))
        nodeCoordinates= self.getNodeCoordinates(**kwargs)
        initRadius=self.getInitRadius(kwargs.get('initRadius'))
        widthControlPoints=kwargs.get('widthControlPoints')

        # Obtains start and end radius:
        r0,r1 = self.getEdgeRadii(source, sink, **kwargs)

        # The width of the curve is described by a beziercurve too (single dimension cubic)
        if widthControlPoints is None:
            widthControlPoints = self.getSegmentWidthControlPoints(r0,r1,source,sink,**kwargs)

        pathForward  =[]
        pathReverse = []
        pathCenter=[]
        pathWidth=[]
        segmentWidthRatio=[]
        segmentWidthAbsolute=[]
        timeRatios=[]

        #x0,y0 = nodeCoordinates[source]
        #x1,y1 = nodeCoordinates[sink]
        (x0,y0),(x1,y1) = self.getNodeConnectionCoordinate(source,sink,**kwargs)
        controlPoints = self.getEdgeInternalControlPoints(x0,y0, x1, y1, wavyness=wavyness)

        timeStart ,timeEnd= source[1],sink[1]

        # Interpolate the center bezier curce
        xyCenter = np.array(
            [interpolateBezier(t=(t-timeStart)/(timeEnd-timeStart), points=controlPoints)
             for t in timeAxis])


        for timePoint,yStart,yEnd in zip(timeAxis,yStarts, yEnds):
            #if yStart<0 or yEnd<0 or yStart>1 or yEnd>1:
            #    raise ValueError('yStarts and yEnds should be values between 0 and 1' )

            timeRatio = (timePoint-timeStart)/(timeEnd-timeStart) # value between 0 and 1, position on the curve
            if timeRatio<0 or timeRatio>1:
                continue

            timeRatios.append(timeRatio)
            bx,by = interpolateBezier(t=timeRatio, points=controlPoints)
            angle = interpolateBezierAngle(controlPoints,timeRatio )
            currentRadius =  interpolateBezier(widthControlPoints,t=timeRatio)[0]

            absoluteYStartRadius = (yStart*currentRadius*2)-currentRadius
            absoluteYEndRadius = (yEnd*currentRadius*2)-currentRadius

            pathForward.append((bx+np.sin(angle)*absoluteYStartRadius, by+np.cos(angle)*absoluteYStartRadius))
            pathReverse.append((bx+np.sin(angle)*absoluteYEndRadius,by+np.cos(angle)*absoluteYEndRadius))

            pathCenter.append((bx,by))
            pathWidth.append(currentRadius)
            segmentWidthAbsolute.append(absoluteYStartRadius-absoluteYEndRadius)
            segmentWidthRatio.append(yEnd-yStart)

        return pathForward, pathReverse, pathWidth,segmentWidthAbsolute, segmentWidthRatio,timeRatios

    def getEdgeInternalControlPoints(self,x0orSource,y0orSink, x1=None, y1=None, **kwargs):

        wavyness = kwargs.get('wavyness')
        if wavyness is None:
            raise ValueError('Supply wavyness')
        trellisOrder = self.getTrellisOrder(kwargs.get('trellisOrder'))
        graph=self.getGraph(kwargs.get('graph'))
        lineageSizes=self.getLineageSizes(kwargs.get('lineageSizes'))
        verticalSpacing=self.getVerticalSpacing(kwargs.get('verticalSpacing'))
        trellisCoordinates = self.getTrellisCoordinates(**kwargs)
        radiusAttribute= self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius = self.getDefaultRadius(kwargs.get('defaultRadius'))
        xDistance=self.getXDistance(kwargs.get('xDistance'))
        nodeCoordinates= self.getNodeCoordinates(**kwargs)
        initRadius=self.getInitRadius(kwargs.get('initRadius'))

        if x1 is None or y1 is None:
            #nodeCoordinates = self.getNodeCoordinates()[(x0orSource,y0orSink)] # @todo params
            #x0,y0 = nodeCoordinates[source]
            #x1,y1 = nodeCoordinates[sink]
            source,sink = x0orSource,y0orSink
            (x0,y0),(x1,y1) = self.getNodeConnectionCoordinate(source,sink, **kwargs)
        else:
            x0  = x0orSource
            y0  = y0orSink

        controlx0 = x0+(x1-x0)*wavyness
        controly0 = y0
        controlx1 = x0+(x1-x0)*(1-wavyness)
        controly1 = y1
        controlx0A=controlx0
        controlx1A=controlx1

        return [
            (x0,y0),
            (controlx0, controly0),
            (controlx1, controly1),
            (x1, y1)
            ]

    def interpolateEdge(self,x0,y0, x1, y1, wavyness=None, stepCount=30,**kwargs):
        kwargs['wavyness'] = wavyness
        kwargs['stepCount'] = stepCount

        if wavyness is None or stepCount is None:
            raise ValueError('Wavyness and stepCount are required')

        controlPoints = self.getEdgeInternalControlPoints(x0,y0, x1, y1, **kwargs)
        interpolatedBezier = interpolateBezier( points=controlPoints, steps=stepCount)
        return [
            (bx,by)
            for bi,(bx,by) in enumerate(interpolatedBezier)
        ]


    def setRadiusAttributes(self,**kwargs ):
        radiusAttribute = self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier = self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius=kwargs.get('defaultRadius')
        graph = self.getGraph(kwargs.get('graph'))
        if defaultRadius is not None:
            defaultRadius =defaultRadius/radiusMultiplier
        else:
            defaultRadius=self.getDefaultRadius()
        self.defaultRadius = defaultRadius
        self.radiusAttribute = radiusAttribute
        self.radiusMultiplier = radiusMultiplier
        self.lineageSizes = self.getSizePerLineage(**kwargs)


    def optimizeGraph(self,graph,
                  trellisOrder, radiusAttribute = 'radius', xDistance =100,
                  verticalSpacing=30, defaultRadius=1, radiusMultiplier=1, initRadius=0.1, minRadius=0.1,plot=True, ax=None, quitAtScore=None,
                      storeParams=False):


        kwargs = {
            'trellisOrder':trellisOrder,
            'graph':graph,
            'radiusAttribute':radiusAttribute,
            'xDistance':xDistance,
            'verticalSpacing':verticalSpacing,
            'defaultRadius':defaultRadius,
            'radiusMultiplier':radiusMultiplier,
            'initRadius':initRadius,
            'minRadius':minRadius,
            'plot':plot,
            'ax':ax,
            'quitAtScore':quitAtScore,
            'storeParams':storeParams
        }

        if defaultRadius is not None:
            defaultRadius =defaultRadius/radiusMultiplier
            kwargs['defaultRadius'] = defaultRadius
        if storeParams:
            self.graph = graph
            self.trellisOrder = trellisOrder
            self.setRadiusAttributes(graph=graph,
                                     defaultRadius=defaultRadius,
                                     radiusAttribute=radiusAttribute,
                                     radiusMultiplier=radiusMultiplier )
            self.xDistance = xDistance
            self.verticalSpacing = verticalSpacing
            self.minRadius=minRadius
            self.initRadius=initRadius


        overlaps = 0
        # Calculate maximum size per lineage:
        lineageSizes = self.getSizePerLineage(
            graph=graph, radiusAttribute=radiusAttribute,
            defaultRadius=defaultRadius,
            radiusMultiplier=radiusMultiplier,
            minRadius=minRadius)
        kwargs['lineageSizes'] = lineageSizes

        # Calculate lineage existences: (At what timepoint does what clone exist)
        lineageExistences = self.getLineageExistences(graph=graph)
        kwargs['lineageExistences'] = lineageExistences

        if storeParams:
            self.lineageSizes = lineageSizes
            self.lineageExistences = lineageExistences

        # Calculate trellis y coordinates:
        trellisCoordinates = self.getTrellisCoordinates(**kwargs)
        kwargs['trellisCoordinates'] = trellisCoordinates
        nodeCoordinates = self.getNodeCoordinates(**kwargs)
        kwargs['nodeCoordinates'] = nodeCoordinates
        # Display the edges:
        # This set of edges should not overlap:
        checkOverlap = []
        edges=[]
        lengths = 0
        for source,sink,d in graph.edges(data=True):
            sourceLineage, sourceTimePoint = source
            sinkLineage, sinkTimePoint = sink
            x0,y0 = nodeCoordinates[source]
            x1,y1 =  nodeCoordinates[sink]
            (x0,y0),(x1,y1) = self.getNodeConnectionCoordinate(source,sink,**kwargs)
            radius0 = graph[source].get(radiusAttribute,defaultRadius)*radiusMultiplier
            radius1 = graph[sink].get(radiusAttribute,defaultRadius)*radiusMultiplier
            lengths += np.sqrt( np.power(y1-y0,2) + np.power(x1-x0,2) )
            p1,q1 = expandStart(x0,y0,x1,y1 )
            checkOverlap.append(((p1, q1), source, sink))
            edges.append((x0,y0, x1, y1))
            if quitAtScore is not None and overlaps+lengths*0.0001>quitAtScore:
                return None


        edgeIntersections = {}
        intersectingClones = set()
        for edgeIndexA,(lineA, sourceA, sinkA) in enumerate(checkOverlap):
            for edgeIndexB,(lineB, sourceB, sinkB) in enumerate(checkOverlap):
                if lineA!=lineB and sourceA[0]!=sourceB[0]:
                    # Check overlap
                    if do_intersect(*(lineA+lineB)):
                        edgeIntersections[edgeIndexA] = True
                        edgeIntersections[edgeIndexB] = True
                        intersectingClones.add(sourceA[0])
                        intersectingClones.add(sourceB[0])
                        intersectingClones.add(sinkA[0])
                        intersectingClones.add(sinkB[0])
                        overlaps+=1
                        if quitAtScore is not None and overlaps+lengths*0.0001>quitAtScore:
                            return None

        lines=[]
        for edgeIndex, edge in enumerate(edges):
            radius=1
            color='r' if edgeIndex in edgeIntersections else 'g'
            #lines.append([*edge, color, (radius1+radius0)/2])
            lines.append(
                [*(checkOverlap[edgeIndex][0][0]+checkOverlap[edgeIndex][0][1]),
                 color,
                 (radius1+radius0)/2])

        if plot:
            if ax is None:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                fig.canvas.draw()
            for i, (x1,y1,x2,y2,color,radius) in enumerate(lines):
                ax.plot( [x1,x2],[y1,y2],color, lw=radius)
            for (lineage, timePoint) in graph:
                y = trellisCoordinates[lineage]
                x = timePoint*xDistance
                radius = graph.nodes[(lineage, timePoint)].get(radiusAttribute, defaultRadius)*radiusMultiplier
                ax.add_patch(plt.Circle((x, y), radius, color='r', alpha=0.5))
                plt.plot((x,x),(y,y+radius), color='k' )

        return overlaps+lengths*0.0001, intersectingClones

    def getXCoordinate(self, x):
        return x*self.xDistance

    def drawNodeRadii(self, ax,**kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        trellisCoordinates = self.getTrellisCoordinates(kwargs.get('trellisCoordinates'))
        radiusAttribute = self.getRadiusAttribute(kwargs.get('radiusAttribute'))
        radiusMultiplier=self.getRadiusMultiplier(kwargs.get('radiusMultiplier'))
        defaultRadius=self.getDefaultRadius(kwargs.get('defaultRadius'))
        xDistance = self.getXDistance(kwargs.get('xDistance'))

        for (lineage, timePoint) in graph:
            y = trellisCoordinates[lineage]
            x = timePoint*xDistance

            radius =  self.getNodeRadius(
                                         (lineage, timePoint),
                                            **kwargs)

            ax.add_patch(plt.Circle((x, y), radius, color='r', alpha=0.5))

    # Returns true if  node is a leaf
    def nodeIsLeaf(self, node, **kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        isLeaf = self.graph.nodes[node].get('leaf')
        return isLeaf

    # Returns True if the clone died out after this timepoint
    def nodeIsExtinct(self, node, **kwargs):
        graph = self.getGraph(kwargs.get('graph'))
        isExtinct = self.graph.nodes[node].get('extinct')
        return isExtinct

    def annotateNodes(self, ax, x_offset=0.5, nodesToAnnotate=None, **kwargs):
        plotArgs = kwargs.get('plotArgs',{})
        coordinates = self.getNodeCoordinates(**kwargs)
        for node,(x,y) in coordinates.items():
            if nodesToAnnotate is not None and not node in nodesToAnnotate:
                continue
            label = None
            if self.nodeIsExtinct(node):
                label = f'{node[0]} †'
            elif self.nodeIsLeaf(node):
                label = f'{node[0]}'
            if label is not None:
                ax.text( x+x_offset,
                        y,
                        label,
                        verticalalignment='center',
                        horizontalalignment='left',
                        **plotArgs)

    """
    Plot patches (edges of the graph) over the tree
        ax : axis to plot to
        facecolor : default color to fill the patches with,
                    None to not plot if no data is supplied in patchData
        stepCount : precision
        linewidth : width of edges around the patches
        edgecolor : edge color
        zorder : default z-order
        wavyness
        patchData : dataframe like:
                    pd.DataFrame(
                    {
                        ((2, 10), (3, 22)):{'color':'red'},
                        ((2, 10), (5, 22)):{'color':'red'}

                    }).T

    """
    def plotPatches(self, ax,
        facecolor=(0.7,0.7,0.7,1),
        stepCount=30,
        linewidth=0.0, edgecolor='b',wavyness=0.4,linestyle='-',
        zorder=2,
        patchData=None ):
        for source, sink in self.graph.edges():
            fc = facecolor
            ec = edgecolor
            l = linewidth
            ls = linestyle
            z = zorder
            assignedFromPatchData = False
            if patchData is not None and (source,sink) in list(patchData.index):
                # Find its postion @slow
                assignedFromPatchData = True
                index = list(patchData.index).index((source,sink))
                try:
                    fc =  patchData.iloc[index]['color']
                    if np.isnan(fc):
                        fc = facecolor
                except Exception as e:
                    pass
                try:
                    z =  patchData.iloc[index]['zorder']
                    if np.isnan(z):
                        z = zorder
                except Exception as e:
                    pass

                try:
                    ec =  patchData.iloc[index]['edgecolor']
                    if np.isnan(ec):
                        ec = edgecolor
                except Exception as e:
                    pass
                try:
                    l =  patchData.iloc[index]['linewidth']
                    if np.isnan(l):
                        l = linewidth
                except Exception as e:
                    pass
                try:
                    ls =  patchData.iloc[index]['linestyle']
                    if np.isnan(ls):
                        ls = linestyle
                except Exception as e:
                    pass
            if assignedFromPatchData or ( not assignedFromPatchData and facecolor is not None):
                ax.add_patch(
                    plt.Polygon(
                        self.getSegmentOutline( source, sink , wavyness=wavyness,stepCount=stepCount ) ,
                               facecolor=fc, linestyle=ls, lw=l, edgecolor=ec, zorder=z  ) )


    def getEmptyPlot(self):
        fig, ax = plt.subplots()
        return fig, ax


    def plot_vertical_lines(self, ax, coords, linewidth=1, linestyle=':', **kwargs):
        for coord in coords:
            ax.axvline(self.getXCoordinate(coord), lw=linewidth, linestyle=linestyle,**kwargs)


    def plot_xticks(self, ax, coords, **kwargs):
        xticks = np.insert(coords, 0, 0)
        ax.set_xticks(self.getXCoordinate(xticks), **kwargs)


    def plot_stack(self, cellCnvBarcode, tf, ax, fig):
        """
        Plot stack of subclones on top of main clones
        Args:
            cellCnvBarcode(dict) : {cell_name : clone_id }
            tf(pd.DataFrame) : dataframge containing abundance of subclones over time
            ax : axis to plot onto
            fig : handle to figure
        Returns:
            cloneColors (dict) : colors/markers assigned to every subclone
        """

        random.seed(3)
        graph=self.graph
        colorsList = list(plt.get_cmap('tab20').colors)
        markersList = 'os^<>'
        #sns.palplot(colorsList)
        del colorsList[-6:-4] # remove grays
        #sns.palplot(colorsList)
        colorsList = colorsList[::2] + colorsList[1::2]
        random.shuffle(colorsList)
        colors =  list(itertools.product(markersList,colorsList))
        random.shuffle(colors)
        color_index = 0
        wavyness=0.5
        self.xDistance=10
        self.verticalSpacing = 0.2
        #lg.initRadius=0.5
        #lg.setRadiusAttributes(defaultRadius=0.5)

        cellJitter=0.75*self.xDistance

        ax.clear()


        for source, sink in graph.edges():
            ax.add_patch( plt.Polygon( self.getSegmentOutline( source, sink , wavyness=wavyness,stepCount=30 ) ,
                                      facecolor=(0.7,0.7,0.7,1), lw=0.0, edgecolor='b'  ) )

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_yaxis().set_visible(False)
        # This plots the internals of the curves
        self.plotEdges(ax, bezier=True,wavyness=wavyness,stepCount=30,plotArgs={'lw':0}, offsetCentroid=True)


        fig.canvas.draw()
        #raise KeyboardInterrupt()

        #lg.drawNodeRadii(ax)

        cloneColors = {}
        #Iterate first edges where sourceLineage!=sinkLineage
        for source, sink in [
            (u,v) for u,v in graph.edges()
                if u[0]!=v[0]]+[
            (u,v) for u,v in graph.edges() if u[0]==v[0]]:

            timeStart,timeEnd= source[1],sink[1]
            internalControl = self.getEdgeInternalControlPoints(source,sink,wavyness=wavyness)
            if not sink[0] in tf.columns:
                continue

            timeAxis = [tp for tp in tf[sink[0]].columns if tp>=timeStart and tp<=timeEnd]
            if len(timeAxis)==0:
                continue

            #Obtain the abundance of sub copy number clones:
            q = tf[sink[0]].copy()
            #Interpolate the start and end value:
            if not timeStart in q.columns:
                q[timeStart] = np.nan
            if not timeEnd in q.columns:
                q[timeEnd] = np.nan


            # Resort the columns so they are in order with time (Otherwise the interpolation fails)
            q = q.reindex(sorted(q.columns), axis=1)
            q =q.T.interpolate('index').T

            timeAxis = list(sorted(list(set(timeAxis+[timeStart,timeEnd]))))

            currentY=None
            cloneAbundance = q.loc[:,timeAxis]
            cloneAbundance = cloneAbundance[cloneAbundance.sum(1)>0]


            yCenter = np.array([
                interpolateBezier(t=(t-timeStart)/(timeEnd-timeStart), points=internalControl)[1]
                for t in timeAxis])


            cellsAtThisTimepoint = cellCnvBarcode.loc[sink[1]]
            cellsAtThisTimepointAndCluster = cellsAtThisTimepoint[cellsAtThisTimepoint['cluster']==sink[0]]

            prevColors = collections.deque(maxlen=4)

            cloneOrder = [clone for clone,row in cloneAbundance.iterrows()]
            for cloneIndex,(clone,row) in enumerate(cloneAbundance.iterrows()):

                if currentY is None:
                    currentY=np.zeros(row.values.shape[0])
                    # Prevent a mixture of missing and non-missing data, initialise nan if no data:
                    currentY[np.isnan(row.values)] = np.nan
                nextY=currentY+row.values
                #print(currentY,nextY)
                r = self.interpolateSubSegment(source,sink, wavyness=wavyness,
                              timeAxis=timeAxis,
                              yStarts=currentY, yEnds=nextY)


                forward = r[0]
                reverse  = r[1]
                if len(cloneOrder)>1:
                    nextColor,nextMarker = cloneColors.get(cloneOrder[(cloneIndex+1)%len(cloneOrder)],(None,None))
                else:
                    nextColor,nextMarker= (None,None)
                if not clone in cloneColors:
                    marker,color = colors[color_index]
                    color_index = color_index+1 if color_index<len(colors)-1 else 0

                    while color in prevColors or color == nextColor:
                        marker,color =  colors[color_index]
                        color_index = color_index+1 if color_index<len(colors)-1 else 0
                    cloneColors[clone] = (color, marker)
                else:
                    color, marker = cloneColors[clone]

                prevColors.append(color)
                ax.add_patch(
                        plt.Polygon(
                           forward+list(reversed(reverse)),facecolor=color #,lw=0.1,  edgecolor='#000000'
                    )
                )
                if True:
                    # Add cells:
                    cellX =[]
                    cellY =[]
                    cellLabels=[]
                    cellsWithCurrentBarcode = cellsAtThisTimepointAndCluster[cellsAtThisTimepointAndCluster['barcode']==clone]
                    for cell,_ in cellsWithCurrentBarcode.iterrows():
                        #print(cell, clone)
                        xForward = np.array([x for x,y in forward])
                        yForward = np.array([y for x,y in forward])
                        xReverse = np.array([x for x,y in reverse])
                        yReverse = np.array([y for x,y in reverse])

                        isLeaf = (self.graph.out_degree(sink)==0)
                        nx,ny = self.getNodeCoordinates()[sink]
                        if isLeaf:
                            low = -cellJitter*1.5
                            high = 0
                        else:
                            low = -cellJitter
                            high = cellJitter
                        x = np.random.uniform(low=low,high=high) + nx

                        if True:
                            yMin = np.interp(x,xForward[np.isfinite(xForward)],yForward[np.isfinite(yForward)])
                            yMax = np.interp(x,xReverse[np.isfinite(xReverse)],yReverse[np.isfinite(yReverse)])
                        else:
                            yMin = xForward[1]
                            yMax = xReverse[1]


                        y = np.random.uniform(
                            low=yMin,
                            high=yMax)

                        #x.add_patch(plt.Circle((x, y), 1, color=color, alpha=1))
                        cellX.append(x)
                        cellY.append(y)
                        cellLabels.append(cell)
                    if len(cellX)>0:
                        size = 14
                        plt.scatter(cellX,cellY, s=size*2, marker=marker,color='k', zorder=6,alpha=0.1)
                        plt.scatter(cellX,cellY, s=size*1.5,marker=marker, color='k', zorder=7,alpha=0.5)
                        plt.scatter(cellX,cellY, s=14, color=color, marker=marker, zorder=10, edgecolor='w',linewidth=0.3)
                currentY = nextY


        fig.canvas.draw()
        return cloneColors


    def plot_stack_X(self,
        cell_to_clone_mapping,
        clone_abundance_frame,
        ax,
        wavyness=0.8,
        cellJitter=10,
        ):


        random.seed(4)

        # Assign clones to unique marker/color combinations,
        #this is only possible _given_ a ordering
        colorsList = list(plt.get_cmap('tab20').colors)
        markersList = 'os^<>'
        del colorsList[-6:-4] # remove grays
        random.shuffle(colorsList)
        colors =  itertools.product(markersList,colorsList)

        cloneColors = {}

        graph=self.graph
        tf = clone_abundance_frame
        cellCnvBarcode = cell_to_clone_mapping

        #Iterate first edges where sourceLineage!=sinkLineage
        for source, sink in [
            (u,v) for u,v in graph.edges()
                if u[0]!=v[0]]+[
            (u,v) for u,v in graph.edges() if u[0]==v[0]]:

            timeStart,timeEnd= source[1],sink[1]

            internalControl = self.getEdgeInternalControlPoints(source,sink,wavyness=wavyness)

            #internalControl = (x0,y0),(x1,y1) = self.getNodeConnectionCoordinate(source,sink)
            #controlPoints = self.getEdgeInternalControlPoints(x0,y0, x1, y1, wavyness=wavyness)
            if not sink[0] in tf.columns:
                continue

            timeAxis = [tp for tp in tf[sink[0]].columns if tp>=timeStart and tp<=timeEnd]
            if len(timeAxis)==0:
                continue

            #Obtain the abundance of sub copy number clones:
            q = tf[sink[0]].copy()
            #Interpolate the start and end value:
            if not timeStart in q.columns:
                q[timeStart] = np.nan
            if not timeEnd in q.columns:
                q[timeEnd] = np.nan


            # Resort the columns so they are in order with time (Otherwise the interpolation fails)
            q = q.reindex(sorted(q.columns), axis=1)
            q = q.T.interpolate('index').T

            timeAxis = list(sorted(list(set(timeAxis+[timeStart,timeEnd]))))

            currentY=None
            cloneAbundance = q.loc[:,timeAxis]
            cloneAbundance = cloneAbundance[cloneAbundance.sum(1)>0]


            yCenter = np.array([
                interpolateBezier(t=(t-timeStart)/(timeEnd-timeStart), points=internalControl)[1]
                for t in timeAxis])


            cellsAtThisTimepoint = cellCnvBarcode.loc[sink[1]]
            cellsAtThisTimepointAndCluster = cellsAtThisTimepoint[cellsAtThisTimepoint['cluster']==sink[0]]

            prevColors = collections.deque(maxlen=4)

            cloneOrder = [clone for clone,row in cloneAbundance.iterrows()]
            for cloneIndex,(clone,row) in enumerate(cloneAbundance.iterrows()):

                if currentY is None:
                    currentY=np.zeros(row.values.shape[0])
                    # Prevent a mixture of missing and non-missing data, initialise nan if no data:
                    currentY[np.isnan(row.values)] = np.nan
                nextY=currentY+row.values
                #print(currentY,nextY)
                r = self.interpolateSubSegment(source,sink, wavyness=wavyness,
                              timeAxis=timeAxis,
                              yStarts=currentY, yEnds=nextY)


                forward = r[0]
                reverse  = r[1]
                if len(cloneOrder)>1:
                    nextColor,nextMarker = cloneColors.get(cloneOrder[(cloneIndex+1)%len(cloneOrder)],(None,None))
                else:
                    nextColor,nextMarker= (None,None)
                if not clone in cloneColors:

                    marker,color = colors[color_index]
                    color_index = color_index+1 if color_index<len(colors)-1 else 0
                    while color in prevColors or color == nextColor:
                        marker,color = colors[color_index]
                        color_index = color_index+1 if color_index<len(colors)-1 else 0
                    cloneColors[clone] = (color, marker)
                else:
                    color, marker = cloneColors[clone]

                prevColors.append(color)
                ax.add_patch(
                        plt.Polygon(
                           forward+list(reversed(reverse)),facecolor=color ,lw=1,  edgecolor='#000000'
                    )
                )

                # Add cells:
                cellX =[]
                cellY =[]
                cellLabels=[]
                cellsWithCurrentBarcode = cellsAtThisTimepointAndCluster[cellsAtThisTimepointAndCluster['barcode']==clone]
                for cell in cellsWithCurrentBarcode:
                    xForward = np.array([x for x,y in forward])
                    yForward = np.array([y for x,y in forward])
                    xReverse = np.array([x for x,y in reverse])
                    yReverse = np.array([y for x,y in reverse])

                    x = np.random.uniform(
                        low=np.nanmin(xForward)*(1+cellJitter),
                        high=np.nanmax(xForward)*(1-cellJitter))

                    high = np.nanmin((np.nanmax(xForward), (np.nanmean(xForward)+(cellJitter))))
                    low = np.nanmin((np.nanmax(xForward), (np.nanmean(xForward)-(cellJitter))))


                    x = np.random.uniform(
                        low=low,
                        high=high)


                    if True:
                        yMin = np.interp(x,xForward[np.isfinite(xForward)],yForward[np.isfinite(yForward)])
                        yMax = np.interp(x,xReverse[np.isfinite(xReverse)],yReverse[np.isfinite(yReverse)])
                    else:
                        yMin = xForward[1]
                        yMax = xReverse[1]


                    y = np.random.uniform(
                        low=yMin,
                        high=yMax)

                    #x.add_patch(plt.Circle((x, y), 1, color=color, alpha=1))
                    cellX.append(x)
                    cellY.append(y)
                    cellLabels.append(cell)
                if len(cellX)>0:
                    size = 14
                    plt.scatter(cellX,cellY, s=size*2, marker=marker,color='k', zorder=6,alpha=0.1)
                    plt.scatter(cellX,cellY, s=size*1.5,marker=marker, color='k', zorder=7,alpha=0.5)
                    plt.scatter(cellX,cellY, s=14, color=color, marker=marker, zorder=10, edgecolor='w',linewidth=0.3)
            currentY = nextY


    def plotSingleCells(
            self,
            cellData,
            ax = None,
            fig =None,
            cloneAttribute = 'cluster',
            timeAttribute = 'tp',
            colorAttribute = 'color',
            labelAttribute = 'label',
            markerAttribute= 'marker',
            sizeAttribute = 'size',
            edgeWidthAttribute='edge_width',
            edgeColorAttribute='edge_color',
            wavyness = 0.4,
            defaultCellColor = 'k',
            defaultCellSize = 20,
            cellJitter=10,
            plotPatches=False,
            defaultCellMarker = 'o',
            enableShadow=True,
            shadow_x_offset = 0.1,
            shadow_y_offset = -0.1,
            shadow_alpha = 0.3,
            shadow_size=1.3,
            shadow_color = (0.1,0.1,0.1),
            defaultLinewidth = 1,
            defaultEdgecolor = (1,1,1),
            defaultZorder=0,
            zorderAttribute='z-order',
            xOffsetAttribute='x',
            yOffsetAttribute='y',
            defaultXOffset = 0,
            defaultYOffset = None, # Value between 0 and 1
            yJitterRatio = 1 # yJitter in respect to the size of the clone
            ):
        if ax is None:
            fig, ax = self.getEmptyPlot()


        for clone, tp in self.graph:
            node = (clone, tp)
            isLeaf = (self.graph.out_degree(node)==0)

            # Dataframe with all cells which are present at this timepoint and clone
            cellsToDraw = cellData[
                    (cellData[cloneAttribute]==clone) &
                    ( cellData[timeAttribute]==tp)
            ]
            radius = self.getNodeRadius(node=node)

            x = []
            y = []
            color = []
            marker = []
            size = []
            annotation = []
            edge_color = []
            edge_width = []
            zOrders=[]
            for cellIndex, metaData in cellsToDraw.iterrows():

                # Cell position calculation
                nx,ny = self.getNodeCoordinates()[node]
                if isLeaf:
                    low = -cellJitter*1.5
                    high = 0
                else:
                    low = -cellJitter
                    high = cellJitter

                # Cell X offset:
                try:
                    x_offset = metaData[xOffsetAttribute]
                except Exception as e:
                    x_offset = defaultXOffset

                if x_offset is not None:
                    x_offset*=self.xDistance
                else:
                    x_offset=0

                try:
                    # y offset should be a value between 0 and 1
                    y_offset = metaData[yOffsetAttribute]
                except Exception as e:
                    y_offset = defaultYOffset

                if y_offset is not None:
                    if y_offset<0 or y_offset>1:
                        raise ValueError(f'{yOffsetAttribute} (y_offset) should be between 0 and 1')
                    cell_y = -radius + radius*2*y_offset + ny
                else:
                    cell_y = np.random.uniform(low=-radius*yJitterRatio,high=radius*yJitterRatio) + ny

                cell_x = np.random.uniform(low=low,high=high) + nx + x_offset
                x.append(cell_x)
                y.append(cell_y)

                # Cell color:
                try:
                    c = metaData[colorAttribute]
                except Exception as e:
                    c = defaultCellColor
                color.append(c)

                # Cell Marker:
                try:
                    m = metaData[markerAttribute]
                except Exception as e:
                    m = defaultCellMarker
                marker.append(m)

                # Cell marker size
                try:
                    s = metaData[sizeAttribute]
                except Exception as e:
                    s = defaultCellSize
                size.append(s)

                # Edge color
                try:
                    edgecolor = metaData[edgeColorAttribute]
                except Exception as e:
                    edgecolor = defaultEdgecolor
                edge_color.append(edgecolor)

                # Edge width
                try:
                    linewidth = metaData[edgeWidthAttribute]
                except Exception as e:
                    linewidth = defaultLinewidth
                edge_width.append(linewidth)

                # Z-order
                try:
                    zorder = metaData[zorderAttribute]
                except Exception as e:
                    zorder = defaultZorder
                zOrders.append(zorder)

                ax.scatter([cell_x],[cell_y],c=[c],s=[s],
                    marker=m,
                    zorder=8+ zorder,
                    linewidth=linewidth,
                    edgecolor=edgecolor
                    )
                if enableShadow:
                    ax.scatter([cell_x+shadow_x_offset],
                                [cell_y+shadow_y_offset],
                                c=[shadow_color],
                                s=[ (s+linewidth)*shadow_size],
                                alpha= shadow_alpha, marker=m, zorder=7)


        if plotPatches:
            self.plotEdges(ax, bezier=True,wavyness=wavyness,stepCount=30,plotArgs={'linewidth':0}, offsetCentroid=True)
            self.plotPatches(ax=ax,wavyness=wavyness)

        # Remove plot spines:
        despine(ax)
        # Scale labels to plot size:
        format_x_axis_labels(ax)
        fig.canvas.draw()
